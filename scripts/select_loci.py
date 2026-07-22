"""Fresh RT-atlas search + control-locus selection (FDZ004.2).

Reads all_rts_metadata.parquet + all_meta_ncbi_gff.parquet from ERDA (read-only) and:
  1. counts transcriptomic vs metatranscriptomic RTs, split short/long by platform;
  2. per registry subclass (rt_subclass_registry), orders repr_c90 clusters by the longest-contig
     representative and picks the top-N reps whose RT gene has >= ORF_MIN ORFs upstream AND downstream;
  3. emits control_loci.tsv (consolidated) + selection_report.md (counts + per-subclass yield + ORF dist).

ERDA base URL from $ERDA_BASE (never written to any output). Usage:
  ERDA_BASE=$(cat ~/.erda_base) select_loci.py --out <dir> [--orf-min 3] [--n-short 3] [--n-long 1] [--pool 25]
"""
import os
import sys
import csv
import argparse
import urllib.request
import urllib.error
from collections import defaultdict
from pathlib import Path
import duckdb

# registry subclass -> (role, common_name). Verbatim rt_subclass_registry labels.
PANEL = [
    ("UG2", "positive", "DRT2"),
    ("UG28", "positive", "DRT9"),
    ("UG26", "positive", "UG26"),
    ("GII", "positive", "GII intron"),
    ("Retron Type VI-A", "positive", "retron VI-A"),
    ("Retron Type VI-B/UG11", "positive", "retron VI-B"),
    ("DGR", "positive", "DGR"),
    ("G2L4", "positive", "G2L4"),
    ("G2L cluster 3", "positive", "G2Lc3"),
    ("UG1", "negative", "DRT1"),
    ("UG15", "negative", "DRT4"),
    ("UG10", "negative", "DRT7"),
    ("CRISPR-associated", "negative", "RT-CRISPR"),
    ("UG16", "negative", "DRT5"),
]
LONG_PLATFORMS = ("PACBIO_SMRT", "OXFORD_NANOPORE")
READTYPE_SQL = (f"CASE WHEN platform IN {LONG_PLATFORMS} THEN 'long' ELSE 'short' END")
LONG_READ_BP = 500   # avg read length (bp) above which a SINGLE-end library is long — matches rnaseq_pipeline.LONG_READ_BP
ENA_URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run"
           "&fields=read_count,base_count,instrument_platform,instrument_model,library_layout&format=tsv")


def sanitize(subclass):
    return subclass.replace("Retron Type ", "retron_").replace("/", "_").replace(" ", "_")


def ena_read_stats(acc, timeout=20):
    """Average read length + library layout from ENA run metadata (dependency-free HTTP), or None on any failure.
    The atlas `platform` column labels merged/synthetic-long ILLUMINA (single-end, ~1 kb reads) as short; ENA's
    base_count/read_count exposes their true length so we can route them to minimap2. Cheap (~1 s/accession) and
    one-time at selection; the pipeline's map-time measured-length check stays the ground-truth override."""
    try:
        with urllib.request.urlopen(ENA_URL.format(acc=acc), timeout=timeout) as resp:
            rows = resp.read().decode().splitlines()
        if len(rows) < 2:
            return None
        f = rows[1].split("\t")   # run_accession, read_count, base_count, platform, model, library_layout
        rc, bc = int(f[1]), int(f[2])
        return {"avg_read_len": bc // rc if rc else None, "library_layout": f[5],
                "ena_platform": f[3], "ena_model": f[4]}
    except (urllib.error.URLError, ValueError, IndexError, OSError):
        return None


def confirm_read_type(platform_rt, stats):
    """Truth of read length from metadata, checked in this order:
      1. atlas `platform` says long (PacBio/ONT) -> long, no length test.
      2. ENA `instrument_platform` says long -> long, no length test (catches a native long-read run the atlas
         platform column mislabelled short — a short-read-*length* ONT/PacBio run would slip a length-only test).
      3. PAIRED -> short (paired-end is short-read Illumina, always; its per-spot avg can exceed the cutoff, so
         layout is checked BEFORE length — e.g. 2x300 bp computes to ~600 bp but is still short).
      4. SINGLE and avg read length >= LONG_READ_BP -> long (catches merged/synthetic-long ILLUMINA).
    Falls back to the platform label if ENA has no record; the pipeline's map-time measured median is the final
    override. Length is only ever tested for a single-end, non-long-platform library."""
    if platform_rt == "long":
        return "long"
    if stats:
        if (stats.get("ena_platform") or "").upper() in LONG_PLATFORMS:
            return "long"
        if (stats.get("library_layout") or "").upper() == "PAIRED":
            return "short"
        al = stats.get("avg_read_len")
        if al is not None and al >= LONG_READ_BP:
            return "long"
    return platform_rt


def orf_neighbors(genes, rt_id, rt_start, rt_end):
    """Return (#CDS upstream, #CDS downstream) of the RT gene on its contig, or None if RT not locatable."""
    cds = sorted((g for g in genes if g["type"] == "CDS"), key=lambda g: g["start"])
    if not cds:
        return None
    idx = next((i for i, g in enumerate(cds) if g["protein_id"] == rt_id), None)
    if idx is None:  # fall back to coordinate overlap
        idx = next((i for i, g in enumerate(cds) if g["start"] <= rt_end and g["end"] >= rt_start), None)
    if idx is None:
        return None
    return idx, len(cds) - 1 - idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--orf-min", type=int, default=3, help="strict tier-1 ORF cutoff (up AND down)")
    ap.add_argument("--orf-fallback", type=int, default=1, help="tier-2 fallback ORF cutoff for short buckets short of quota")
    ap.add_argument("--n-short", type=int, default=3)
    ap.add_argument("--n-long", type=int, default=1)
    ap.add_argument("--pool", type=int, default=25, help="candidate reps per (subclass,read_type) to ORF-check")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    base = os.environ["ERDA_BASE"].rstrip("/")
    M = f"read_parquet('{base}/all_rts_metadata.parquet')"
    G = f"read_parquet('{base}/all_meta_ncbi_gff.parquet')"
    con = duckdb.connect()
    con.execute("INSTALL httpfs; LOAD httpfs; SET enable_progress_bar=false;")

    subclasses = [p[0] for p in PANEL]
    inlist = ",".join("'" + s.replace("'", "''") + "'" for s in subclasses)

    # ---- 1a. counts ----
    def scalar(sql):
        return con.execute(sql).fetchone()[0]
    n_or = scalar(f"SELECT count(DISTINCT repr_c90) FROM {M} WHERE librarysource IN ('TRANSCRIPTOMIC','METATRANSCRIPTOMIC')")
    n_t = scalar(f"SELECT count(DISTINCT repr_c90) FROM {M} WHERE librarysource='TRANSCRIPTOMIC'")
    n_mt = scalar(f"SELECT count(DISTINCT repr_c90) FROM {M} WHERE librarysource='METATRANSCRIPTOMIC'")
    sl = con.execute(
        f"SELECT {READTYPE_SQL} rt, count(*) n_rows, count(DISTINCT repr_c90) n_clusters "
        f"FROM {M} WHERE librarysource='TRANSCRIPTOMIC' GROUP BY 1 ORDER BY 1").fetchall()
    sl = {r[0]: (r[1], r[2]) for r in sl}

    # ---- 1b. candidate pool: longest-contig rep per (subclass, read_type, repr_c90), top `pool` per bucket ----
    pool = con.execute(f"""
        WITH cand AS (
          SELECT id, accession, contig_id, contig_length, protein_start, protein_end, strand,
                 repr_c90, rt_subclass_registry AS subclass, {READTYPE_SQL} AS read_type
          FROM {M}
          WHERE rt_subclass_registry IN ({inlist}) AND librarysource='TRANSCRIPTOMIC'
                AND repr_c90 IS NOT NULL AND contig_id IS NOT NULL
        ),
        rep AS (
          SELECT *, row_number() OVER (PARTITION BY subclass, read_type, repr_c90
                                       ORDER BY contig_length DESC, id) rk FROM cand
        ),
        ranked AS (
          SELECT *, row_number() OVER (PARTITION BY subclass, read_type
                                       ORDER BY contig_length DESC, id) pool_rk
          FROM rep WHERE rk=1
        )
        SELECT subclass, read_type, pool_rk, id, accession, contig_id, contig_length,
               protein_start, protein_end, strand, repr_c90
        FROM ranked WHERE pool_rk <= {args.pool} ORDER BY subclass, read_type, pool_rk
    """).fetchall()
    cols = [d[0] for d in con.description]
    cand = [dict(zip(cols, r)) for r in pool]

    # ---- batched GFF fetch for all candidate contigs (one scan) ----
    contigs = sorted({c["contig_id"] for c in cand})
    genes_by_contig = defaultdict(list)
    if contigs:
        cin = ",".join("'" + c.replace("'", "''") + "'" for c in contigs)
        rows = con.execute(
            f'SELECT contig, "start", "end", strand, type, protein_id FROM {G} '
            f'WHERE contig IN ({cin})').fetchall()
        for contig, start, end, strand, typ, pid in rows:
            genes_by_contig[contig].append({"start": start, "end": end, "strand": strand,
                                            "type": typ, "protein_id": pid})

    # ---- ORF-gate + pick top-N per (subclass, read_type) ----
    picks = []                       # chosen loci rows
    orf_dist = defaultdict(list)     # subclass -> list of (up, down) over all candidates
    yield_ct = defaultdict(lambda: {"short": (0, 0), "long": (0, 0)})   # (n_strict, n_fallback)
    flags = []
    targets = {"short": args.n_short, "long": args.n_long}
    by_bucket = defaultdict(list)
    for c in cand:
        by_bucket[(c["subclass"], c["read_type"])].append(c)
    for sub, role, common in PANEL:
        for rt in ("short", "long"):
            scored = []              # (candidate, up, down) in contig-length-desc order
            for c in by_bucket.get((sub, rt), []):
                nb = orf_neighbors(genes_by_contig.get(c["contig_id"], []),
                                   c["id"], c["protein_start"], c["protein_end"])
                if nb is None:
                    continue
                orf_dist[sub].append(nb)
                scored.append((c, nb[0], nb[1]))
            target = targets[rt]
            chosen, used = [], set()
            for c, up, down in scored:                       # tier 1: strict orf_min/orf_min
                if len(chosen) >= target:
                    break
                if up >= args.orf_min and down >= args.orf_min:
                    chosen.append((c, up, down, f"{args.orf_min}/{args.orf_min}")); used.add(id(c))
            if rt == "short":                                # tier 2 (short only): relax to orf_fallback both sides
                for c, up, down in scored:
                    if len(chosen) >= target:
                        break
                    if id(c) not in used and up >= args.orf_fallback and down >= args.orf_fallback:
                        chosen.append((c, up, down, f"{args.orf_fallback}/{args.orf_fallback}")); used.add(id(c))
            for rank, (c, up, down, tier) in enumerate(chosen, 1):
                picks.append({"role": role, "subclass": sub, "common_name": common, "read_type": rt,
                              "rank": rank, "rt_id": c["id"], "sra_accession": c["accession"],
                              "contig_id": c["contig_id"], "contig_length": c["contig_length"],
                              "rt_start": c["protein_start"], "rt_end": c["protein_end"],
                              "strand": c["strand"], "n_orf_up": up, "n_orf_down": down,
                              "pass_tier": tier, "repr_c90": c["repr_c90"]})
            n_strict = sum(1 for _, _, _, t in chosen if t == f"{args.orf_min}/{args.orf_min}")
            yield_ct[sub][rt] = (n_strict, len(chosen) - n_strict)
            if len(chosen) < target:
                flags.append(f"{sub} ({common}) {rt}: {len(chosen)}/{target} "
                             f"(pool={len(by_bucket.get((sub, rt), []))}, {len(scored)} locatable)")

    # ---- confirm read_type from ENA run metadata, ONLY for the subset that can benefit ----
    # A PacBio/ONT pick is already long by platform and gains nothing from the lookup, so we skip it and query ENA
    # only for platform-SHORT (Illumina-family) picks — the only ones that can hide a merged/synthetic-long library.
    # (We can't pre-restrict further to single-end: the atlas parquet has no library_layout column, so we learn
    # single/paired from the ENA response itself — a PAIRED result just resolves to short via confirm_read_type.)
    for p in picks:
        p["read_type_platform"] = p["read_type"]                          # provenance: the platform-based label
        stats = (ena_read_stats(p["sra_accession"]) or {}) if p["read_type"] == "short" else {}
        p["avg_read_len"] = stats.get("avg_read_len", "")
        p["library_layout"] = stats.get("library_layout", "")
        p["read_type"] = confirm_read_type(p["read_type"], stats)         # truth used for routing + render
    reclassified = [p for p in picks if p["read_type"] != p["read_type_platform"]]

    # ---- emit control_loci.tsv ----
    cols_out = ["role", "subclass", "common_name", "read_type", "read_type_platform", "avg_read_len",
                "library_layout", "rank", "rt_id", "sra_accession", "contig_id", "contig_length",
                "rt_start", "rt_end", "strand", "n_orf_up", "n_orf_down", "pass_tier", "repr_c90"]
    with open(args.out / "control_loci.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols_out, delimiter="\t")
        w.writeheader()
        w.writerows(picks)

    # ---- emit selection_report.md ----
    L = ["# FDZ004.2 selection report", "",
         "## Search counts (distinct repr_c90 clusters)",
         f"- transcriptomic **OR** metatranscriptomic: **{n_or:,}**",
         f"- transcriptomic (has ≥1 T occurrence): **{n_t:,}**",
         f"- metatranscriptomic (has ≥1 MT occurrence): **{n_mt:,}**",
         f"- metatranscriptomic-only (OR − transcriptomic): **{n_or - n_t:,}**",
         f"- transcriptomic-only (OR − metatranscriptomic): **{n_or - n_mt:,}**",
         f"- both T and MT (T + MT − OR): **{n_t + n_mt - n_or:,}**", "",
         "## Short vs long within transcriptomic",
         "| read type | occurrence rows | distinct clusters |", "|---|---|---|"]
    for k in ("short", "long"):
        r, cl = sl.get(k, (0, 0))
        L.append(f"| {k} | {r:,} | {cl:,} |")
    L += ["", "_Bucketing uses the platform label (long = PACBIO_SMRT / OXFORD_NANOPORE; short = everything else). "
          "The emitted `read_type` column is then CONFIRMED against ENA run metadata (`avg_read_len`, `library_layout`): "
          "a single-end library measuring ≥ %d bp is re-labelled long (catches merged/synthetic-long ILLUMINA). "
          "The pipeline's map-time measured-length check remains the ground-truth override._" % LONG_READ_BP]
    if reclassified:
        L += ["", "### read_type reclassified by ENA metadata (platform-short → measured-long)",
              "| subclass | accession | platform | avg_read_len | layout |", "|---|---|---|---|---|"]
        L += [f"| {p['subclass']} | {p['sra_accession']} | {p['read_type_platform']} | "
              f"{p['avg_read_len']} | {p['library_layout']} |" for p in reclassified]
    L += [
          "", f"## Per-subclass selection (strict {args.orf_min}/{args.orf_min}, else {args.orf_fallback}/{args.orf_fallback} fallback for short; target {args.n_short} short + {args.n_long} long)",
          "| role | subclass | common | short (strict + fb) | long |", "|---|---|---|---|---|"]
    for sub, role, common in PANEL:
        ss, sf = yield_ct[sub]["short"]; ls, lf = yield_ct[sub]["long"]
        L.append(f"| {role} | {sub} | {common} | {ss + sf}/{args.n_short} ({ss} strict + {sf} fb) | {ls + lf}/{args.n_long} |")
    L += ["", "## ORF-neighbourhood distribution (all transcriptomic candidates checked, up/down)"]
    for sub, role, common in PANEL:
        d = orf_dist.get(sub, [])
        pass33 = sum(1 for u, dn in d if u >= args.orf_min and dn >= args.orf_min)
        summ = ", ".join(f"{u}/{dn}" for u, dn in sorted(d, reverse=True)[:12]) or "none"
        L.append(f"- **{sub}** ({common}): {len(d)} candidates checked, {pass33} pass {args.orf_min}/{args.orf_min}; top (up/down): {summ}")
    if flags:
        L += ["", "## FLAGS (buckets short of target)"] + [f"- {f}" for f in flags]
    L += ["", f"Total loci selected: **{len(picks)}** ({sum(1 for p in picks if p['read_type']=='short')} short, "
          f"{sum(1 for p in picks if p['read_type']=='long')} long).", ""]
    (args.out / "selection_report.md").write_text("\n".join(L) + "\n")
    print(f"wrote {args.out/'control_loci.tsv'} ({len(picks)} loci) and {args.out/'selection_report.md'}")
    print(f"counts: OR={n_or} T={n_t} MT={n_mt} | short={sl.get('short')} long={sl.get('long')}")


if __name__ == "__main__":
    main()
