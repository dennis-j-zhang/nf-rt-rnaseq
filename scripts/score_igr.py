"""FDZ004.2 — peak-IGR IGR/RT scoring (reanalysis of the aligner run's coverage).

For each control locus, over the IGRs in its refined neighbourhood window (from reference/regions.bed), compute
per-base depth statistics — mean / p70 / p80 / p90 / max — each divided by the RT-ORF mean depth, then take the
MAX across the neighbourhood's IGRs (the "peak-IGR") as the locus score, separately for each statistic. The
statistic gradient (mean -> p70 -> p80 -> p90 -> max) trades flat-averaging for peak-sensitivity: a localized
ncRNA is a sharp coverage spike a mean would dilute. Short loci use bowtie2 coverage, long loci minimap2
(labelled). Records which IGR won so a dot traces to one candidate contig+SRA neighbourhood.

Usage:  python scripts/score_igr.py <run_tree> [control_loci.tsv] [out.tsv]
"""
import sys, csv, gzip, json
from pathlib import Path
import numpy as np

TREE = Path(sys.argv[1])
LOCI = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("01_pilot_aligner/control_loci.tsv")
OUT = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("01_pilot_aligner/scoring_benchmark/igr_scores.tsv")
STATS = ["mean", "p70", "p80", "p90", "max"]
EXCLUDE = {"DGR", "CRISPR-associated"}      # REMOVED from all analysis (ambiguous control — pos/neg unclear)
RT_MIN_DEPTH = 5.0                          # below this the RT is ~unexpressed -> IGR/RT unreliable; keep but FLAG
ALN_FOR = {"short": "bowtie2", "long": "minimap2"}   # coverage source per read_type (short verdict / long structure)


def sani(s):                                # matches select_loci.sanitize (inlined to avoid its duckdb import)
    return s.replace("Retron Type ", "retron_").replace("/", "_").replace(" ", "_")


def load_depth(cov_dir):
    """(contig, pos1based) -> depth, from `samtools depth -a` output. Keyed by (contig,pos) so a
    MULTI-contig reference (an SRA whose RT loci span >1 contig, all mapped to one combined reference)
    does not collide positions across contigs."""
    d = {}
    with gzip.open(cov_dir / "depth.tsv.gz", "rt") as fh:
        for line in fh:
            c, p, x = line.rstrip("\n").split("\t")
            d[(c, int(p))] = int(x)
    return d


def igr_intervals(regions_bed, rt_id):
    """(contig, name, start0, end) for the IGR rows of THIS rt_id (col3 = 'IGR_<side><n>|<rt_id>');
    BED is 0-based half-open. Filtering by rt_id scopes scoring to the locus's own contig + flanks even
    when the reference carries other RT loci (multi-contig)."""
    out = []
    for line in open(regions_bed):
        f = line.rstrip("\n").split("\t")
        parts = f[3].split("|")
        if parts[0].startswith("IGR") and len(parts) > 1 and parts[1] == rt_id:
            out.append((f[0], parts[0], int(f[1]), int(f[2])))
    return out


def igr_stat(depths, s):
    a = np.asarray(depths, float)
    if a.size == 0:
        return 0.0
    if s == "mean":
        return float(a.mean())
    if s == "max":
        return float(a.max())
    return float(np.percentile(a, int(s[1:])))     # p70/p80/p90


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out_rows, skipped = [], 0
    for r in csv.DictReader(open(LOCI), delimiter="\t"):
        if r["subclass"] in EXCLUDE:                 # drop ambiguous controls entirely
            continue
        rt_type = (r.get("read_type") or "short").lower()
        base = TREE / sani(r["role"]) / sani(r["subclass"]) / f"{rt_type}_read" / r["sra_accession"]
        regions = base / "reference" / "regions.bed"
        covroot = base / "coverage"
        aln = ALN_FOR.get(rt_type, "bowtie2")
        if not (covroot / aln / "depth.tsv.gz").exists() and covroot.exists():   # fall back to whichever aligner ran
            present = sorted(p.name for p in covroot.iterdir() if (p / "depth.tsv.gz").exists())
            if present:
                aln = present[0]
        cov = covroot / aln
        cs = cov / "coverage_summary.json"
        if not (cov / "depth.tsv.gz").exists() or not cs.exists() or not regions.exists():
            skipped += 1
            continue                                    # locus didn't complete for any aligner
        per_rt = json.load(open(cs)).get("per_rt", {}).get(r["rt_id"])
        if not per_rt:
            skipped += 1
            continue
        rt_depth = float(per_rt.get("rt_mean_depth", 0.0))
        depth = load_depth(cov)
        igrs = igr_intervals(regions, r["rt_id"])
        rec = {"system": r["subclass"], "role": r["role"], "sra": r["sra_accession"], "contig": r["contig_id"],
               "rt_id": r["rt_id"], "aligner": aln, "read_type": rt_type, "rt_mean_depth": round(rt_depth, 2),
               "low_rt": int(rt_depth < RT_MIN_DEPTH)}
        for s in STATS:
            best_ratio, best_igr = None, ""
            for contig, name, s0, e in igrs:            # BED 0-based half-open [s0,e) -> 1-based depth (s0, e]
                v = igr_stat([depth.get((contig, p), 0) for p in range(s0 + 1, e + 1)], s)
                ratio = (v / rt_depth) if rt_depth > 0 else None
                if ratio is not None and (best_ratio is None or ratio > best_ratio):
                    best_ratio, best_igr = ratio, name
            rec[f"score_{s}"] = round(best_ratio, 3) if best_ratio is not None else ""
            rec[f"igr_{s}"] = best_igr
        out_rows.append(rec)
    cols = (["system", "role", "sra", "contig", "rt_id", "aligner", "read_type", "rt_mean_depth", "low_rt"]
            + [f"score_{s}" for s in STATS] + [f"igr_{s}" for s in STATS])
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(out_rows)
    print(f"wrote {OUT} ({len(out_rows)} loci scored, {skipped} skipped for missing coverage)")
    for s in ("mean", "p90", "max"):
        vals = [float(x[f"score_{s}"]) for x in out_rows if x[f"score_{s}"] != ""]
        if vals:
            print(f"  {s}: n={len(vals)} median IGR/RT={np.median(vals):.2f} max={max(vals):.2f}")


if __name__ == "__main__":
    main()
