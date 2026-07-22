#!/usr/bin/env python3
"""FDZ004.2 RNA-seq pilot pipeline — transcriptomic evidence for RT loci.

For each SRA accession behind a pilot locus:
  fetch_sra -> extract_fastq -> build_reference -> map_reads -> postprocess_bam -> coverage

Design (from FDZ004_brainstorm.md, confirmed 2026-07-13):
  * Input = transcriptomic contigs from the RT Atlas that carry an RT of interest.
    Each contig is "perfectly paired" to its source SRA run, so we map that run's
    reads back onto the contig and read out coverage over the RT + its IGRs.
  * Map the FULL FASTQ directly (no BLAST/k-mer pre-subset — the aligner's own
    k-mer seeding handles that).
  * Per-accession aligner routing: short reads -> coverage (FDZ004_SHORT_ALIGNERS), long reads -> minimap2 splice.
  * Resumable: per-accession state.json records each step's status/timing/outputs;
    a rerun skips completed steps and resumes from the first unfinished one.
  * Multi-mapping is annotated, never filtered; only *unmapped* reads are dropped
    from the BAM to keep it light. CIGAR / clipping / indels are preserved.

Coding rules honoured: deterministic ordering (sorted ids/contigs, coordinate-sorted
BAM); fail-loud (missing lookups raise, every subprocess is check=True); provenance
(each artifact + its command/params/versions recorded in state.json). The ERDA base
URL is a bearer credential, read from the environment, never written to any output.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import erda_extract as ee  # noqa: E402  (local sibling helper)
from select_loci import sanitize  # noqa: E402  (organized output path <role>/<subclass>)

PIPELINE_VERSION = "0.2.0"
STEPS = ["fetch_sra", "extract_fastq", "build_reference", "map_reads", "postprocess_bam", "coverage"]


# --------------------------------------------------------------------------- config
@dataclass(frozen=True)
class Config:
    erda_base: str            # bearer; reference-building only
    work_dir: Path
    out_dir: Path
    max_sra_gb: float = 20.0
    threads: int = 8
    igr_max: int = 3
    igr_min: int = 1
    window_fallback_bp: int = 5000
    sra_s3_prefix: str = "s3://sra-pub-run-odp/sra"
    cleanup: bool = True      # after coverage, delete .sra/FASTQ/raw-BAM (keep coverage+mapped BAM)
    short_aligners: tuple = ("bowtie2",)    # coverage arm (short reads); verdict = bowtie2. benchmark: bwa-mem2,bowtie2,minimap2-sr,minibwa
    min_aln_frac: float = 0.0               # coverage filter OFF by default (transcriptomic, aligner-invariant); metaT: ~0.90
    min_mapq: int = 0                       # coverage filter OFF by default; metaT specificity: ~20 (drops multimappers)
    minimal: bool = False                   # deployment/scale-up: depth over the RT-neighborhood window only; skip the unused calmd/features/aligner-comparison artifacts

    @staticmethod
    def from_env(work_dir: Path, out_dir: Path) -> "Config":
        base = os.environ.get("RT_ATLAS_ERDA_BASE", "").rstrip("/")
        if not base:
            raise RuntimeError("RT_ATLAS_ERDA_BASE is unset — source config.local.sh first")
        return Config(
            erda_base=base,
            work_dir=work_dir,
            out_dir=out_dir,
            max_sra_gb=float(os.environ.get("FDZ004_MAX_SRA_GB", 20)),
            threads=int(os.environ.get("FDZ004_THREADS", 8)),
            igr_max=int(os.environ.get("FDZ004_IGR_MAX", 3)),
            igr_min=int(os.environ.get("FDZ004_IGR_MIN", 1)),
            window_fallback_bp=int(os.environ.get("FDZ004_WINDOW_FALLBACK_BP", 5000)),
            sra_s3_prefix=os.environ.get("FDZ004_SRA_S3_PREFIX", "s3://sra-pub-run-odp/sra"),
            cleanup=os.environ.get("FDZ004_CLEANUP", "1") not in ("0", "false", "False", ""),
            short_aligners=tuple(a.strip() for a in os.environ.get("FDZ004_SHORT_ALIGNERS", "bowtie2").split(",") if a.strip()),
            min_aln_frac=float(os.environ.get("FDZ004_MIN_ALN_FRAC", 0.0)),
            min_mapq=int(os.environ.get("FDZ004_MIN_MAPQ", 0)),
            minimal=os.environ.get("FDZ004_MINIMAL", "0") not in ("0", "false", "False", ""),
        )


class SkipOversized(Exception):
    """Raised when an .sra archive exceeds max_sra_gb; the accession is skipped, not failed."""


# --------------------------------------------------------------------------- utilities
def log(msg: str) -> None:
    print(msg, flush=True)


def _run(cmd: list[str], *, log_path: Path | None = None) -> float:
    """Run cmd with check=True; append stdout+stderr to log_path; return wall seconds.

    A non-zero exit raises CalledProcessError and stops the step (fail loud) — so a
    failed aligner can never leave a stale BAM that the next step silently consumes.
    """
    t0 = time.time()
    if log_path is not None:
        with open(log_path, "ab") as fh:
            fh.write(f"\n$ {' '.join(map(str, cmd))}\n".encode())
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, check=True)
    return round(time.time() - t0, 2)


def _need(tool: str) -> str:
    p = shutil.which(tool)
    if p is None:
        raise RuntimeError(f"required tool not on PATH: {tool}")
    return p


def _tool_version(tool: str) -> str:
    try:
        out = subprocess.run([tool, "--version"], capture_output=True, text=True)
        return (out.stdout or out.stderr).splitlines()[0].strip()
    except Exception:
        return "unknown"


def _gnu_time() -> str | None:
    """GNU time (not the shell builtin) for CPU%/peak-RSS capture, if present."""
    cand = shutil.which("gtime") or ("/usr/bin/time" if Path("/usr/bin/time").exists() else None)
    return cand


def _parse_gnu_time(path: Path) -> dict:
    rss_kb = cpu_pct = None
    for line in path.read_text().splitlines():
        if "Maximum resident set size" in line:
            rss_kb = int(line.rsplit(":", 1)[-1].strip())
        elif "Percent of CPU this job got" in line:
            cpu_pct = line.rsplit(":", 1)[-1].strip()
    return {"max_rss_kb": rss_kb, "cpu_pct": cpu_pct}


def _proc_subtree(root_pid: int) -> list[int]:
    """PIDs in root_pid's subtree via /proc PPID links (Linux only; [] elsewhere).

    Needed because GNU `time -v` on a launcher script (e.g. the `bowtie2` Perl wrapper)
    reports the wrapper's RSS, not the real `bowtie2-align-s` child that dominates memory.
    """
    proc = Path("/proc")
    if not proc.exists():
        return []
    children: dict[int, list[int]] = {}
    for entry in proc.iterdir():
        if not entry.name.isdigit():
            continue
        try:  # stat field 4 = ppid; comm (field 2) may hold spaces/parens -> split after ')'
            data = (entry / "stat").read_text()
            ppid = int(data[data.rindex(")") + 1:].split()[1])
        except (OSError, ValueError):
            continue
        children.setdefault(ppid, []).append(int(entry.name))
    out, stack = [], [root_pid]
    while stack:
        p = stack.pop()
        out.append(p)
        stack.extend(children.get(p, []))
    return out


def _sample_peak_rss(root_pid: int, stop: threading.Event, result: dict, interval: float = 0.2) -> None:
    """Best-effort: poll summed RSS of root_pid's whole subtree, record the peak (KB).

    Non-fatal by contract — any error just yields a lower/absent estimate, never crashes
    the mapping it measures. Captures the aligner's real child processes (see _proc_subtree).
    """
    page_kb = os.sysconf("SC_PAGE_SIZE") // 1024
    peak = 0
    while not stop.is_set():
        total = 0
        for p in _proc_subtree(root_pid):
            try:
                rss_pages = int((Path("/proc") / str(p) / "statm").read_text().split()[1])
                total += rss_pages * page_kb
            except (OSError, ValueError, IndexError):
                continue
        peak = max(peak, total)
        stop.wait(interval)
    result["peak_rss_kb"] = peak


def _duckdb_con(erda_base: str):
    import duckdb  # lazy: only reference-building needs it, not arg-parsing/loci-loading
    con = duckdb.connect()
    con.execute("INSTALL httpfs; LOAD httpfs; SET enable_progress_bar=false;")
    return con


# --------------------------------------------------------------------------- state
def load_state(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {"pipeline_version": PIPELINE_VERSION, "steps": {}}


def save_state(path: Path, state: dict) -> None:
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(state, indent=2, default=str) + "\n")
    tmp.replace(path)  # atomic


def step_done(state: dict, step: str) -> bool:
    return state["steps"].get(step, {}).get("status") == "done"


# --------------------------------------------------------------------------- loci
def load_loci(loci_tsv: Path) -> dict[str, list[dict]]:
    """control_loci.tsv -> {sra_accession: [locus_row, ...]} (deterministic)."""
    by_sra: dict[str, list[dict]] = {}
    with open(loci_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_sra.setdefault(row["sra_accession"], []).append(row)
    for sra in by_sra:
        by_sra[sra].sort(key=lambda r: (r["contig_id"], int(r["rt_start"])))
    return dict(sorted(by_sra.items()))


# --------------------------------------------------------------------------- steps
def fetch_sra(cfg: Config, sra: str, dest: Path, log_path: Path) -> dict:
    dest.mkdir(parents=True, exist_ok=True)
    src = f"{cfg.sra_s3_prefix}/{sra}/{sra}"
    _need("aws")
    # size check first (public bucket -> --no-sign-request)
    ls = subprocess.run(
        ["aws", "s3", "ls", "--no-sign-request", src],
        check=False, capture_output=True, text=True,
    )
    sra_path = dest / f"{sra}.sra"
    used = "s3"
    if ls.returncode == 0 and ls.stdout.split():
        nbytes = int(ls.stdout.split()[2])
        if nbytes > cfg.max_sra_gb * 1e9:
            raise SkipOversized(f"{sra}: .sra is {nbytes/1e9:.1f} GB > {cfg.max_sra_gb} GB cap")
        _run(["aws", "s3", "cp", "--no-sign-request", src, str(sra_path)], log_path=log_path)
    else:  # fall back to NCBI prefetch
        used = "prefetch"
        _need("prefetch")
        _run(["prefetch", "--max-size", f"{int(cfg.max_sra_gb)}G", "-O", str(dest), sra], log_path=log_path)
        found = sorted(dest.rglob(f"{sra}.sra"))
        if not found:
            raise RuntimeError(f"{sra}: prefetch produced no .sra under {dest}")
        sra_path = found[0]
    return {"sra_path": str(sra_path), "bytes": sra_path.stat().st_size, "source": used}


def _median_read_len(fastq: str, sample: int = 2000) -> int:
    """Median length of the first `sample` reads — used to route actual long reads to minimap2
    even when the SRA is platform-labeled ILLUMINA (merged/synthetic-long libraries)."""
    lens = []
    with open(fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                lens.append(len(line.rstrip()))
                if len(lens) >= sample:
                    break
    lens.sort()
    return lens[len(lens) // 2] if lens else 0


def extract_fastq(cfg: Config, sra: str, sra_path: Path, dest: Path, log_path: Path) -> dict:
    dest.mkdir(parents=True, exist_ok=True)
    _need("fasterq-dump")
    # --split-3 keeps mate1/mate2 SYNCED (orphan/single reads go to a separate <sra>.fastq) — avoids the
    # --split-files desync that makes bwa-mem2 reject the pair ("reads have different names"). Paired
    # mapping uses only _1/_2; the orphan file is ignored.
    _run(["fasterq-dump", str(sra_path), "--split-3", "-e", str(cfg.threads),
          "-O", str(dest), "--temp", str(dest)], log_path=log_path)
    r1, r2 = dest / f"{sra}_1.fastq", dest / f"{sra}_2.fastq"
    if r1.exists() and r2.exists():
        fastqs, paired = [str(r1), str(r2)], True          # ignore orphan <sra>.fastq for paired mapping
    else:
        fastqs, paired = sorted(str(p) for p in dest.glob(f"{sra}*.fastq")), False
    if not fastqs:
        raise RuntimeError(f"{sra}: fasterq-dump produced no FASTQ in {dest}")
    total = sum(Path(p).stat().st_size for p in fastqs)
    return {"fastqs": fastqs, "paired": paired, "total_bytes": total,
            "read_len_median": _median_read_len(fastqs[0])}


def _pick_igrs(genes: list[dict], rt: dict, side: str, n: int) -> list[dict]:
    """Return up to n intergenic regions on `side` ('up'/'down') of the RT gene.
    Genes are start-sorted; an IGR is the gap between two consecutive genes."""
    idx = genes.index(rt)
    igrs: list[dict] = []
    if side == "up":
        for j in range(idx, 0, -1):
            gap_start, gap_end = genes[j - 1]["end"], genes[j]["start"]
            if gap_end > gap_start:
                igrs.append({"start": gap_start, "end": gap_end, "side": "up", "n": len(igrs) + 1})
            if len(igrs) >= n:
                break
    else:
        for j in range(idx, len(genes) - 1):
            gap_start, gap_end = genes[j]["end"], genes[j + 1]["start"]
            if gap_end > gap_start:
                igrs.append({"start": gap_start, "end": gap_end, "side": "down", "n": len(igrs) + 1})
            if len(igrs) >= n:
                break
    return igrs


def _load_genes(con, gff_url, cid):
    """Genes for a contig from a local cache (FDZ004_GENES_CACHE/<contig>.tsv) if present, else ERDA.
    The cache avoids 16 concurrent accessions each full-scanning the remote GFF parquet (httpfs timeouts)."""
    cache = os.environ.get("FDZ004_GENES_CACHE")
    if cache:
        p = Path(cache) / f"{cid}.tsv"
        if p.exists():
            out = []
            for line in p.read_text().splitlines()[1:]:
                f = line.split("\t")
                out.append({"contig": f[0], "start": int(f[1]), "end": int(f[2]),
                            "strand": f[3], "type": f[4], "protein_id": f[5]})
            return out
    return ee.gff_rows_for_contig(con, gff_url, cid)


def build_reference(cfg: Config, sra: str, loci: list[dict], dest: Path, con) -> dict:
    """Reference = the contig(s) carrying this SRA's RT loci; annotate RT gene +
    flanking IGRs (>=igr_min, <=igr_max per side) into regions.bed. If a side lacks
    enough flanking genes, pad the RT by window_fallback_bp and flag it."""
    dest.mkdir(parents=True, exist_ok=True)
    contigs = sorted({r["contig_id"] for r in loci})
    ref_fasta = dest / "reference.fasta"
    fna = f"{cfg.erda_base}/all_contigs.fna"
    fai = f"{cfg.erda_base}/all_contigs.fna.fai.parquet"
    gff = f"{cfg.erda_base}/all_meta_ncbi_gff.parquet"
    lengths = ee.extract_contigs(con, fna, fai, contigs, ref_fasta)
    genes_by_contig = {c: _load_genes(con, gff, c) for c in contigs}  # local cache if set (avoids ERDA contention), else ERDA

    bed_lines: list[str] = []
    per_locus: list[dict] = []
    warnings: list[str] = []
    for loc in sorted(loci, key=lambda r: (r["contig_id"], int(r["rt_start"]))):
        cid, rt_id = loc["contig_id"], loc["rt_id"]
        rt_start, rt_end = int(loc["rt_start"]), int(loc["rt_end"])
        genes = genes_by_contig[cid]
        # locate the RT gene: by protein_id, else by coordinate overlap
        rt_gene = next((g for g in genes if str(g["protein_id"]) == rt_id), None)
        if rt_gene is None:
            rt_gene = next((g for g in genes if not (g["end"] < rt_start or g["start"] > rt_end)), None)
        if rt_gene is None:
            raise RuntimeError(f"{sra}/{cid}: RT gene {rt_id} not found in GFF (by id or overlap)")
        up = _pick_igrs(genes, rt_gene, "up", cfg.igr_max)
        down = _pick_igrs(genes, rt_gene, "down", cfg.igr_max)
        if len(up) < cfg.igr_min or len(down) < cfg.igr_min:
            warnings.append(f"{cid}: only {len(up)} up / {len(down)} down IGR(s) "
                            f"(need >={cfg.igr_min}); padding RT by {cfg.window_fallback_bp} bp")
            pad = cfg.window_fallback_bp
            up = up or [{"start": max(0, rt_gene["start"] - pad), "end": rt_gene["start"],
                         "side": "up", "n": 1, "fallback": True}]
            down = down or [{"start": rt_gene["end"], "end": min(lengths[cid], rt_gene["end"] + pad),
                             "side": "down", "n": 1, "fallback": True}]
        # regions.bed: RT locus + IGRs in GFF-native coords (see coverage caveat in README)
        bed_lines.append(f"{cid}\t{rt_gene['start']}\t{rt_gene['end']}\tRT|{rt_id}\t0\t{rt_gene['strand']}")
        for igr in up + down:
            bed_lines.append(f"{cid}\t{igr['start']}\t{igr['end']}\tIGR_{igr['side']}{igr['n']}|{rt_id}\t0\t.")
        per_locus.append({"rt_id": rt_id, "contig_id": cid, "contig_len": lengths[cid],
                          "rt_gene": {k: rt_gene[k] for k in ("start", "end", "strand")},
                          "igr_up": len(up), "igr_down": len(down)})

    (dest / "regions.bed").write_text("\n".join(bed_lines) + "\n")
    # persist the FULL per-contig gene list so the neighborhood plots read it here (no separate ERDA fetch)
    glines = ["\t".join(["contig", "start", "end", "strand", "type", "protein_id"])]
    for c in contigs:
        for g in genes_by_contig[c]:
            glines.append(f"{c}\t{g['start']}\t{g['end']}\t{g['strand']}\t{g['type']}\t{g['protein_id']}")
    (dest / "genes.tsv").write_text("\n".join(glines) + "\n")
    manifest = {"contigs": {c: lengths[c] for c in contigs}, "loci": per_locus, "warnings": warnings}
    (dest / "reference_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return {"ref_fasta": str(ref_fasta), "regions_bed": str(dest / "regions.bed"),
            "genes_tsv": str(dest / "genes.tsv"), "n_contigs": len(contigs), "warnings": warnings}


def _map_one(cfg: Config, aligner: str, ref_fasta: Path, fastqs: list[str], paired: bool,
             dest: Path, log_path: Path) -> dict:
    dest.mkdir(parents=True, exist_ok=True)
    _need("samtools")
    raw = dest / "aligned.sorted.bam"
    if aligner == "bowtie2":
        _need("bowtie2"); _need("bowtie2-build")
        idx = dest / "bt2_index"
        idx_s = _run(["bowtie2-build", "--threads", str(cfg.threads), str(ref_fasta), str(idx)], log_path=log_path)
        reads = (["-1", fastqs[0], "-2", fastqs[1]] if paired else ["-U", ",".join(fastqs)])
        # --local: soft-clip read ends (capture junction/circular/repetitive-product boundaries),
        # the fair comparison against bwa-mem2's default soft-clipping.
        map_cmd = ["bowtie2", "--local", "-x", str(idx), *reads, "-p", str(cfg.threads)]
    elif aligner == "bwa-mem2":
        _need("bwa-mem2")
        _run(["bwa-mem2", "index", str(ref_fasta)], log_path=log_path)
        idx_s = 0.0  # index time folded into the command above; measured separately below if needed
        map_cmd = ["bwa-mem2", "mem", "-t", str(cfg.threads), str(ref_fasta), *fastqs]
    elif aligner == "minimap2":
        _need("minimap2")  # splice-aware long-read mapping (STRUCTURE arm; ONT/PacBio-style libraries)
        idx_s = 0.0
        map_cmd = ["minimap2", "-ax", "splice", "-t", str(cfg.threads), str(ref_fasta), *fastqs]
    elif aligner == "minimap2-sr":
        _need("minimap2")  # short-read preset (COVERAGE arm candidate in the 3-way benchmark)
        idx_s = 0.0
        map_cmd = ["minimap2", "-ax", "sr", "-t", str(cfg.threads), str(ref_fasta), *fastqs]
    elif aligner == "minibwa":
        _need("minibwa")  # bwa-mem successor (BWT + minimap2-style chaining); COVERAGE arm candidate.
        # No spliced alignment -> short arm only. `map ref r1 [r2]` (bwa CLI); no extra tags, for a fair cost
        # comparison vs the others (none compute MD in their map step).
        idx_s = _run(["minibwa", "index", str(ref_fasta)], log_path=log_path)
        map_cmd = ["minibwa", "map", "-t", str(cfg.threads), str(ref_fasta), *fastqs]
    else:
        raise RuntimeError(f"unknown aligner: {aligner}")
    # wrap the aligner with GNU time (-o file) to capture CPU% + peak RSS without
    # polluting the aligner's own stderr log. map | sort -> coordinate-sorted BAM.
    gtime = _gnu_time()
    timefile = dest / "map.time.txt"
    run_cmd = ([gtime, "-v", "-o", str(timefile)] + map_cmd) if gtime else map_cmd
    samp: dict = {}
    stop = threading.Event()
    sampler = None
    t0 = time.time()
    with open(log_path, "ab") as lg:
        lg.write(f"\n$ {' '.join(run_cmd)} | samtools view -F4 | samtools sort\n".encode())
        aln = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=lg)
        if Path("/proc").exists():  # process-tree peak-RSS sampler (captures the aligner's children)
            sampler = threading.Thread(target=_sample_peak_rss, args=(aln.pid, stop, samp), daemon=True)
            sampler.start()
        # Drop unmapped reads BEFORE sort: on a locus where few reads map (e.g. a mis-paired contig) the
        # aligner still streams every unmapped read; sorting those wastes memory and crashed bwa-mem2's sort
        # on one junk locus. `-m` caps the sort buffer (spills to disk). NB: the aligner's OWN memory is
        # separate (sampled above) and is NOT bounded by this — a hard per-map cap would be needed for that.
        vf = subprocess.Popen(["samtools", "view", "-u", "-F", "4", "-"], stdin=aln.stdout,
                              stdout=subprocess.PIPE, stderr=lg)
        aln.stdout.close()
        subprocess.run(["samtools", "sort", "-m", "768M", "-@", str(cfg.threads), "-o", str(raw), "-"],
                       stdin=vf.stdout, check=True, stderr=lg)
        vf.stdout.close()
        vfrc = vf.wait()
        rc = aln.wait()
        stop.set()
        if sampler is not None:
            sampler.join(timeout=2)
        if rc != 0:
            raise subprocess.CalledProcessError(rc, map_cmd)
        if vfrc != 0:
            raise subprocess.CalledProcessError(vfrc, ["samtools", "view", "-F4"])
    map_s = round(time.time() - t0, 2)
    _run(["samtools", "index", str(raw)], log_path=log_path)
    bench = {"index_build_s": idx_s, "map_s": map_s}
    if gtime and timefile.exists():
        gt = _parse_gnu_time(timefile)
        bench["gnu_max_rss_kb"] = gt.get("max_rss_kb")  # wrapper-level; under-reports wrapped aligners
        bench["cpu_pct"] = gt.get("cpu_pct")
    # authoritative peak RSS = process-tree sampler if it ran, else GNU time's (weaker) figure
    bench["peak_rss_kb"] = samp.get("peak_rss_kb") or bench.get("gnu_max_rss_kb")
    return {"raw_bam": str(raw), **bench}


LONG_READ_BP = 500   # median read length above which we treat a library as long-read → minimap2

def aligners_for_loci(cfg: Config, loci: list[dict], read_len_median: int | None = None) -> list[str]:
    """Per-accession aligner routing (coverage/structure pivot). Actual long reads (median > LONG_READ_BP —
    catches platform-labeled ILLUMINA merged/synthetic-long) → minimap2 (splice, the STRUCTURE arm). Else by
    control_loci read_type: long → minimap2; short → cfg.short_aligners (the COVERAGE arm; one aligner in
    production, the benchmark set during a bake-off). read_type is REQUIRED — a missing one raises (fail loudly)."""
    if read_len_median and read_len_median > LONG_READ_BP:
        return ["minimap2"]
    if not all(l.get("read_type") for l in loci):
        raise RuntimeError(f"locus missing read_type (regenerate control_loci.tsv via select_loci.py): "
                           f"{loci[0].get('sra_accession')}")
    if "long" in {(l.get("read_type") or "").lower() for l in loci}:
        return ["minimap2"]
    return list(cfg.short_aligners)


def map_reads(cfg: Config, aligners: list[str], ref_fasta: Path, fastqs: list[str], paired: bool,
              dest: Path, log_path: Path) -> dict:
    """Map with each requested aligner, isolating per-aligner failure. If one aligner
    crashes (e.g. bowtie2 --local OOM-aborts under the ulimit guard), the others are kept
    so the benchmark still has data; only if ALL fail does the step fail."""
    out: dict = {}
    for a in aligners:
        try:
            out[a] = _map_one(cfg, a, ref_fasta, fastqs, paired, dest / a, log_path)
        except subprocess.CalledProcessError as exc:
            log(f"  map_reads: aligner {a!r} FAILED ({exc!r}) — continuing; other aligners preserved")
    if not out:
        raise RuntimeError("all configured aligners failed during mapping")
    return out


def _flagstat(bam: Path) -> dict:
    out = subprocess.run(["samtools", "flagstat", "-O", "tsv", str(bam)],
                         check=True, capture_output=True, text=True).stdout
    stats = {}
    for line in out.splitlines():
        passed, _failed, name = line.split("\t")
        stats[name.strip()] = int(passed) if passed.isdigit() else passed  # % / N/A stay strings
    return stats


def postprocess_bam(cfg: Config, mapresult: dict, dest: Path, log_path: Path) -> dict:
    """Drop unmapped reads (-F 4), keep everything else (CIGAR/clipping/tags/multimap).
    Record mapping stats per aligner."""
    out: dict = {}
    for aligner, mr in mapresult.items():
        adest = dest / aligner
        adest.mkdir(parents=True, exist_ok=True)
        raw = Path(mr["raw_bam"])
        mapped = adest / "mapped.bam"
        _run(["samtools", "view", "-b", "-F", "4", "-@", str(cfg.threads), "-o", str(mapped), str(raw)],
             log_path=log_path)
        _run(["samtools", "index", str(mapped)], log_path=log_path)
        (adest / "flagstat.tsv").write_text(
            subprocess.run(["samtools", "flagstat", "-O", "tsv", str(mapped)],
                           check=True, capture_output=True, text=True).stdout)
        _run(["bash", "-c", f"samtools idxstats {mapped} > {adest/'idxstats.tsv'}"], log_path=log_path)
        fs = _flagstat(mapped)
        total = fs.get("primary", fs.get("total", 0))
        secondary = fs.get("secondary", 0)
        out[aligner] = {"mapped_bam": str(mapped),
                        "primary_mapped": fs.get("primary mapped", 0),
                        "secondary": secondary,
                        "mapping_rate_pct": fs.get("primary mapped %", None),
                        "total_alignments": total,
                        # map-step benchmark, carried from _map_one for the aligner comparison
                        "peak_rss_kb": mr.get("peak_rss_kb"),
                        "gnu_max_rss_kb": mr.get("gnu_max_rss_kb"),
                        "cpu_pct": mr.get("cpu_pct"),
                        "map_s": mr.get("map_s"),
                        "index_build_s": mr.get("index_build_s")}
    return out


def _filter_for_coverage(bam: Path, out_bam: Path, min_aln_frac: float, min_mapq: int) -> tuple[int, int]:
    """Coverage-specificity pre-filter (trim-free metatranscriptomic specificity): keep a primary mapped
    read only if MAPQ >= min_mapq AND aligned-fraction (1 - soft-clip/read-len, from CIGAR — tool-agnostic)
    >= min_aln_frac. Emulates end-to-end specificity (rejects partial / cross-mapped / ambiguous reads)
    without a trim step. COVERAGE DEPTH ONLY — the structural pileup uses the unfiltered BAM. Returns
    (reads_kept, reads_considered)."""
    import pysam
    src = pysam.AlignmentFile(str(bam), "rb")
    dst = pysam.AlignmentFile(str(out_bam), "wb", template=src)
    kept = tot = 0
    for r in src.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        tot += 1
        if r.mapping_quality < min_mapq:
            continue
        ql = r.query_length or 0                                        # incl soft-clips, excl hard-clips
        soft = sum(l for op, l in (r.cigartuples or []) if op == 4)
        if ql and (ql - soft) / ql >= min_aln_frac:
            dst.write(r); kept += 1
    dst.close(); src.close()
    return kept, tot


def coverage(cfg: Config, postresult: dict, ref_manifest: dict, regions_bed: Path,
             ref_fasta: Path, dest: Path, log_path: Path) -> dict:
    """Per-aligner: depth quantification (enriched IGRs) + RT-biology feature extraction
    (clip/split/mismatch/indel/multimap via features.py) + a per-RT aligner comparison.
    In minimal mode (cfg.minimal, for the scale-up deployment): depth is restricted to the padded
    RT-neighborhood window (smaller depth.tsv.gz, less RAM; scores unchanged — see _write_depth_window),
    and the metric-irrelevant calmd/features/aligner-comparison artifacts are skipped."""
    regions = _load_regions(regions_bed)
    out: dict = {}
    region_arg = ""
    if cfg.minimal:
        dest.mkdir(parents=True, exist_ok=True)
        contig_lens = {c: int(n) for c, n in (ref_manifest.get("contigs") or {}).items()}
        window_bed = dest / "depth_window.bed"
        _write_depth_window(regions, contig_lens, DEPTH_WINDOW_MARGIN_BP, window_bed)
        region_arg = f"-b {window_bed} "
    for aligner, pr in postresult.items():
        adest = dest / aligner
        adest.mkdir(parents=True, exist_ok=True)
        bam = Path(pr["mapped_bam"])
        # MD-tagged BAM so features.py can call mismatches (A→N for DGR) — features unused in minimal mode
        md_bam = adest / "mapped.md.bam"
        if not cfg.minimal:
            _run(["bash", "-c", f"samtools calmd -b {bam} {ref_fasta} > {md_bam} 2>/dev/null"], log_path=log_path)
            _run(["samtools", "index", str(md_bam)], log_path=log_path)
        # depth quantification (enriched IGR calls) — COVERAGE-SPECIFICITY pre-filter first:
        # keep only reads with MAPQ >= min_mapq AND aligned-fraction >= min_aln_frac (trim-free
        # metatranscriptomic specificity). Structural analysis below uses the UNfiltered md_bam.
        depth_tsv = adest / "depth.tsv"
        filt_bam = adest / "coverage_filtered.bam"
        # long-read (splice) coverage: MAPQ-only (af=0) — long reads legitimately carry big soft-clips, so the
        # short-read aligned-fraction cut would drop ~half of them. Short reads: full AF+MAPQ filter.
        af = 0.0 if aligner == "minimap2" else cfg.min_aln_frac
        if af <= 0 and cfg.min_mapq <= 0:   # out-of-the-box: no coverage filter — count ALL primary mapped reads
            _run(["bash", "-c", f"samtools depth -a -G 0x900 {region_arg}{bam} > {depth_tsv}"], log_path=log_path)  # -G: drop secondary+supplementary
            kept = tot = None
        else:
            kept, tot = _filter_for_coverage(bam, filt_bam, af, cfg.min_mapq)
            _run(["bash", "-c", f"samtools depth -a {region_arg}{filt_bam} > {depth_tsv}"], log_path=log_path)
        depth = _read_depth(depth_tsv)
        summary = _quantify(regions, depth)
        summary["coverage_filter"] = {"min_aln_frac": af, "min_mapq": cfg.min_mapq,
                                      "reads_kept": kept, "reads_total": tot}
        (adest / "coverage_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        # feature extraction (clip/split/mismatch/indel/multimap) — JSON only, NO default PNG (the custom
        # 2-profile viz renders separately). md_bam is UNfiltered so structural reads are preserved.
        if not cfg.minimal:
            import features as feat
            try:
                summary["features"] = feat.analyze(str(md_bam), str(ref_fasta), regions,
                                                   adest / "features.json", None)
            except Exception as exc:
                summary["features_error"] = repr(exc)
        _run(["bash", "-c", f"gzip -f {depth_tsv}"], log_path=log_path)
        filt_bam.unlink(missing_ok=True)   # coverage filter is an intermediate; depth.tsv.gz is the artifact
        out[aligner] = summary
    if not cfg.minimal:   # single-aligner deployment has nothing to compare
        out["_comparison"] = _compare_aligners(out, regions, postresult)
        (dest / "aligner_comparison.json").write_text(json.dumps(out["_comparison"], indent=2) + "\n")
    return out


PRIMARY_METRICS = ["rt_mean_depth", "split_SA_reads", "pct_soft_clipped", "pct_multimapped", "peak_rss_mb"]


def _compare_aligners(out: dict, regions: list[dict], postresult: dict | None = None) -> dict:
    """Per RT gene, benchmark metrics side-by-side across aligners → the winner-picking table.

    Primary axes (Dennis, 2026-07-14): RT mean depth, split_SA reads, % soft-clipped,
    % multimapped, peak RAM. These decide the aligner. The remaining biology metrics
    (mismatch / A-specific-frac / indel) + mapping rate / timing are retained under
    `secondary` — demoted, never dropped. Peak RAM is a per-aligner property of the whole
    mapping run (not per-RT), so it is echoed into every RT row for that aligner.
    """
    postresult = postresult or {}
    aligners = [a for a in out if not a.startswith("_")]

    def peak_rss_mb(a: str):
        kb = postresult.get(a, {}).get("peak_rss_kb")
        return round(kb / 1024, 1) if kb else None

    cmp: dict = {"_primary_metrics": PRIMARY_METRICS, "per_rt": {}}
    for r in [r for r in regions if r["name"].startswith("RT|")]:
        rid = r["name"].split("|", 1)[1]
        cmp["per_rt"][rid] = {}
        for a in aligners:
            cov = out[a].get("per_rt", {}).get(rid, {})
            fr = out[a].get("features", {}).get("per_region", {}).get(r["name"], {})
            pr = postresult.get(a, {})
            cmp["per_rt"][rid][a] = {
                # --- primary benchmark axes (order = PRIMARY_METRICS) ---
                "rt_mean_depth": cov.get("rt_mean_depth"),
                "split_SA_reads": fr.get("split_SA_reads"),
                "pct_soft_clipped": fr.get("pct_soft_clipped"),
                "pct_multimapped": fr.get("pct_multimapped"),
                "peak_rss_mb": peak_rss_mb(a),
                # --- secondary (retained, not decisive) ---
                "secondary": {
                    "mismatch_rate_pct": fr.get("mismatch_rate_pct"),
                    "a_specific_frac": fr.get("a_specific_frac"),
                    "pct_indel_reads": fr.get("pct_indel_reads"),
                    "mapping_rate_pct": pr.get("mapping_rate_pct"),
                    "map_s": pr.get("map_s"),
                    "cpu_pct": pr.get("cpu_pct"),
                },
            }
    return cmp


DEPTH_WINDOW_MARGIN_BP = 200   # minimal mode: pad the per-contig RT-neighborhood depth window. Wide enough that
                               # every scored RT/IGR base is inside it, so quantified depths equal whole-contig.


def _write_depth_window(regions: list[dict], contig_lens: dict[str, int], margin: int, out_path: Path) -> None:
    """One BED interval per contig spanning all its RT/IGR regions, padded by `margin` and clamped to the
    contig length. `samtools depth -b` over this shrinks depth.tsv from whole-contig to the neighborhood;
    _region_depths() defaults uncovered positions to 0 and the margin keeps every scored base covered, so
    the quantified per-region depths are identical to the whole-contig computation."""
    span: dict[str, list[int]] = {}
    for r in regions:
        if r["contig"] in span:
            span[r["contig"]][0] = min(span[r["contig"]][0], r["start"])
            span[r["contig"]][1] = max(span[r["contig"]][1], r["end"])
        else:
            span[r["contig"]] = [r["start"], r["end"]]
    lines = []
    for c, (lo, hi) in sorted(span.items()):
        s = max(0, lo - margin)
        e = hi + margin
        if c in contig_lens:
            e = min(e, contig_lens[c])
        lines.append(f"{c}\t{s}\t{e}")
    out_path.write_text("\n".join(lines) + "\n")


def _load_regions(bed: Path) -> list[dict]:
    regs = []
    for line in bed.read_text().splitlines():
        c, s, e, name, _score, strand = line.split("\t")
        regs.append({"contig": c, "start": int(s), "end": int(e), "name": name, "strand": strand})
    return regs


def _read_depth(tsv: Path) -> dict[str, dict[int, int]]:
    d: dict[str, dict[int, int]] = {}
    with open(tsv) as fh:
        for line in fh:
            c, pos, dp = line.split("\t")
            d.setdefault(c, {})[int(pos)] = int(dp)
    return d


def _region_depths(region: dict, depth: dict[str, dict[int, int]]) -> list[int]:
    cd = depth.get(region["contig"], {})
    # regions.bed carries GFF-native coords (assumed 1-based inclusive; confirm on first
    # run); samtools depth positions are 1-based -> inclusive range here. A consistent
    # 1-bp convention offset affects RT and IGR equally, so relative enrichment is robust.
    return [cd.get(p, 0) for p in range(region["start"], region["end"] + 1)]


def _quantify(regions: list[dict], depth: dict[str, dict[int, int]]) -> dict:
    """Threshold = mean depth over the RT gene; report IGR bases/runs above it."""
    per_rt: dict[str, dict] = {}
    rt_regions = [r for r in regions if r["name"].startswith("RT|")]
    for rt in rt_regions:
        rt_id = rt["name"].split("|", 1)[1]
        rt_dp = _region_depths(rt, depth)
        thr = (sum(rt_dp) / len(rt_dp)) if rt_dp else 0.0
        igrs = [r for r in regions if r["name"].endswith(f"|{rt_id}") and r["name"].startswith("IGR")]
        igr_out = []
        for igr in igrs:
            dp = _region_depths(igr, depth)
            enriched = [1 if x > thr else 0 for x in dp]
            runs, cur = [], 0
            for x in enriched:
                if x:
                    cur += 1
                elif cur:
                    runs.append(cur); cur = 0
            if cur:
                runs.append(cur)
            igr_out.append({"name": igr["name"], "len": len(dp),
                            "mean_depth": round(sum(dp) / len(dp), 2) if dp else 0,
                            "bases_above_thr": sum(enriched),
                            "max_contiguous_run": max(runs) if runs else 0})
        per_rt[rt_id] = {"rt_mean_depth": round(thr, 2), "rt_len": len(rt_dp), "igrs": igr_out}
    return {"threshold_rule": "mean RT-gene depth", "per_rt": per_rt}


# _plot_depth removed post-pivot: the default per-contig coverage PNG is retired. Coverage is now the
# custom single-track profile (gene_plot.py); the depth JSON/tsv.gz data is unchanged. (README v00005)


# --------------------------------------------------------------------------- driver
def _accession_subdir(loci: list[dict]) -> Path:
    """Organized output subpath <role>/<subclass>/<short|long>_read (sanitized registry subclass) from the
    accession's loci, so the tree is final DURING the run, not rearranged at the end. One SRA spanning >1
    role/subclass is rare; file under the first row's and warn."""
    roles = sorted({(l.get("role") or "").strip() for l in loci if l.get("role")})
    subs = sorted({(l.get("subclass") or "").strip() for l in loci if l.get("subclass")})
    if len(roles) > 1 or len(subs) > 1:
        log(f"  [warn] accession spans multiple role/subclass {roles}/{subs}; filing under the first")
    role = (loci[0].get("role") or "misc").strip() or "misc"
    sub = (loci[0].get("subclass") or "misc").strip() or "misc"
    rt = (loci[0].get("read_type") or "short").strip().lower() or "short"
    return Path(sanitize(role)) / sanitize(sub) / f"{rt}_read"


def run_accession(cfg: Config, sra: str, loci: list[dict], con) -> dict:
    """Run all steps for one accession, resuming from state.json. Returns final state."""
    adir = cfg.out_dir / _accession_subdir(loci) / sra   # organized <role>/<subclass>/<read_type>_read/<SRA>, written during the run
    adir.mkdir(parents=True, exist_ok=True)
    state_path = adir / "state.json"
    state = load_state(state_path)
    state.setdefault("sra", sra)
    state.setdefault("n_loci", len(loci))
    log_path = adir / "run.log"

    ctx: dict = {}  # carries prior step outputs within a single run
    for step in STEPS:
        if step_done(state, step):
            log(f"  [{sra}] {step}: done (skip)")
            ctx[step] = state["steps"][step].get("result", {})
            continue
        log(f"  [{sra}] {step}: running")
        rec = {"status": "running", "started": time.strftime("%Y-%m-%dT%H:%M:%S")}
        state["steps"][step] = rec
        save_state(state_path, state)
        t0 = time.time()
        try:
            result = _dispatch(cfg, step, sra, loci, adir, ctx, con, log_path)
        except SkipOversized as exc:
            rec.update(status="skipped", reason=str(exc), elapsed_s=round(time.time() - t0, 2))
            save_state(state_path, state)
            log(f"  [{sra}] {step}: SKIPPED — {exc}")
            state["skipped"] = True
            save_state(state_path, state)
            return state
        except Exception as exc:
            rec.update(status="failed", error=repr(exc), elapsed_s=round(time.time() - t0, 2))
            save_state(state_path, state)
            log(f"  [{sra}] {step}: FAILED — {exc!r}")
            raise
        rec.update(status="done", elapsed_s=round(time.time() - t0, 2), result=result)
        ctx[step] = result
        save_state(state_path, state)
    state["complete"] = True
    save_state(state_path, state)
    if cfg.cleanup:
        state["cleaned"] = _cleanup_accession(adir)
        save_state(state_path, state)
        log(f"  [{sra}] cleanup: removed {len(state['cleaned'])} large intermediate(s)")
    return state


def _cleanup_accession(adir: Path) -> list[str]:
    """After completion, delete large regenerable intermediates (.sra, FASTQ, raw
    sorted BAM incl. unmapped). Keep reference/, mapped BAM, coverage, state, log."""
    removed = []
    for p in (adir / "sra", adir / "fastq"):
        if p.exists():
            shutil.rmtree(p)
            removed.append(str(p))
    mapdir = adir / "map"
    if mapdir.exists():
        for raw in sorted(mapdir.glob("*/aligned.sorted.bam*")):
            raw.unlink()
            removed.append(str(raw))
    return removed


def _dispatch(cfg, step, sra, loci, adir, ctx, con, log_path):
    if step == "fetch_sra":
        return fetch_sra(cfg, sra, adir / "sra", log_path)
    if step == "extract_fastq":
        return extract_fastq(cfg, sra, Path(ctx["fetch_sra"]["sra_path"]), adir / "fastq", log_path)
    if step == "build_reference":
        return build_reference(cfg, sra, loci, adir / "reference", con)
    if step == "map_reads":
        ex = ctx["extract_fastq"]
        return map_reads(cfg, aligners_for_loci(cfg, loci, ex.get("read_len_median")),
                         Path(ctx["build_reference"]["ref_fasta"]), ex["fastqs"], ex["paired"], adir / "map", log_path)
    if step == "postprocess_bam":
        return postprocess_bam(cfg, ctx["map_reads"], adir / "bam", log_path)
    if step == "coverage":
        ref_manifest = json.loads((adir / "reference" / "reference_manifest.json").read_text())
        return coverage(cfg, ctx["postprocess_bam"], ref_manifest,
                        Path(ctx["build_reference"]["regions_bed"]),
                        Path(ctx["build_reference"]["ref_fasta"]), adir / "coverage", log_path)
    raise RuntimeError(f"unknown step: {step}")


def main() -> None:
    ap = argparse.ArgumentParser(description="FDZ004.2 RNA-seq pilot pipeline")
    ap.add_argument("--loci", required=True, type=Path, help="control_loci.tsv")
    ap.add_argument("--work-dir", type=Path, default=Path(os.environ.get("FDZ004_WORK", "./rnaseq-work")))
    ap.add_argument("--out-dir", type=Path, default=None, help="default: <work-dir>/out")
    ap.add_argument("--only-sra", nargs="*", help="restrict to these accessions (default: all in loci)")
    args = ap.parse_args()

    work = args.work_dir.resolve()
    out = (args.out_dir or work / "out").resolve()
    out.mkdir(parents=True, exist_ok=True)
    cfg = Config.from_env(work, out)

    by_sra = load_loci(args.loci)
    if args.only_sra:
        by_sra = {k: v for k, v in by_sra.items() if k in set(args.only_sra)}
        if not by_sra:
            raise RuntimeError(f"--only-sra matched none of {sorted(load_loci(args.loci))}")

    log(f"[cfg] short_aligners={','.join(cfg.short_aligners)} threads={cfg.threads} "
        f"max_sra={cfg.max_sra_gb}GB cov_filter=AF>={cfg.min_aln_frac},MQ>={cfg.min_mapq} "
        f"out={out}\n[cfg] tools: " + ", ".join(
            f"{t}={_tool_version(t)}" for t in ("samtools", "bowtie2", "bwa-mem2", "fasterq-dump")
            if shutil.which(t)))
    log(f"[run] {len(by_sra)} accession(s): {', '.join(by_sra)}")

    con = _duckdb_con(cfg.erda_base)
    results = {}
    for sra, loci in by_sra.items():
        results[sra] = run_accession(cfg, sra, loci, con)
    # per-accession run summary — one file per SRA (NOT a single run_summary.json), so parallel Batch
    # tasks publishing to the same S3 dir never collide/overwrite each other.
    for sra, st in results.items():
        rec = {sra: {"complete": st.get("complete", False), "skipped": st.get("skipped", False)}}
        (out / f"run_summary.{sra}.json").write_text(json.dumps(rec, indent=2, default=str) + "\n")
    log(f"\n[DONE] {len(results)} summary file(s) -> {out}/run_summary.<sra>.json")


if __name__ == "__main__":
    main()
