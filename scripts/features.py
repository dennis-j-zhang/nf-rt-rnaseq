#!/usr/bin/env python3
"""Per-locus alignment-feature extraction for FDZ004.2 — the RT-biology signals.

From a coordinate-sorted, indexed, MD-tagged BAM + reference + regions.bed, quantify
over each RT gene and IGR the features that matter for reverse-transcriptase biology:
  * soft-clipping        — splice / circular-junction / repetitive-product boundaries
  * split / SA reads     — junction-spanning (GII intron splicing, circular RNA/DNA)
  * mismatches + A→N     — DGR adenine-templated mutagenesis (substitution spectrum)
  * indels               — insertions/deletions
  * multimapping         — repetitive cDNA products
Emits a metrics JSON and a feature-coloured coverage PNG: grey depth + one coloured
per-base track per feature (clip=black, SA/split=blue, mismatch=red, indel=magenta,
multimap=green). Requires pysam; reads need MD (run `samtools calmd` upstream).
"""
from __future__ import annotations

import json
from pathlib import Path

CLIP_OPS = {4, 5}        # soft/hard clip
INS_OP, DEL_OP = 1, 2

FEATURE_COLORS = {"clip": "black", "split_SA": "blue", "mismatch": "red",
                  "indel": "magenta", "multimap": "green"}


def _empty_region_metrics() -> dict:
    return {"reads": 0, "primary": 0, "soft_clipped": 0, "split_SA": 0, "secondary": 0,
            "multimapped": 0, "indel_reads": 0, "aligned_bp": 0, "mismatches": 0,
            "mismatch_by_ref": {b: 0 for b in "ACGT"}}


def _read_is_multimap(read) -> bool:
    return read.is_secondary or read.has_tag("XA") or read.has_tag("SA")


def analyze(bam_path: str, ref_fasta: str, regions: list[dict], out_json: Path, out_png: Path | None) -> dict:
    """regions: [{contig,start,end,name,strand}, ...] (RT| and IGR_ names). Returns summary dict."""
    import pysam

    bam = pysam.AlignmentFile(bam_path, "rb")
    per_region: dict[str, dict] = {}
    contigs = sorted({r["contig"] for r in regions})
    # per-base tracks for the plot
    depth = {c: {} for c in contigs}
    track = {c: {f: {} for f in FEATURE_COLORS} for c in contigs}

    def bump(c, feat, pos):
        if pos is not None:
            track[c][feat][pos] = track[c][feat].get(pos, 0) + 1

    # --- per-region metrics ---
    for reg in regions:
        m = _empty_region_metrics()
        c, s, e = reg["contig"], reg["start"], reg["end"]
        for read in bam.fetch(c, max(0, s), e):
            m["reads"] += 1
            if read.is_secondary:
                m["secondary"] += 1
                m["multimapped"] += 1
                continue
            m["primary"] += 1
            cig = read.cigartuples or []
            if any(op in CLIP_OPS and n > 0 for op, n in cig):
                m["soft_clipped"] += 1
            if read.has_tag("SA") or read.is_supplementary:
                m["split_SA"] += 1
            if read.has_tag("XA") or read.has_tag("SA"):
                m["multimapped"] += 1
            if any(op in (INS_OP, DEL_OP) and n > 0 for op, n in cig):
                m["indel_reads"] += 1
            try:
                for qpos, rpos, refbase in read.get_aligned_pairs(with_seq=True):
                    if qpos is None or rpos is None or refbase is None:
                        continue
                    m["aligned_bp"] += 1
                    if refbase.islower():
                        m["mismatches"] += 1
                        rb = refbase.upper()
                        if rb in m["mismatch_by_ref"]:
                            m["mismatch_by_ref"][rb] += 1
            except ValueError:
                pass
        per_region[reg["name"]] = m

    # --- per-base tracks over each contig's region span (for the coloured plot) ---
    for c in contigs:
        spans = [(r["start"], r["end"]) for r in regions if r["contig"] == c]
        lo, hi = min(s for s, _ in spans), max(e for _, e in spans)
        for col in bam.pileup(c, max(0, lo), hi, truncate=True):
            depth[c][col.reference_pos] = col.nsegments
        for read in bam.fetch(c, max(0, lo), hi):
            if _read_is_multimap(read):
                bump(c, "multimap", read.reference_start)
            if read.is_secondary:
                continue
            cig = read.cigartuples or []
            if cig and cig[0][0] in CLIP_OPS:
                bump(c, "clip", read.reference_start)
            if cig and cig[-1][0] in CLIP_OPS:
                bump(c, "clip", read.reference_end)
            if read.has_tag("SA") or read.is_supplementary:
                bump(c, "split_SA", read.reference_start)
            if any(op in (INS_OP, DEL_OP) and n > 0 for op, n in cig):
                bump(c, "indel", read.reference_start)
            try:
                for qpos, rpos, refbase in read.get_aligned_pairs(with_seq=True):
                    if refbase is not None and rpos is not None and refbase.islower() and lo <= rpos < hi:
                        bump(c, "mismatch", rpos)
            except ValueError:
                pass
    bam.close()

    summary = {"per_region": _summarize(per_region)}
    out_json.write_text(json.dumps(summary, indent=2) + "\n")
    if out_png is not None:                      # None -> JSON only (default feature PNG retired post-pivot)
        try:
            _plot(contigs, regions, depth, track, out_png)
            summary["png"] = str(out_png)
        except Exception as exc:
            summary["png_error"] = repr(exc)
    return summary


def _summarize(per_region: dict) -> dict:
    out = {}
    for name, m in per_region.items():
        p = m["primary"] or 1
        abp = m["aligned_bp"] or 1
        out[name] = {
            "reads": m["reads"], "primary": m["primary"], "secondary": m["secondary"],
            "pct_soft_clipped": round(100 * m["soft_clipped"] / p, 2),
            "split_SA_reads": m["split_SA"],
            "pct_multimapped": round(100 * m["multimapped"] / p, 2),
            "pct_indel_reads": round(100 * m["indel_reads"] / p, 2),
            "mismatch_rate_pct": round(100 * m["mismatches"] / abp, 3),
            "mismatch_by_ref": m["mismatch_by_ref"],
            "a_specific_frac": round(m["mismatch_by_ref"]["A"] / (m["mismatches"] or 1), 3),  # DGR signature
        }
    return out


def _plot(contigs, regions, depth, track, png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(contigs), 1, figsize=(12, 3.4 * len(contigs)), squeeze=False)
    for ax, c in zip(axes[:, 0], contigs):
        d = depth[c]
        if d:
            xs = sorted(d)
            ax.fill_between(xs, [d[x] for x in xs], step="mid", color="0.78", linewidth=0, label="depth")
        for r in [r for r in regions if r["contig"] == c]:
            ax.axvspan(r["start"], r["end"], alpha=0.12,
                       color="tab:red" if r["name"].startswith("RT|") else "tab:orange")
        ax2 = ax.twinx()  # feature counts on a second axis
        for feat, colour in FEATURE_COLORS.items():
            t = track[c][feat]
            if not t:
                continue
            xt = sorted(t); yt = [t[x] for x in xt]
            if feat in ("clip", "indel", "split_SA"):
                ax2.vlines(xt, 0, yt, color=colour, lw=0.6, alpha=0.8, label=feat)
            else:  # mismatch, multimap → density line
                ax2.plot(xt, yt, color=colour, lw=0.7, alpha=0.8, label=feat)
        ax.set_title(c, fontsize=9)
        ax.set_ylabel("depth", fontsize=8)
        ax2.set_ylabel("feature count", fontsize=8)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=6, loc="upper right", ncol=3)
    axes[-1, 0].set_xlabel("contig position (bp)", fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=120)
    plt.close(fig)
