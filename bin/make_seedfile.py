"""Build the nf-rt-rnaseq seedfile + per-accession loci TSVs from a control_loci.tsv (e.g. 03_pilot_scaleup/).

The scale-up control_loci is unlabeled (no pos/neg), but rnaseq_pipeline.py builds its output path from
`role`/`subclass`, so we add `role` (default 'scaleup') + `common_name` (= subclass). One loci TSV per
accession lets each Batch task stage only its own loci.

Writes:  <out>/loci/<SRA>.tsv   (one per accession, rnaseq_pipeline schema)
         <out>/seedfile.csv      (header: sra_accession,loci)
Then upload <out>/ to S3 and pass  --seedfile s3://.../seedfile.csv  with --s3-prefix set so the seedfile's
`loci` column points at the uploaded S3 paths (Nextflow stages them per task).

Usage:
  python make_seedfile.py <control_loci.tsv> <out_dir> [--role scaleup] \\
         [--s3-prefix s3://genomics-workflow-core/Results/rt-rnaseq/<project>/input]
"""
import csv
import argparse
from pathlib import Path
from collections import defaultdict

# columns rnaseq_pipeline.py actually reads from the loci TSV:
NEEDED = ["role", "subclass", "common_name", "read_type", "rt_id",
          "sra_accession", "contig_id", "rt_start", "rt_end", "strand"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("control_loci")
    ap.add_argument("out_dir")
    ap.add_argument("--role", default="scaleup")
    ap.add_argument("--s3-prefix", default=None,
                    help="if set, seedfile `loci` = <s3-prefix>/loci/<SRA>.tsv (else the local path)")
    a = ap.parse_args()

    out = Path(a.out_dir)
    (out / "loci").mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(a.control_loci), delimiter="\t"))
    by_sra = defaultdict(list)
    for r in rows:
        by_sra[r["sra_accession"]].append(r)

    seed = []
    for sra, rs in sorted(by_sra.items()):               # sorted() = deterministic
        p = out / "loci" / f"{sra}.tsv"
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=NEEDED, delimiter="\t")
            w.writeheader()
            for r in rs:
                w.writerow({"role": a.role, "subclass": r["subclass"], "common_name": r["subclass"],
                            "read_type": r["read_type"], "rt_id": r["rt_id"], "sra_accession": sra,
                            "contig_id": r["contig_id"], "rt_start": r["rt_start"], "rt_end": r["rt_end"],
                            "strand": r.get("strand", "")})
        loci_ref = (f"{a.s3_prefix.rstrip('/')}/loci/{sra}.tsv") if a.s3_prefix else str(p)
        seed.append({"sra_accession": sra, "loci": loci_ref})

    with open(out / "seedfile.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sra_accession", "loci"])
        w.writeheader()
        w.writerows(seed)
    print(f"wrote {out/'seedfile.csv'} ({len(seed)} accessions, {len(rows)} loci) + {out/'loci'}/")


if __name__ == "__main__":
    main()
