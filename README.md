# nf-rt-rnaseq
**RT-flanking transcription screen on MAF AWS Batch.** For each reverse-transcriptase (RT) locus, map the
contig's own SRA RNA-seq run back onto it and score `(immediate-flank IGR depth-stat ÷ RT-ORF mean depth)_max`
— evidence a flank is transcribed (candidate ncRNA). Patterned on `FischbachLab/nf-ninjamap`.

> **v0.1 — reuses the validated pilot code.** Each Batch task runs `scripts/rnaseq_pipeline.py` (steps 1–6) +
> `scripts/score_igr.py` (multi-contig-corrected). The minimizations from the FDZ004.2 spec — neighborhood-only
> `samtools depth`, drop `features.json`/`calmd` — are **TODOs** (see bottom); v0 is correct, just not yet lean.

## Layout
```
nf-rt-rnaseq/
├── main.nf              # per-accession workflow (RUN_ACCESSION)
├── nextflow.config      # MAF Batch profile (awsbatch, default-maf-pipelines, ECR container, erda_dk secret)
├── Dockerfile           # task image: sra-tools/bowtie2/minimap2/samtools/... + vendored scripts
├── bin/make_seedfile.py # control_loci.tsv -> per-accession loci TSVs + seedfile.csv
├── scripts/             # vendored pipeline code (rnaseq_pipeline.py, score_igr.py, erda_extract.py, ...)
└── README.md
```

## One-time setup

**1. Build + push the container to ECR** (needs Docker + ECR push perms; the account ECR is `458432034220.dkr.ecr.us-west-2.amazonaws.com`):
```bash
cd nf-rt-rnaseq
aws --profile fischbach-lab ecr get-login-password --region us-west-2 \
  | docker login --username AWS --password-stdin 458432034220.dkr.ecr.us-west-2.amazonaws.com
aws --profile fischbach-lab ecr create-repository --repository-name nf-rt-rnaseq --region us-west-2 || true
docker build --platform linux/amd64 -t 458432034220.dkr.ecr.us-west-2.amazonaws.com/nf-rt-rnaseq:latest .
docker push 458432034220.dkr.ecr.us-west-2.amazonaws.com/nf-rt-rnaseq:latest
```

**2. Pipeline repo:** published (public) at **`dennis-j-zhang/nf-rt-rnaseq`** — the `nextflow-production` head runs `nextflow run dennis-j-zhang/nf-rt-rnaseq`. (Transfer into the FischbachLab org later if desired.)

**3. Confirm the TASK role perms** (jobs run as a role, not your user): `secretsmanager:GetSecretValue` on `erda_dk` + write to `s3://genomics-workflow-core/Results/rt-rnaseq/*`. (Ask Xiandong if unsure.)

## Per-run

**4. Make the seedfile + upload input to S3:**
```bash
python bin/make_seedfile.py ../03_pilot_scaleup/control_loci.tsv input \
  --s3-prefix s3://genomics-workflow-core/Results/rt-rnaseq/<project>/input
aws --profile fischbach-lab s3 cp --recursive input s3://genomics-workflow-core/Results/rt-rnaseq/<project>/input/
```

**5. Submit the head job** (start with a 2–3-accession test seedfile before the full 624):
```bash
aws --profile fischbach-lab batch submit-job \
  --job-name nf-rt-rnaseq-<project> \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides command="dennis-j-zhang/nf-rt-rnaseq, --seedfile s3://genomics-workflow-core/Results/rt-rnaseq/<project>/input/seedfile.csv, --project <project>"
```
**Note the returned `jobId`/`jobName`.**

**6. Monitor:** `aws --profile fischbach-lab batch describe-jobs --jobs <jobId> --query 'jobs[0].status'`

## Output
`s3://genomics-workflow-core/Results/rt-rnaseq/<project>/runs/<role>/<subclass>/<read_type>_read/<SRA>/`
— `reference/{genes.tsv,regions.bed}`, `coverage/<aln>/{depth.tsv.gz,coverage_summary.json}`, `state.json`, `igr_scores.<SRA>.tsv`.

## TODO before the full sweep
- [ ] Build + push the ECR container (step 1); push the repo to GitHub (step 2).
- [ ] Confirm task-role secret + S3 perms (step 3).
- [ ] **Test on 2–3 accessions** (a small seedfile) end-to-end before the 624-accession run.
- [ ] **Minimizations** (v0.2): restrict `samtools depth` to the RT-neighborhood window (~18× smaller depth.tsv.gz); drop `features.json` + `calmd`/`mapped.md.bam` from the core; single-aligner only (already via `FDZ004_SHORT_ALIGNERS=bowtie2`).
- [ ] Bulk **genes-cache** prebuild (per-accession GFF fetch was 34 s) → set `FDZ004_GENES_CACHE`.
- [ ] Aggregate `igr_scores.*.tsv` across accessions + the AUC round on updated controls.
