#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * nf-rt-rnaseq — per-accession RT-flanking transcription screen on MAF AWS Batch.
 *
 * Seedfile (CSV, header):  sra_accession,loci
 *   - sra_accession : the SRA run to map
 *   - loci          : path (S3 or local) to THIS accession's loci TSV, in rnaseq_pipeline.py schema
 *                     (role,subclass,common_name,read_type,rt_id,sra_accession,contig_id,rt_start,rt_end,strand)
 *   Build both with bin/make_seedfile.py from a control_loci.tsv (e.g. 03_pilot_scaleup/).
 *
 * One task = one accession (Batch fans out; each task gets its own instance + scratch — the feasibility run
 * showed the pipeline is disk/IO-bound, so per-job isolation beats one shared box). Light outputs publish to
 * ${output_path}/${project}/runs/... ; heavy scratch (.sra/FASTQ/raw BAM) is deleted in-task (FDZ004_CLEANUP=1).
 */

process RUN_ACCESSION {
    tag "${sra}"
    publishDir "${params.output_path}/${params.project}", mode: 'copy', overwrite: true

    input:
    tuple val(sra), path(loci)

    output:
    path "runs/**", emit: results, optional: true

    script:
    """
    set -euo pipefail

    # ERDA atlas bearer URL from Secrets Manager (TASK role needs secretsmanager:GetSecretValue on ${params.erda_secret_id})
    export RT_ATLAS_ERDA_BASE=\$(aws secretsmanager get-secret-value --secret-id ${params.erda_secret_id} --region ${params.aws_region} --query SecretString --output text)
    export FDZ004_SHORT_ALIGNERS='${params.short_aligner}' FDZ004_THREADS=${task.cpus} FDZ004_MAX_SRA_GB=${params.max_sra_gb}
    export FDZ004_IGR_MIN=${params.igr_min} FDZ004_IGR_MAX=${params.igr_max} FDZ004_CLEANUP=1

    # steps 1-6: fetch_sra -> extract_fastq -> build_reference -> map (bowtie2 short / minimap2 long) -> postprocess -> coverage
    # cleanup frees .sra + FASTQ + raw sorted BAM after coverage; keeps reference/ + coverage/ (the light outputs)
    python /pipeline/scripts/rnaseq_pipeline.py --loci "${loci}" --only-sra ${sra} --out-dir runs

    # step 8: score (immediate-flank depth-stat / RT-ORF mean, all 5 stats) — multi-contig-corrected scorer
    python /pipeline/scripts/score_igr.py runs "${loci}" runs/igr_scores.${sra}.tsv || true
    """
}

workflow {
    if( !params.seedfile ) error "Provide --seedfile <CSV with header: sra_accession,loci>"

    channel.fromPath(params.seedfile)
        .splitCsv(header: true)
        .map { row -> tuple(row.sra_accession, file(row.loci)) }
        | RUN_ACCESSION
}
