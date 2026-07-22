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
    // Publish ONLY the minimal artifact set (allow-list) — drops the regenerable/heavy scratch
    // (bt2 index, mapped.bam, reference.fasta, features.json, aligner_comparison.json, md.bam) from S3.
    publishDir "${params.output_path}/${params.project}", mode: 'copy', overwrite: true,
        saveAs: { f ->
            def b = f.tokenize('/').last()
            (b in ['regions.bed', 'genes.tsv', 'depth.tsv.gz', 'coverage_summary.json',
                   'state.json', 'run.log', 'flagstat.tsv', 'idxstats.tsv']
             || b.startsWith('igr_scores.') || b.startsWith('run_summary.')) ? f : null
        }

    input:
    tuple val(sra), path(loci)

    output:
    path "runs/**", emit: results, optional: true

    script:
    """
    set -euo pipefail

    # region for ALL aws calls in this task (per Xiandong; matches nf-basespace bin/aws_secretsmanager*.sh) — also fixes the SRA aws-s3 pull
    export AWS_DEFAULT_REGION=${params.aws_region}
    # ERDA atlas bearer URL from Secrets Manager — plain-string secret, read directly (no jq); \$() strips the trailing newline
    export RT_ATLAS_ERDA_BASE=\$(aws secretsmanager get-secret-value --secret-id ${params.erda_secret_id} --query SecretString --output text)
    export FDZ004_SHORT_ALIGNERS='${params.short_aligner}' FDZ004_THREADS=${task.cpus} FDZ004_MAX_SRA_GB=${params.max_sra_gb}
    export FDZ004_IGR_MIN=${params.igr_min} FDZ004_IGR_MAX=${params.igr_max} FDZ004_CLEANUP=1
    export FDZ004_MINIMAL=1   # depth over the RT-neighborhood window only; skip calmd/features/aligner-comparison

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
