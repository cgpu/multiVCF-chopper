#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "                 Select samples from multisample VCF                 "
log.info "===================================================================="

// set threadmem equal to total memory divided by number of threads
int threads = Runtime.getRuntime().availableProcessors()
threadmem = (((Runtime.getRuntime().maxMemory() * 4) / threads) as nextflow.util.MemoryUnit)

// fasta
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { ch_fasta; fasta_vcf2maf }
}

// fai
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fasta index file not found: ${params.fai}" }
           .into { ch_fai ; fai_vcf2maf  }
}

// dict
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { ch_dict ; dict_vcf2maf }
}

Channel
    .fromPath(params.multiVCF)
    .toSortedList
    .set {  ch_multiVCF}

Channel
    .fromPath(params.multiVCF_index)
    .toSortedList
    .set {  ch_multiVCF_tbi}

Channel
    .fromPath(params.list_folder)
    .set {  ch_subset_lists}

process chop_multiVCF {

    tag "${filtered_vcf}"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/subsampled_multisample_vcf/{$sample_list.simpleName}", mode: 'copy'

    input:
    each file(vcf) from ch_multiVCF
    each file(vcf_index) from ch_multiVCF_tbi
    each file(fasta) from ch_fasta
    each file(fai) from ch_fai
    each file(dict) from ch_dict
    file(sample_list) from ch_subset_lists

    output:
    file("*") into ch_nowhere

    script:
    """
    gatk SelectVariants \
    -R ${fasta} \
    -V $vcf \
    -O ${filtered_vcf.simpleName}.passed.SNPs.vcf \
    --sample_file ${sample_list} \
   """
}