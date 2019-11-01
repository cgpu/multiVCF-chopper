#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "                      SelectVariants: SNPs && PASS                  "
log.info "===================================================================="

// set threadmem equal to total memory divided by number of threads
int threads = Runtime.getRuntime().availableProcessors()
threadmem = (((Runtime.getRuntime().maxMemory() * 4) / threads) as nextflow.util.MemoryUnit)

// fasta
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_select_variants_PASS; fasta_vcf2maf }
}

// fai
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fasta index file not found: ${params.fai}" }
           .into { fai_select_variants_PASS ; fai_vcf2maf  }
}

// dict
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { dict_select_variants_PASS ; dict_vcf2maf }
}

Channel
    .fromPath("${params.inputdir}/*.vcf")
    .set {  vcf_filtered_for_select_variants}

Channel
    .fromPath("${params.inputdir}/*.idx")
    .set {  idx_vcf_filtered_for_select_variants}

// SelectVariants only the ones with PASS:
// Found here: 
// https://gatkforums.broadinstitute.org/gatk/discussion/1742/using-selectvariants-to-output-pass-records
// Should have done with grep but let's gatk much, shall we
//     -select 'vc.isNotFiltered()' : only the ones with flag PASS
//     -select-type SNP             : only SNPs

process SelectSNPsPASS {

    tag "${filtered_vcf}"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/SelectedSomaticSNPs_VCF", mode: 'copy'

    input:
    file(filtered_vcf) from vcf_filtered_for_select_variants
    file(filtered_vcf_idx) from idx_vcf_filtered_for_select_variants
    each file(fasta) from fasta_select_variants_PASS
    each file(fai) from fai_select_variants_PASS
    each file(dict) from dict_select_variants_PASS

    output:
    set file("*vcf") into vcf_SNPs_PASS_for_vcf2maf, vcf_SNP_count_info_channel
    file("*vcf.idx") into idx_vcf_SNPs_PASS_for_vcf2maf

    script:
    """
    gatk SelectVariants \
    -R ${fasta} \
    -V $filtered_vcf \
    -O ${filtered_vcf.simpleName}.passed.SNPs.vcf \
    -select 'vc.isNotFiltered()' \
    -select-type SNP
   """
}

process Vcf2maf {

    tag "${vcf_passed_SNPs}"
    container 'levim/vcf2maf:1.0'
    publishDir "${params.outdir}/SelectedSomaticSNPs_MAF", mode: 'copy'
    echo true

    input:
    file(vcf_passed_SNPs) from vcf_SNPs_PASS_for_vcf2maf
    file(idx_vcf_passed_SNPs) from idx_vcf_SNPs_PASS_for_vcf2maf
    each file(fasta) from fasta_vcf2maf
    each file(fai) from fai_vcf2maf
    each file(dict) from dict_vcf2maf

    output:
    file("*") into vcf2maf_annotated_files_channel
    file("*.maf") into maf_SNP_count_info_channel
 
    script:
    """
    basename=\$(echo ${vcf_passed_SNPs.simpleName})
    tumourID=\$(echo \$basename | cut -f 1 -d '_')
    normalID=\$(echo \$basename | cut -f 4 -d '_')

    echo  \$tumourID
    echo  \$normalID

    perl /opt/vcf2maf/vcf2maf.pl \
    --input-vcf $vcf_passed_SNPs \
    --output-maf "\${basename}.passed.SNPs.maf"  \
    --tumor-id \${tumourID} \
    --normal-id \${normalID} \
    --ref-fasta /vepdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --ncbi-build  GRCh37 \
    --filter-vcf /vepdata/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    --vep-path /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor \
    --vep-data /vepdata/ \
    --vep-forks 2 \
    --buffer-size 200 \
    --species homo_sapiens \
    --cache-version 89
    """
}

process CountSNPs {

    tag "Counting.."
    container 'levim/vcf2maf:1.0'
    publishDir "${params.outdir}/SummaryInfo", mode: 'copy'
    echo true

    input:
    file (vcf) from vcf_SNP_count_info_channel.collect().ifEmpty([])
    file (maf) from maf_SNP_count_info_channel.collect().ifEmpty([])

    output:
    file("*") into counts_of_snps_channel
 
    script:
    """
    grep -vcw '^#' * >> lines_without_comments_per_file.txt
    grep -cw '^#' * >> lines_with_comments_per_file.txt
    wc -l * >> all_lines_per_file.txt
    """
}