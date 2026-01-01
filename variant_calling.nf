#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads       = "${baseDir}/01_Data/01_Raw_Reads/*_{1,2}.fastq.gz"
params.outdir      = "${baseDir}/02_Results"

/*********************************
 * Channels
 *********************************/
paired_reads = Channel.fromFilePairs(params.reads, size: 2, checkIfExists: true)
                  .map { id, files -> tuple(id, files[0], files[1]) }

/*********************************
 * QC
 *********************************/
process fastqc {
    publishDir "${params.outdir}/01_Reads_QC/01_Raw_fastqc", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv7'
    input:
      tuple val(id), path(r1), path(r2)
    output:
      path "*_fastqc.*"
    script:
    """
    fastqc $r1 $r2 --outdir=.
    """
}

process multiqc_raw_reads {
    publishDir "${params.outdir}/01_Reads_QC/02_MultiQC_Raw_reads", mode: 'copy'
    container 'multiqc/multiqc:latest'
    input:
      path fastqc_files
    output:
      path "multiqc_raw_reads_report.html"
    script:
    """
    multiqc . -n multiqc_raw_reads_report.html
    """
}

/*********************************
 * Trimming
 *********************************/
process trimmomatic {
    publishDir "${params.outdir}/02_Trimmed_reads", mode: 'copy'
    container 'staphb/trimmomatic:latest'
    input:
      tuple val(id), path(r1), path(r2)
    output:
      tuple val(id), path("${id}_trimmed_1.fastq.gz"), path("${id}_trimmed_2.fastq.gz")
    script:
    """
    trimmomatic PE -threads $task.cpus \
      $r1 $r2 \
      ${id}_trimmed_1.fastq.gz ${id}_forward_unpaired.fastq.gz \
      ${id}_trimmed_2.fastq.gz ${id}_reverse_unpaired.fastq.gz \
      LEADING:30 TRAILING:30 MINLEN:60
    """
}

process fastqc_trimmed {
    publishDir "${params.outdir}/01_Reads_QC/03_Trimmed_fastqc", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv7'
    input:
      tuple val(id), path(r1), path(r2)
    output:
      path "*_fastqc.*"
    script:
    """
    fastqc $r1 $r2 --outdir=.
    """
}

process multiqc_trim_reads {
    publishDir "${params.outdir}/01_Reads_QC/04_MultiQC_Trimmed_reads", mode: 'copy'
    container 'multiqc/multiqc:latest'
    input:
      path fastqc_files
    output:
      path "multiqc_trimmed_reads_report.html"
    script:
    """
    multiqc . -n multiqc_trimmed_reads_report.html
    """
}

/*********************************
 * Alignment
 *********************************/
process BWA_MEM {
    tag { id }
    container 'pegi3s/bwa:0.7.17'

    publishDir "${params.outdir}/03_Alignment/01_BWA_MEM", mode: 'copy'

    input:
      tuple val(id), path(r1), path(r2)

    output:
      tuple val(id), path("${id}.bam")

    script:
    """
    bwa mem -t ${task.cpus} -R "@RG\\\\tID:${id}\\\\tSM:${id}\\\\tPL:ILLUMINA" /mnt/02_Ref_genome/hg38.fa $r1 $r2 > ${id}.bam
    """
}

process SORT_BAM {
    tag { id }
    container 'staphb/samtools:latest'
    publishDir "${params.outdir}/03_Alignment/02_Sorted", mode: 'copy'
    input:
      tuple val(id), path(bam)
    output:
      tuple val(id), path("${id}.sorted.bam")
    script:
    """
    samtools sort -o ${id}.sorted.bam $bam
    """
}

process MARK_DUPLICATES {
    tag { id }
    container 'broadinstitute/picard:latest'
    publishDir "${params.outdir}/03_Alignment/03_Mark_Dup", mode: 'copy'
    input:
      tuple val(id), path(bam)
    output:
      tuple val(id), path("${id}.md.bam"), path("${id}.metrics.txt")
    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates I=$bam O=${id}.md.bam M=${id}.metrics.txt CREATE_INDEX=true
    """
}

/*********************************
 * Variant Calling
 *********************************/
process BQSR {
    publishDir "${params.outdir}/03_Alignment/04_BQSR", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.bqsr.bam"), path("${id}.bqsr.bai"), path("recal.table")

    script:
    """
    gatk BaseRecalibrator -R /mnt/02_Ref_genome/hg38.fa -I $bam --known-sites /mnt/02_Ref_genome/Homo_sapiens_assembly38.dbsnp138.vcf -O recal.table --create-output-bam-index true
    gatk ApplyBQSR -R /mnt/02_Ref_genome/hg38.fa -I $bam --bqsr-recal-file recal.table -O ${id}.bqsr.bam --create-output-bam-index true
    """
}

process HAPLOTYPE_CALLER {
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/04_Variants/01_Raw", mode: 'copy'
    input:
      tuple val(id), path(bam), path(ref)
    output:
      tuple val(id), path("${id}_raw.vcf")
    script:
    """
    gatk HaplotypeCaller --native-pair-hmm-threads $task.cpus -R /mnt/02_Ref_genome/hg38.fa -I $bam -O ${id}_raw.vcf
    """
}

process SPLIT_VARIANTS {

    publishDir "${params.outdir}/04_Variants/02_Split", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_raw_snps.vcf"), path("${id}_raw_indels.vcf")

    script:
    """
    gatk SelectVariants -R /mnt/02_Ref_genome/hg38.fa -V $vcf --select-type SNP   -O ${id}_raw_snps.vcf
    gatk SelectVariants -R /mnt/02_Ref_genome/hg38.fa -V $vcf --select-type INDEL -O ${id}_raw_indels.vcf
    """
}

/*********************************
 * Filter Variants
 *********************************/

process FILTER_SNPS {
    
    publishDir "${params.outdir}/04_Variants/03_Filtered", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_filtered_snps.vcf")

    script:
    """
    gatk VariantFiltration -R /mnt/02_Ref_genome/hg38.fa -V $vcf -O ${id}_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
        -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"
    """
}

process FILTER_INDELS {
    
    publishDir "${params.outdir}/04_Variants/03_Filtered", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_filtered_indels.vcf")

    script:
    """
    gatk VariantFiltration -R /mnt/02_Ref_genome/hg38.fa -V $vcf -O ${id}_filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" \
        -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"
    """
}

process SELECT_PASS_SNPS {
    
    publishDir "${params.outdir}/04_Variants/04_Pass", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_snps_pass.vcf")

    script:
    """
    gatk SelectVariants --exclude-filtered -V $vcf -O ${id}_snps_pass.vcf
    """
}

process SELECT_PASS_INDELS {
    
    publishDir "${params.outdir}/04_Variants/04_Pass", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_indels_pass.vcf")

    script:
    """
    gatk SelectVariants --exclude-filtered -V $vcf -O ${id}_indels_pass.vcf
    """
}

process REMOVE_GT_FILTERS_SNPS {
    
    publishDir "${params.outdir}/04_Variants/05_GT_Filtered", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_snps_gt_clean.vcf")

    script:
    """
    grep -v -E "DP_filter|GQ_filter" $vcf > ${id}_snps_gt_clean.vcf
    """
}

process REMOVE_GT_FILTERS_INDELS {
    
    publishDir "${params.outdir}/04_Variants/05_GT_Filtered", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_indels_gt_clean.vcf")

    script:
    """
    grep -v -E "DP_filter|GQ_filter" $vcf > ${id}_indels_gt_clean.vcf
    """
}

/*********************************
 * FUNCOTATOR
 *********************************/

process FUNCOTATOR_SNPS {
    
    publishDir "${params.outdir}/04_Variants/06_Funcotator", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_snps_funcotated.vcf")

    script:
    """
    gatk Funcotator \
        --variant $vcf \
        --reference /mnt/02_Ref_genome/hg38.fa \
        --ref-version hg38 \
        --data-sources-path /mnt/03_Funcotator/funcotator_dataSources.v1.8.hg38.20230908g \
        --output ${id}_snps_funcotated.vcf \
        --output-file-format VCF
    """
}

process FUNCOTATOR_INDELS {
    
    publishDir "${params.outdir}/04_Variants/06_Funcotator", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_indels_funcotated.vcf")

    script:
    """
    gatk Funcotator \
        --variant $vcf \
        --reference /mnt/02_Ref_genome/hg38.fa \
        --ref-version hg38 \
        --data-sources-path /mnt/03_Funcotator/funcotator_dataSources.v1.8.hg38.20230908g \
        --output ${id}_indels_funcotated.vcf \
        --output-file-format VCF
    """
}

process VARIANTS_TO_TABLE_SNPS {
    
    publishDir "${params.outdir}/04_Variants/07_Tables", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    path "${id}_snps.table"

    script:
    """
    gatk VariantsToTable -V $vcf -F AC -F AN -F DP -F AF -F FUNCOTATION -O ${id}_snps.table
    """
}

process VARIANTS_TO_TABLE_INDELS {
    
    publishDir "${params.outdir}/04_Variants/07_Tables", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(id), path(vcf)

    output:
    path "${id}_indels.table"

    script:
    """
    gatk VariantsToTable -V $vcf -F AC -F AN -F DP -F AF -F FUNCOTATION -O ${id}_indels.table
    """
}

/*********************************
 * END POINT
 *********************************/

process end_point {

    input:
    path raw_report
    path trim_report
    path snps_tables
    path indels_tables

    script:
    """
    echo "Variants tables generation completed successfully."
    echo "Raw MultiQC report(s):"
    ls -lh ${raw_report}

    echo "Trimmed MultiQC report(s):"
    ls -lh ${trim_report}

    echo "SNPs tables:"
    ls -lh ${snps_tables}

    echo "INDELs tables:"
    ls -lh ${indels_tables}
    """
}

/*********************************
 * Workflow
 *********************************/
workflow {
    
    // QC and trimming steps
    fastqc_out = fastqc(paired_reads)
    multiqc_raw_report = multiqc_raw_reads(fastqc_out)
    trimmed_reads = trimmomatic(paired_reads)
    fastqc_trim_out = fastqc_trimmed(trimmed_reads)
    multiqc_trim_report = multiqc_trim_reads(fastqc_trim_out)

    // Alignment
    aligned = BWA_MEM(trimmed_reads)
    sorted  = SORT_BAM(aligned)
    marked  = MARK_DUPLICATES(sorted)

    // Base recalibration and variant calling
    bqsr    = BQSR(marked)
    raw_vcf = HAPLOTYPE_CALLER(bqsr)
    split   = SPLIT_VARIANTS(raw_vcf)

    // Filter variants
    snps_f  = FILTER_SNPS(split.map { id,s,i -> tuple(id,s) })
    ind_f   = FILTER_INDELS(split.map { id,s,i -> tuple(id,i) })

    snps_p  = SELECT_PASS_SNPS(snps_f)
    ind_p   = SELECT_PASS_INDELS(ind_f)

    snps_c  = REMOVE_GT_FILTERS_SNPS(snps_p)
    ind_c   = REMOVE_GT_FILTERS_INDELS(ind_p)

    // FUNCOTATOR
    snps_a  = FUNCOTATOR_SNPS(snps_c)
    ind_a   = FUNCOTATOR_INDELS(ind_c)

    // Generate final variant tables
    snps_t = VARIANTS_TO_TABLE_SNPS(snps_a)
    ind_t  = VARIANTS_TO_TABLE_INDELS(ind_a)

    // END_POINT process takes the SNP and INDEL tables
    end_point(multiqc_raw_report, multiqc_trim_report, snps_t, ind_t)
}