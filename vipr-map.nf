#!/usr/bin/env nextflow
println "   ViPR3 (bwa)  "
println "===================="

params.input = "reads" //Input folder
params.patterns = "*_R{1,2}.fastq.gz" //Patterns of reads
params.read_pairs = params.input + "/" + params.patterns
params.output = "results"
params.ref = "FJ639679.fasta"
params.threads = 4
params.memory = 8
params.max_depth = 1000
params.qc = false

threads = params.threads
memory = params.memory
ref = file(params.ref)
max_depth = params.max_depth

PLOTPY="$PWD/aux/plot.py"

params.help = false
// Display help menu
if(params.help) {
    log.info ''
    log.info '¯\\_(ツ)_/¯ ViPR3 Pipeline (bwa) ¯\\_(ツ)_/¯'
    log.info ''
    log.info ''
    log.info 'Notes: ••••Using docker, check nextflow.config••••'
    log.info 'Notes: ••••This pipeline aims to generate genome/consensus for viral NGS sequence••••'
    log.info 'Notes: ••••Ported from ViPR snakemake pipelines by Andreas Wilm••••'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow vip-spades.nf [options]'
    log.info ''
    log.info 'Required params: '
    log.info '    --input           DIR     Directory of paired FASTQ files'
    log.info '    --patterns        FILE    Pattern for scanning the reads [*_R{1,2}.fastq.gz]'
    log.info '    --ref             FILE    Fasta file of the reference genome'
    log.info '    --threads         INT     Number of threads to use for each process [4]'
    log.info '    --memory          INT     Memory allocated for each process [8]'
    log.info '    --output          DIR     Directory to write output files to [results]'
    log.info ''
    log.info ''
    return
}

Channel
    .fromFilePairs(params.read_pairs, flat: true)
    .ifEmpty { exit 1, "Read pairs could not be found: ${params.read_pairs}" }
    .set { read_pairs }

if (params.qc) {
    process reads_qc {
    publishDir "$params.output/qc", mode: "copy"

    tag { read_id }

    input:
    set read_id, file(forward), file(reverse) from read_pairs

    output:
    set read_id, file("${read_id}_R1.trimmed.fastq.gz"), file("${read_id}_R2.trimmed.fastq.gz") into (qc_tadpole, ch_map)

    """
    fastp -c \
    -i $forward \
    -I $reverse \
    -o ${read_id}_R1.trimmed.fastq.gz \
    -O ${read_id}_R2.trimmed.fastq.gz \
    -h ${read_id}_fastp.html
    """
    }
} else {
    read_pairs.into{ qc_tadpole; ch_map }
}



process map_align {

    tag { read_id }

    input:
    set read_id, file(forward), file(reverse) from ch_map
    file(ref)
    output:
    set read_id, file("${read_id}.bam") into ch_map_align

    """
    bwa index $ref

    bwa mem -M -t $threads \
    -R '@RG\tID:${read_id}\tSM:null\tLB:null\tCN:null' \
    $ref $forward $reverse | lofreq viterbi -f $ref - | samtools fixmate - - | \
    lofreq alnqual -u - $ref | \
    lofreq indelqual --dindel -f $ref - | \
    samtools sort -@ ${threads} -o ${read_id}.bam -T ${read_id}.bam.tmp >& map.log

    """
}

process lofreq_call {
    publishDir "$params.output/$read_id/map", mode: "copy"

    tag { read_id }

    input:
    file(ref)
    set read_id, file("${read_id}.bam") from ch_map_align

    output:
    file("lofreq_call.log")
    file("${read_id}.bam") into ch_coverage
    set read_id, file("${read_id}.vcf.gz") into vcf_lofreq

    """
    samtools faidx $ref

    samtools index ${read_id}.bam

    lofreq call-parallel --pp-threads $threads --call-indels \
    -f $ref -o ${read_id}.vcf -d $max_depth ${read_id}.bam >& lofreq_call.log

    gzip ${read_id}.vcf
    """
}

process coverage {
    publishDir "$params.output/$read_id/map", mode: "copy"

    tag { read_id }

    input:
    set read_id, file("${read_id}.vcf.gz") from vcf_lofreq
    file("${read_id}.bam") from ch_coverage

    output:
    file("${read_id}.png")

    """
    bedtools genomecov -d -ibam ${read_id}.bam | gzip > ${read_id}.cov.gz

    $PLOTPY --vcf ${read_id}.vcf.gz --cov ${read_id}.cov.gz --plot ${read_id}.png
    """
}

workflow.onComplete {
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command Line:     $workflow.commandLine"
    log.info "Container:        $workflow.container"
    log.info "Duration:     $workflow.duration"
    log.info "Output Directory: $params.output"
}
