#!/usr/bin/env python3

import pathlib2

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

##############
##CONTAINERS##
##############

matlock_container = 'shub://TomHarrop/seq-utils:matlock_9fe3fdd'
tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'

#########
# RULES #
#########

rule target:
    input:
    	##sampe output:
    	'output/mh_hic_matlock/matlock_viral_scaffolds.csv'

rule filter_matlock_viral_scaffolds:
	input:
		matlock_output = 'output/mh_hic_matlock/matlock_links.out',
		viral_scaffold_list = 'data/mh_genome/viral_scaffold_ids.txt'
	output:
		viral_matlock = 'output/mh_hic_matlock/matlock_viral_scaffolds.csv'
	singularity:
		tidyverse_container
	threads:
		10
	log:
		'output/logs/filter_matlock_viral.log'
	script:
		'src/filter_matlock_viral.R'

rule matlock_juicer:
	input:
		bam = 'data/mh_hic_genome/microctonus_hyperodae_harrop.bam'
	output:
		matlock_juicer = 'output/mh_hic_matlock/matlock_links.out'
	log:
		'output/logs/matlock_juicer.log'
	singularity:
		matlock_container
	threads:
		10
	shell:
		'matlock bam2 juicer '
		'{input.bam} '
		'{output.matlock_juicer} '
		'2> {log}'

########

rule view_bam:
	input:
		bam = 'output/viral_scaffolds_aln/sampe/sorted.bam'
	output:
		'output/viral_scaffolds_aln/sampe/view_bam.out'
	log:
		'output/logs/view_bam.log'
	threads:
		10
	shell:
		'samtools view '
		' -c -F 260 '
		'{input.bam} '
		'> {output} '
		'2> {log}'

rule index_bam:
	input:
		sorted_bam = 'output/viral_scaffolds_aln/sampe/sorted.bam'
	output:
		index = 'output/viral_scaffolds_aln/sampe/sorted.bam.bai'
	log:
		'output/logs/index_bam.log'
	shell:
		'samtools index '
		'{input.sorted_bam} '
		'> {output.index} '
		'2> {log}'

rule samtools_sort:
    input:
        bam = 'output/mh_hic_aln/mh_hic_aln_sampe.bam'
    output:
        sorted_bam = 'output/mh_hic_aln/sorted.bam'
    log:
        'output/logs/samtools_sort.log'
    threads:
        10
    shell:
        'samtools sort -n '
        '{input.bam} '
        '-f {output.sorted_bam} '
        '2> {log}'

rule sam_to_bam:
	input:
		sampe_sam = 'output/mh_hic_aln/mh_hic_aln_sampe.sam'
	output:
		bam = 'output/mh_hic_aln/mh_hic_aln_sampe.bam'
	log:
		'output/logs/sam_to_bam.log'
	shell:
		'samtools view '
		'-bS '
		'{input.sampe_sam} '
		'> {output.bam} '
		'2> {log}'

##run sampe to determine optimal placement of each read pair (combines reads from both fastq files)
rule bwa_sampe:
	input:
		sai_r1 = 'output/mh_hic_aln/HiC_R1.sai',
		sai_r2 = 'output/mh_hic_aln/HiC_R2.sai',
		hic_r1 = 'data/hic_reads/microctonus_hyperodae_harrop_S3HiC_R1.fastq.gz',
		hic_r2 = 'data/hic_reads/microctonus_hyperodae_harrop_S3HiC_R2.fastq.gz'
	output:
		temp('output/mh_hic_aln/mh_hic_aln_sampe.sam')
	params:
		index_dir = 'output/mh_hic_index/mh_hic_genome'
	threads:
		10
	log:
		'output/logs/bwa_sampe.log'
	shell:
		'bwa sampe '
		'{params.index_dir} '
		'{input.sai_r1} '
		'{input.sai_r2} '
		'{input.hic_r1} '
		'{input.hic_r2} '
		'> {output} '
		'2> {log}'

##align each fastq file separately
rule bwa_aln:
	input:
		hic_reads = 'data/hic_reads/microctonus_hyperodae_harrop_S3HiC_{direction}.fastq.gz',
		index = 'output/mh_hic_index/mh_hic_genome.bwt'
	output:
		'output/mh_hic_aln/HiC_{direction}.sai'
	params:
		index_dir = 'output/mh_hic_index/mh_hic_genome',
		hic_reads = lambda wildcards, input: resolve_path(input.hic_reads)
	threads:
		10
	log:
		'output/logs/bwa_aln_{direction}.log'
	shell:
		'bwa aln '
		'{params.index_dir} '
		'{params.hic_reads} '
		'-t {threads} '
		'> {output} '
		'2> {log}'

##index draft assembly - If your genome is large (>10 Mb), you will need to use the flag -a bwtsw. This will create a file <assembly>.fasta.bwt.
rule bwa_index:
	input:
		genome = 'data/mh_hic_genome/Mh_Hi-C_PGA_assembly.fasta'
	output:
		index = 'output/mh_hic_index/mh_hic_genome.bwt'
	params:
		index_dir = 'output/mh_hic_index/mh_hic_genome'
	log:
		'output/logs/bwa_index.log'
	shell:
		'bwa index '
		'{input.genome} '
		'-p {params.index_dir} '
		'2> {log}'
