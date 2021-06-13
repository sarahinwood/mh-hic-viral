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
samblaster_container = 'shub://TomHarrop/align-utils:samblaster_0.1.24'

#########
# RULES #
#########

rule target:
    input:
    	##sampe output:
    	'output/mh_hic_matlock/matlock_viral_scaffolds.csv',
    	'output/hic_genome_aln/interaction_matrix.csv',
    	'output/hic_genome_aln/interaction_locations.csv'

##as per phase genomics reccomendations - https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html

#################
## look at res ##
#################

rule hic_genome_interaction_matrix:
	input:
		bam = 'output/hic_genome_aln/matlock_links.out'
	output:
		interaction_matrix = 'output/hic_genome_aln/interaction_matrix.csv',
		interaction_locations = 'output/hic_genome_aln/interaction_locations.csv'
	log:
		'output/logs/hic_genome_interaction_matrix.log'
	singularity:
		tidyverse_container
	threads:
		10
	script:
		'src/hic_genome_interaction_matrix.R'

#################################
## generate interaction matrix ##
#################################

rule matlock_juicer_hic:
	input:
		'output/hic_genome_aln/hic_align.bam'
	output:
		'output/hic_genome_aln/matlock_links.out'
	log:
		'output/logs/matlock_juicer_hic.log'
	singularity:
		matlock_container
	threads:
		10
	shell:
		'matlock bam2 juicer '
		'{input} '
		'{output} '
		'2> {log}'

################
## create bam ##
################

##-F 2316 excludes reads which were unmapped or whose mate was unmapped,
##as well as supplementary & not primary alignments
rule samtools_hic_align:
	input:
		'output/hic_genome_aln/samblaster.sam'
	output:
		'output/hic_genome_aln/hic_align.bam'
	threads:
		10
	log:
		'output/logs/samtools_hic_align.log'
	shell:
		'samtools view '
		'-S -h -b -F 2316 '
		'{input} '
		'> {output} '
		'2> {log}'

####################################
## flag and remove PCR duplicates ##
####################################

rule samblaster_hic_align:
	input:
		'output/hic_genome_aln/bwa_mem.sam'
	output:
		temp('output/hic_genome_aln/samblaster.sam')
	log:
		'output/logs/samblaster.log'
	singularity:
		samblaster_container
	shell:
		'samblaster '
		'-i {input} '
		'-o {output} '
		'2> {log}'

#####################################
## align hi-c reads to hi-c genome ##
#####################################

	input:
		hic_genome = 'data/mh_hic_genome/Mh_Hi-C_PGA_assembly.fasta',
		hic_r1 = 'data/hic_reads/microctonus_hyperodae_harrop_S3HiC_R1.fastq.gz',
		hic_r2 = 'data/hic_reads/microctonus_hyperodae_harrop_S3HiC_R2.fastq.gz'
	output:
		temp('output/hic_genome_aln/bwa_mem.sam')
	log:
		'output/logs/bwa_mem.log'
	params:
		index_dir = 'output/mh_hic_index/mh_hic_genome'
	threads:
		20
	shell:
		'bin/bwa/bwa mem -5SP -t '
		'{threads} {params.index_dir} '
		'{input.hic_r1} {input.hic_r2} '
		'> {output} '
		'2> {log}'

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