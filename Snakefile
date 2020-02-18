#!/usr/bin/env python3

import pathlib2

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

#########
# RULES #
#########

rule target:
    input:
    	expand('output/lbfv_aln/lbfv_aln_{direction}.bam', direction=['R1', 'R2']),
    	expand('output/lbfv_aln/view_bam_{direction}.out', direction=['R1', 'R2']),
    	expand('output/lbfv_aln/lbfv_aln_{direction}.bam.bai', direction=['R1', 'R2'])

rule view_bam:
	input:
		bam = 'output/lbfv_aln/lbfv_aln_{direction}.bam'
	output:
		'output/lbfv_aln/view_bam_{direction}.out'
	params:
		bam = lambda wildcards, input: resolve_path(input.bam)
	log:
		'output/logs/view_bam_{direction}.log'
	threads:
		10
	shell:
		'samtools view '
		' -c -F 260 '
		'{params.bam} '
		'> {output} '
		'2> {log}'

rule index_bam:
	input:
		bam = 'output/lbfv_aln/lbfv_aln_{direction}.bam'
	output:
		index = 'output/lbfv_aln/lbfv_aln_{direction}.bam.bai'
	params:
		bam = lambda wildcards, input: resolve_path(input.bam)
	log:
		'output/logs/index_bam_{direction}.log'
	shell:
		'samtools index '
		'{params.bam} '
		'> {output.index} '
		'2> {log}'

rule sam_to_bam:
	input:
		sam = 'output/lbfv_aln/lbfv_aln_{direction}.sam'
	output:
		bam = 'output/lbfv_aln/lbfv_aln_{direction}.bam'
	params:
		sam = lambda wildcards, input: resolve_path(input.sam)
	log:
		'output/logs/sam_to_bam_{direction}.log'
	shell:
		'samtools view '
		'-bS '
		'{params.sam} '
		'> {output.bam} '
		'2> {log}'

##running samse for now as wouldn't expect read pairs to map at all
rule bwa_samse:
	input:
		sai = 'output/lbfv_aln/HiC_{direction}.sai',
		hic_reads = 'data/microctonus_hyperodae_harrop_S3HiC_{direction}.fastq.gz'
	output:
		temp('output/lbfv_aln/lbfv_aln_{direction}.sam')
	params:
		index_dir = 'output/lbfv_index/lbfv_genome',
		sai = lambda wildcards, input: resolve_path(input.sai),
		hic_reads = lambda wildcards, input: resolve_path(input.hic_reads)
	threads:
		10
	log:
		'output/logs/bwa_samse_{direction}.log'
	shell:
		'bwa samse '
		'{params.index_dir} '
		'{params.sai} '
		'{params.hic_reads} '
		'> {output} '
		'2> {log}'

rule bwa_aln:
	input:
		hic_reads = 'data/microctonus_hyperodae_harrop_S3HiC_{direction}.fastq.gz'
	output:
		'output/lbfv_aln/HiC_{direction}.sai'
	params:
		index_dir = 'output/lbfv_index/lbfv_genome',
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
		viral_genome = 'data/lbfv_genome.fasta'
	output:
		lbfv_index = 'output/lbfv_index/lbfv_genome.bwt'
	params:
		index_dir = 'output/lbfv_index/lbfv_genome'
	log:
		'output/logs/bwa_index.log'
	shell:
		'bwa index '
		'{input.viral_genome} '
		'-p {params.index_dir} '
		'2> {log}'
