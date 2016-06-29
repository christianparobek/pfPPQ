## A Snakefile
## Started 27 June 2016
## To make by-population AFS plots for CP1-4


##########################################################################################

workdir: '/proj/julianog/users/ChristianP/pfPPQ/'
REF = '/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta'
data='/proj/julianog/users/ChristianP/pfPPQ/dadi/data'
orig_vcf='/proj/julianog/users/ChristianP/pfPPQ/variants/slimmed.recode.84fix.vcf'
#readWD = '/proj/julianog/users/ChristianP/cambodiaWGS/pf/'
#DATEDSAMPS, = glob_wildcards('/proj/julianog/users/ChristianP/cambodiaWGS/pf/symlinks/{ds}_R1.fastq.gz')
#SAMPLES, = glob_wildcards('/proj/julianog/users/ChristianP/cambodiaWGS/pf/aln/{sample}.merged.bam')

##########################################################################################


####### Target #######
rule all:
	#input: 'dadi/afsPlots/vcfs/cp1.recode.vcf'
	#input: 'dadi/afsPlots/vcfs/cp2.recode.vcf'
	#input: 'dadi/afsPlots/vcfs/cp3.recode.vcf'
	#input: 'dadi/afsPlots/vcfs/cp4.recode.vcf'
	#input: expand('dadi/afsPlots/data/{cp}/dadi_input/{cp}.syn.dadi', cp = 'cp1 cp2 cp3 cp4 slimmed84fix'.split())
	input: expand('dadi/afsPlots/data/{cp}/dadi_output_afs/{cp}.afs', cp = 'cp1 cp2 cp3 cp4 slimmed84fix'.split())


#rule plot_afs:
#	input: pf_all = 'dadi/data/pf_all/dadi_output_afs/pf_all.afs', pf_cp2 = 'dadi/data/pf_cp2/dadi_output_afs/pf_cp2.afs', pv_all = 'dadi/data/pv_all/dadi_output_afs/pv_all.afs', pv_mono = 'dadi/data/pv_mono/dadi_output_afs/pv_mono.afs'
#	output: full_datasets = 'dadi/data/pfall_pvall.svg', subset_datasets = 'dadi/data/pfcp2_pvmono.svg'
#	shell: 'Rscript dadi/afsPlotter.r {input.pf_all} {input.pf_cp2} {input.pv_all} {input.pv_mono} {output.full_datasets} {output.subset_datasets}'

rule dadi_calc_afs:
	input: annovcf = 'dadi/afsPlots/data/{cp}/vcfs/{cp}.anno.vcf', syn = 'dadi/afsPlots/data/{cp}/dadi_input/{cp}.syn.dadi', nonsyn = 'dadi/afsPlots/data/{cp}/dadi_input/{cp}.nonsyn.dadi', genic = 'dadi/afsPlots/data/{cp}/dadi_input/{cp}.genic.dadi', intergenic = 'dadi/afsPlots/data/{cp}/dadi_input/{cp}.intergenic.dadi'
	output: 'dadi/afsPlots/data/{cp}/dadi_output_afs/{cp}.afs'
	shell: 'bash dadi/afsPlots/afsPrinter.sh {input.annovcf} {input.syn} {input.nonsyn} {input.genic} {input.intergenic} {output}'

rule dadi_format:
	input: vcf = 'dadi/afsPlots/vcfs/{cp}.recode.vcf'
	output: syn = 'dadi/afsPlots/data/{cp}/dadi_input/{cp}.syn.dadi', nonsyn = 'dadi/afsPlots/data/{cp}/dadi_input/{cp}.nonsyn.dadi', genic = 'dadi/afsPlots/data/{cp}/dadi_input/{cp}.genic.dadi', intergenic = 'dadi/afsPlots/data/{cp}/dadi_input/{cp}.intergenic.dadi'
	params: "{cp}"
	shell: 'bash dadi/afsPlots/vcf2dadi.sh {input.vcf} {params} pf'


######################
## MAKE SUBSET VCFS ##
######################

rule make_cp1_vcf :
	input: group = 'dadi/data/cp_groups/cp1.txt'
	output: 'dadi/afsPlots/vcfs/cp1.recode.vcf'
	shell: 'vcftools --vcf {orig_vcf} --keep {input.group} --non-ref-af 0.01 --recode --out dadi/afsPlots/vcfs/cp1'
		# the maf flag excludes all the ghost variants after we filter out CP1, CP3, CP4
		# but it turns out we want to use --non-ref-af instead, because --maf eliminates all fixed non ref variants

rule make_cp2_vcf :
	input: group = 'dadi/data/cp_groups/cp2.txt'
	output: 'dadi/afsPlots/vcfs/cp2.recode.vcf'
	shell: 'vcftools --vcf {orig_vcf} --keep {input.group} --non-ref-af 0.01 --recode --out dadi/afsPlots/vcfs/cp2'

rule make_cp3_vcf :
	input: group = 'dadi/data/cp_groups/cp3.txt'
	output: 'dadi/afsPlots/vcfs/cp3.recode.vcf'
	shell: 'vcftools --vcf {orig_vcf} --keep {input.group} --non-ref-af 0.01 --recode --out dadi/afsPlots/vcfs/cp3'

rule make_cp4_vcf :
	input: group = 'dadi/data/cp_groups/cp4.txt'
	output: 'dadi/afsPlots/vcfs/cp4.recode.vcf'
	shell: 'vcftools --vcf {orig_vcf} --keep {input.group} --non-ref-af 0.01 --recode --out dadi/afsPlots/vcfs/cp4'
