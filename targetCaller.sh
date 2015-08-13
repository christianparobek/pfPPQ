## A tool to variant call in potential PPQ-associated genes
## Started 13 August 2015
## Christian parobek


##########################################################################
############################# DEFINE PATHS ###############################
##########################################################################

ref=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta
picard=/nas02/apps/picard-1.88/picard-tools-1.88
gatk=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar


##########################################################################
########################### CALL HAPLOTYPES ##############################
##########################################################################

## TRYING OUT HAPLOTYPE CALLER
java -jar $gatk \
	-T HaplotypeCaller \
	-R $ref \
	-L targets.intervals \
	-I our_goods_5x@60%.list \
	-o variants/targets.vcf
		# gatk.intervals includes just the chromosomes and mitochondria
