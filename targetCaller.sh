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
########################### CALL GENOTYPES ###############################
##########################################################################

## TRYING OUT HAPLOTYPE CALLER
for indiv in `cat our_goods_5x@60%.list`
do

java -Xmx48g -jar $gatk \
	-T HaplotypeCaller \
	-R $ref \
	-L targets.intervals \
	-I $indiv \
	-ploidy 1 \
	-o variants/${indiv:51:5}.vcf
		# gatk.intervals includes just the chromosomes and mitochondria

done


##########################################################################
######################### CYTOCHROME B MUTS ##############################
##########################################################################

for indiv in `cat our_goods_5x@60%.list`
do

echo ${indiv:51:5}
grep "M76611" variants/${indiv:51:5}.vcf

done


##########################################################################
############################ PfMDR1 MUTS #################################
##########################################################################

cat *.vcf |\
	grep "Pf3D7_05_v3\|CHROM" |\
	grep -v "##contig" |\
	cut -d$'\t' -f1,2,4,5,10 |\
	grep -v "958103\|958218\|958440\|959583\|960137\|960240\|960436\|960723\|961093\|961149\|961830"
## Comparing my population vcf (UG, in my cambodiaWGS repo) to the
## individual HC-called files
## Catting in all the individual HC vcfs, then grepping out the variants
## identified by UG. Left us with an indel in the Asp/Asn rich region.
## Exported the shell output and cleaned up in Excel/gedit

##########################################################################
########################## MIOTTO BACKBONE ###############################
##########################################################################

for indiv in `cat our_goods_5x@60%.list`
do

echo ${indiv:51:5} | tr '\n' '\t' # remove newline
grep "1956225\|2481070\|405362\|405600\|748395\|490720" variants/${indiv:51:5}.vcf | cut -d$'\t' -f2 | tr '\n' '\t'
echo # add a newline

done > backbone.txt

