## To split VCF and GFF files and put them where they belong
## for reading into PopGenome
## Adapted 19 April 2016

##########
### SET USEFUL VARIABLES
##########

for group in cp3_mine cp3_miotto cp4_mine cp4_miotto
do

	## Make a subsetted VCF
	vcftools --vcf ../joint.vcf --keep groups/$group.txt --recode --out $group

	## Split VCF by Chromosome
	for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
	do

		grep "^#\|^Pf3D7_$chr\_v3" $group.recode.vcf > $group/vcf/chr$chr.vcf

	done

done
