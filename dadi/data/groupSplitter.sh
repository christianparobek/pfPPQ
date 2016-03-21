## A VCFtools script to split Pf VCF by CP group
## Started 04 October 2015
## Christian Parobek

## for dadi
sed -i -e '/#CHROM/ s/_20[0-9]*//g' vcf/our.pf.syn.vcf

module add vcftools
vcftools --vcf vcf/our.pf.syn.vcf --keep cp_groups/cp1.txt --recode --out vcf/our.pf.syn.cp1
vcftools --vcf vcf/our.pf.syn.vcf --keep cp_groups/cp2.txt --recode --out vcf/our.pf.syn.cp2
vcftools --vcf vcf/our.pf.syn.vcf --keep cp_groups/cp3.txt --recode --out vcf/our.pf.syn.cp3
vcftools --vcf vcf/our.pf.syn.vcf --keep cp_groups/cp4.txt --recode --out vcf/our.pf.syn.cp4
