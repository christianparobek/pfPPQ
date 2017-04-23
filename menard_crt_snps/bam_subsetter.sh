## Subset the BAM files
## To see bring to local computer and examine for Menard SNPs
## Started 23 April 2017

names=/proj/ideel/julianog/users/ChristianP/cambodiaWGS/pf/names/our_goods_5x@60%.txt

for name in `cat $names`
do

bedtools intersect -a /proj/ideel/julianog/users/ChristianP/cambodiaWGS/pf/aln/$name.realn.bam -b crt.bed > crt_bams/$name.crt.bam

samtools index crt_bams/$name.crt.bam # index it

done
