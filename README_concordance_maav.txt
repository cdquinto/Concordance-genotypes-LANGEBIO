###########CONCORDANCE BETWEEN MAAVP1V1 (maavp1v1.95.snps.rs.map) RUN WITH GLOBAL CLUSTER AND 1KG GENOTYPE DATA (AFFY 6.0)
#####################################################################################

##In /home/cdquinto/test_langebio/analysis_maria_avila, files: maavp1v1.95.snps.rs.ped/map

##Extract PEL sample from 1KG:
HG01935	PEL
##in file: hg01935.txt

###Extract those samples from ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz (only autosomes, and sites with no missing data)
plink2 --vcf /data/Consuelo/1KG_IBS_MXL_Affy6.0_data/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz --keep hg01935.txt --not-chr 0,X,Y,XY,MT --geno 0 --recode --out 1KG_Affy6.0_HG01935
#867, 865 variants

###Extract the same individual from maavp1v1.95.snps.rs.map/ped
48	HG01935_1
96	HG01935_2
##in file: hg01935_maav.txt
plink2 --bfile maavp1v1.95.snps.rs --keep hg01935_maav.txt --not-chr 0,X,Y,XY,MT --geno 0 --recode --make-bed -out maav_hg01935
##1,680,599 SNPs

##Find duplicates in the IDs of the SNPs of the training sample
cut -f2 maav_hg01935.map > id_snps_maav_hg01935.txt
sort id_snps_maav_hg01935.txt | uniq -d > id_snps_maav_hg01935_duplicates.txt

###Remove the duplicates
plink2 --file maav_hg01935 --exclude id_snps_maav_hg01935_duplicates.txt --recode --out maav_hg01935_noduplicates

###Find more duplicates with plink 
plink2 --file maav_hg01935_noduplicates --list-duplicate-vars

sed 1d plink.dupvar | cut -f4 > duplicate_maav_hg01935_noduplicates.txt
cut -f1 -d" " duplicate_maav_hg01935_noduplicates.txt > out1
cut -f2 -d" " duplicate_maav_hg01935_noduplicates.txt > out2
cat out1 out2 > duplicate_maav_hg01935_noduplicates.txt_2.txt
rm out1 out2
sort duplicate_maav_hg01935_noduplicates.txt_2.txt | uniq -d > out

###Remove those SNPs from the MEGA genotypes
plink2 --file maav_hg01935_noduplicates --exclude duplicate_maav_hg01935_noduplicates.txt_2.txt --recode --out maav_hg01935_noduplicates2
plink2 --file maav_hg01935_noduplicates2 --list-duplicate-vars 

###Find common snps
awk 'NR==FNR{a[$2]=$4;next}{if ($2 in a) print $1,a[$2],$4,$2}' OFS="\t" 1KG_Affy6.0_HG01935.map maav_hg01935_noduplicates2.map > overlap_snps.txt
cut -f4 overlap_snps.txt > id_overlap_snps.txt

###
plink2 --file 1KG_Affy6.0_HG01935 --extract id_overlap_snps.txt --recode --make-bed --out 1KG_Affy6.0_HG01935_overlap
plink2 --file maav_hg01935_noduplicates2 --extract id_overlap_snps.txt --recode --make-bed --out maav_hg01935_noduplicates2_overlap

##Merge
plink2 --bfile 1KG_Affy6.0_HG01935_overlap --bmerge maav_hg01935_noduplicates2_overlap.bed maav_hg01935_noduplicates2_overlap.bim maav_hg01935_noduplicates2_overlap.fam --recode --out merged_data
plink2 --bfile maav_hg01935_noduplicates2_overlap --flip merged_data.missnp --make-bed --out source2

##Merge again
plink2 --bfile 1KG_Affy6.0_HG01935_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data2

##Remove pb SNPs
plink2 --bfile 1KG_Affy6.0_HG01935_overlap --exclude-snp rs5525 --make-bed --out temp
mv temp.bim 1KG_Affy6.0_HG01935_overlap.bim
mv temp.bed 1KG_Affy6.0_HG01935_overlap.bed
mv temp.fam 1KG_Affy6.0_HG01935_overlap.fam

plink2 --bfile source2 --exclude-snp rs5525 --make-bed --out temp
mv temp.bim maav_hg01935_noduplicates2_overlap.bim
mv temp.bed maav_hg01935_noduplicates2_overlap.bed
mv temp.fam maav_hg01935_noduplicates2_overlap.fam

###Overlap: 138,145 SNPs

##Merge again
plink2 --bfile 1KG_Affy6.0_HG01935_overlap --bmerge maav_hg01935_noduplicates2_overlap.bed maav_hg01935_noduplicates2_overlap.bim maav_hg01935_noduplicates2_overlap.fam --recode --out merged_data3

####
plink2 --file merged_data3 --distance square allele-ct --out distance_merged_data3

###Results of distance (calculated by allele counts)
48	HG01935_1
96	HG01935_2
HG01935	HG01935

0	0	100100
0	0	100100
100100	100100	0




