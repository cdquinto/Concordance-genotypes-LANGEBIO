###########CONCORDANCE BETWEEN THE MEGA TRAINING RUN AND 1KG GENOTYPE DATA (AFFY 6.0)
#####################################################################################

###Samples genotyped during the MEGA training (run without cluster file):
(I modify the ID to match the ones in the sample sheet from 1KG Affy 6.0 genotype data):
Sample	Sample_ID_1KG	FamilyID	Population
NA-21402 => not included in 1KG 
NA-12878	NA12878	CEU	CEU
NA-21405 => not included in 1KG 
NA-21732 => not included in 1KG 
HG-01938	HG01938	PEL PEL
NA-21685 => not included in 1KG
HG-01941	HG01941 PEL	PEL
NA-19088	NA19088	JPT	JPT

###The information of these samples are in the file: training_samples.txt
NA12878	NA12878
HG01938	HG01938
HG01941	HG01941
NA19088	NA19088

###Extract those samples from ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz (only autosomes, and sites with no missing data)
plink2 --vcf /data/Consuelo/1KG_IBS_MXL_Affy6.0_data/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz --keep training_samples.txt --not-chr 0,X,Y,XY,MT --geno 0 --recode --out 1KG_Affy6.0_training_samples
#863597 variants

###Find intersect of SNPs in the MEGA training genotyping to rename the SNPs
##I used the files (from Pavel): entrenamiento_16sample (map/ped)
awk 'NR==FNR{if($2 != "." && $2 !~ /,/){a[$1]=$2};next}{if($2 in a)sub($2,a[$2]);print}' OFS="\t" MEGA_Consortium_v2_15070954_A1_b138_rsids.txt entrenamiento_16sample.map > entrenamiento_16sample_renamed.map 
mv entrenamiento_16sample_renamed.map entrenamiento_16sample.map

###Keep only autosomes and SNPs on the autosomes from the training sample
plink2 --file entrenamiento_16sample --geno 0 --not-chr 0,X,Y,XY,MT --recode --out entrenamiento_16sample_nomissing 

##Find duplicates in the IDs of the SNPs of the training sample
cut -f2 entrenamiento_16sample_nomissing.map > id_snps_entrenamiento_16sample_nomissing.txt
sort id_snps_entrenamiento_16sample_nomissing.txt | uniq -d > id_snps_entrenamiento_16sample_nomissing_duplicates.txt

###Remove the duplicates
plink2 --file entrenamiento_16sample_nomissing --exclude id_snps_entrenamiento_16sample_nomissing_duplicates.txt --recode --out entrenamiento_16sample_nomissing_noduplicates

###Find more duplicates with plink 
plink2 --file entrenamiento_16sample_nomissing_noduplicates --list-duplicate-vars

###Get the id of the duplicate SNPs
sed 1d plink.dupvar | cut -f4 > duplicate_entrenamiento_16sample_nomissing_noduplicates.txt
cut -f1 -d" " duplicate_entrenamiento_16sample_nomissing_noduplicates.txt > out1
cut -f2 -d" " duplicate_entrenamiento_16sample_nomissing_noduplicates.txt > out2
cat out1 out2 > duplicate_entrenamiento_16sample_nomissing_noduplicates_2.txt
rm out1 out2
sort duplicate_entrenamiento_16sample_nomissing_noduplicates_2.txt | uniq -d > out

###Remove those SNPs from the MEGA genotypes
plink2 --file entrenamiento_16sample_nomissing_noduplicates --exclude duplicate_entrenamiento_16sample_nomissing_noduplicates_2.txt --recode --out entrenamiento_16sample_nomissing_noduplicates2
plink2 --file entrenamiento_16sample_nomissing_noduplicates2 --list-duplicate-vars 

###Check if there are any sites with name "."
cut -f2 entrenamiento_16sample_nomissing_noduplicates2.map > id_snps_entrenamiento_16sample_nomissing_noduplicates2.txt 
cut -f2 1KG_Affy6.0_training_samples.map > id_1KG_Affy6.0_training_samples.txt 

grep -c "\." id_snps_entrenamiento_16sample_nomissing_noduplicates2.txt 
grep -c "\." id_1KG_Affy6.0_training_samples.txt
###There are some snps with that name, but should not come up in the next step 

###Find common snps
awk 'NR==FNR{a[$2]=$4;next}{if ($2 in a) print $1,a[$2],$4,$2}' OFS="\t" 1KG_Affy6.0_training_samples.map entrenamiento_16sample_nomissing_noduplicates2.map > overlap_snps.txt
sort -n -k1,1n overlap_snps.txt > overlap_snps_sorted.txt

cut -f4 overlap_snps_sorted.txt > id_overlap_snps_sorted.txt

plink2 --file 1KG_Affy6.0_training_samples --extract id_overlap_snps_sorted.txt --recode --out 1KG_Affy6.0_training_samples_overlap
plink2 --file entrenamiento_16sample_nomissing_noduplicates2 --extract id_overlap_snps_sorted.txt --recode --out entrenamiento_16sample_nomissing_noduplicates2_overlap

plink2 --file 1KG_Affy6.0_training_samples_overlap --recode --make-bed --out 1KG_Affy6.0_training_samples_overlap
plink2 --file entrenamiento_16sample_nomissing_noduplicates2_overlap --recode --make-bed --out entrenamiento_16sample_nomissing_noduplicates2_overlap

###Merge the datasets (in binary format)
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge entrenamiento_16sample_nomissing_noduplicates2_overlap.bed entrenamiento_16sample_nomissing_noduplicates2_overlap.bim entrenamiento_16sample_nomissing_noduplicates2_overlap.fam --recode --out merged_data
###Some of the SNPs are in the complementary, so I have to change it with plink, I have to flip the strands in the entrenamiento dataset
plink2 --bfile entrenamiento_16sample_nomissing_noduplicates2_overlap --flip merged_data.missnp --make-bed --out source2

###Merge again 
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data2

##Remove 3 SNPs that appear to be triallelic (in merged_data2.missnp)
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --exclude-snps rs3170863,rs3948464,rs5525 --make-bed --out temp
mv temp.bim 1KG_Affy6.0_training_samples_overlap.bim
mv temp.bed 1KG_Affy6.0_training_samples_overlap.bed
mv temp.fam 1KG_Affy6.0_training_samples_overlap.fam

plink2 --bfile entrenamiento_16sample_nomissing_noduplicates2_overlap --exclude-snps rs3170863,rs3948464,rs5525 --make-bed --out temp
mv temp.bim entrenamiento_16sample_nomissing_noduplicates2_overlap.bim
mv temp.bed entrenamiento_16sample_nomissing_noduplicates2_overlap.bed
mv temp.fam entrenamiento_16sample_nomissing_noduplicates2_overlap.fam

plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge entrenamiento_16sample_nomissing_noduplicates2_overlap.bed entrenamiento_16sample_nomissing_noduplicates2_overlap.bim entrenamiento_16sample_nomissing_noduplicates2_overlap.fam --recode --out merged_data
###Some of the SNPs are in the complementary, so I have to change it with plink, I have to flip the strands in the entrenamiento dataset
plink2 --bfile entrenamiento_16sample_nomissing_noduplicates2_overlap --flip merged_data.missnp --make-bed --out source2

###Merge again 
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data2

##Remove SNPs with weird positions
plink2 --bfile merged_data2 --exclude-snps rs9480186,rs9522257 --make-bed --out merged_data3
plink2 --bfile merged_data3 --recode --out merged_data3

###Remove unnecessary individuals
plink2 --bfile merged_data3 --remove remove_ind.txt --recode --out merged_data3_ok
plink2 --file merged_data3_ok --make-bed --out merged_data3_ok

###FINAL FAM FILE: 

10 NA-12878 0 0 0 -9
13 HG-01938 0 0 0 -9
15 HG-01941 0 0 0 -9
16 NA-19088 0 0 0 -9
HG01938 HG01938 0 0 0 -9
HG01941 HG01941 0 0 0 -9
NA12878 NA12878 0 0 0 -9
NA19088 NA19088 0 0 0 -9

###Final number of SNPs in common between the two datasets: 128, 587 SNPs

####Run cal_concordance.pl (this script makes a barplot)
NA12878	242311	257174	0.94220644388624	0.0577935561137596
HG01938	241967	257174	0.940868828108596	0.0591311718914043
HG01941	242138	257174	0.94153374757946	0.0584662524205402
NA19088	242078	257174	0.941300442501964	0.0586995574980363

######################################################
###########CONCORDANCE BETWEEN THE MEGA TRAINING RUN AND HapMap (plink format)
#####################################################################################

####in tikal 
####/data/HighDensitySNPData/hapmap3_r2_b36_fwd.consensus.qc.poly_autosomal.fam(bed/bim)

###The information of these samples are in the file: training_samples_HapMap.txt
2588 NA21402
1463 NA12878
2589 NA21405
2571 NA21685
NA19088	NA19088

###Extract those individuals from plink files 
plink2 --bfile /data/HighDensitySNPData/HapMap_Plink_hapmap3_r2_b36/hapmap3_r2_b36_fwd.consensus.qc.poly_autosomal --keep training_samples_HapMap.txt --geno 0 --recode --out hapmap3_b36_autosomes_training_samples 

###Liftover of positions from b36 to b37
perl test.pl 
liftOver hapmap3_b36_autosomes_to_liftOver.txt /data/LiftoverFiles/hg18ToHg19.over.chain hapmap3_b37_autosomes_to_liftOver.txt unMapped

###Get from unMapped the SNPs whose positions were deleted in the new build
awk 'NR%2==0' unMapped | cut -f4 > snps_to_remove_hapmap.txt

###Remove those SNPs from the plink files
plink2 --file hapmap3_b36_autosomes_training_samples --exclude snps_to_remove_hapmap.txt --recode --make-bed --out hapmap3_b37_autosomes_training_samples

###Update map 
cut -f4 hapmap3_b37_autosomes_to_liftOver.txt  > temp1
cut -f3 hapmap3_b37_autosomes_to_liftOver.txt > temp2
rm temp1 temp2
paste temp1 temp2 > new_pos_snps_hapmap3_b37.txt

plink2 --file hapmap3_b37_autosomes_training_samples --update-map new_pos_snps_hapmap3_b37.txt --make-bed --recode --out temp
mv temp.bed hapmap3_b37_autosomes_training_samples.bed 
mv temp.bim hapmap3_b37_autosomes_training_samples.bim 
mv temp.fam hapmap3_b37_autosomes_training_samples.fam 
mv temp.ped hapmap3_b37_autosomes_training_samples.ped
mv temp.map hapmap3_b37_autosomes_training_samples.map

###Find intersect with MEGA run 
awk 'NR==FNR{a[$2]=$4;next}{if ($2 in a) print $1,a[$2],$4,$2}' OFS="\t" hapmap3_b37_autosomes_training_samples.map entrenamiento_16sample_nomissing_noduplicates2.map > overlap_snps_hapmap_mega.txt
cut -f4 overlap_snps_hapmap_mega.txt > id_overlap_snps_hapmap_mega.txt

###
plink2 --file hapmap3_b37_autosomes_training_samples --extract id_overlap_snps_hapmap_mega.txt --recode --make-bed --out hapmap3_b37_autosomes_training_samples_overlap
plink2 --file entrenamiento_16sample_nomissing_noduplicates2 --extract id_overlap_snps_hapmap_mega.txt --recode --make-bed --out entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap

###Merge the datasets
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap.bed entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap.bim entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap.fam --recode --out merged_data_hapmap
plink2 --bfile entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap --flip merged_data_hapmap.missnp --make-bed --out source2

plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data_hapmap_2

###Remove again SNPs that have different genotypes
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --exclude merged_data_hapmap_2.missnp --make-bed --out temp
mv temp.bed hapmap3_b37_autosomes_training_samples_overlap.bed 
mv temp.fam hapmap3_b37_autosomes_training_samples_overlap.fam
mv temp.bim hapmap3_b37_autosomes_training_samples_overlap.bim

plink2 --bfile source2 --exclude merged_data_hapmap_2.missnp --make-bed --out entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap

##Merge again
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap.bed entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap.bim entrenamiento_16sample_nomissing_noduplicates2_overlap_hapmap.fam --out merged_data_hapmap_3

###Keep individuals
plink2 --bfile merged_data_hapmap_3 --keep keep_ind_merged_data_hapmap_3.txt --make-bed --recode --out merged_data_hapmap_3_ok

##FINAL FAM FILE
9 NA-21402 0 0 0 -9
10 NA-12878 0 0 0 -9
11 NA-21405 0 0 0 -9
14 NA-21685 0 0 0 -9
16 NA-19088 0 0 0 -9
1463 NA12878 NA12891 NA12892 2 -9
2571 NA21685 0 0 1 -9
2588 NA21402 0 0 1 -9
2589 NA21405 0 0 1 -9
NA19088 NA19088 0 0 1 -9

###Final number of SNPs in common between the two datasets: 282,831 SNPs

####Run cal_concordance.pl
2588	0.961344407084089	0.0386555929159109
1463	0.961420424210925	0.0385795757890754
2589	0.961347942764407	0.038652057235593
2571	0.960757484151313	0.039242515848687
NA19088	0.960985535531819	0.0390144644681807

#####I moved the results to this directory: analysis_entrenamiento_nocluster 

