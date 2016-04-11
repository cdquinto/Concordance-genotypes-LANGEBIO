###########CONCORDANCE BETWEEN THE MEGA TRAINING RUN WITH PAGE CLUSTER FILE AND 1KG GENOTYPE DATA (AFFY 6.0)
#####################################################################################

###In analysis_entrenamiento_clusterPAGE

###Keep only autosomes and SNPs on the autosomes from the training sample
plink2 --file entrenamiento_clusterPAGE.95.snps.rs --geno 0 --not-chr 0,X,Y,XY,MT --recode --out entrenamiento_clusterPAGE.95.snps.rs_autosomes 

##Find duplicates in the IDs of the SNPs of the training sample
cut -f2 entrenamiento_clusterPAGE.95.snps.rs_autosomes.map > id_snps_entrenamiento_clusterPAGE.95.snps.rs_autosomes.txt
sort id_snps_entrenamiento_clusterPAGE.95.snps.rs_autosomes.txt | uniq -d > id_snps_entrenamiento_clusterPAGE.95.snps.rs_autosomes_duplicates.txt

###Remove the duplicates
plink2 --file entrenamiento_clusterPAGE.95.snps.rs_autosomes --exclude id_snps_entrenamiento_clusterPAGE.95.snps.rs_autosomes_duplicates.txt --recode --out entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates

###Find more duplicates with plink 
plink2 --file entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates --list-duplicate-vars

sed 1d plink.dupvar | cut -f4 > duplicate_entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates.txt
cut -f1 -d" " duplicate_entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates.txt > out1
cut -f2 -d" " duplicate_entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates.txt > out2
cat out1 out2 > duplicate_entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates_2.txt
rm out1 out2
sort duplicate_entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates_2.txt | uniq -d > out

###Remove those SNPs from the MEGA genotypes
plink2 --file entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates --exclude duplicate_entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates_2.txt --recode --out entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2
plink2 --file entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2 --list-duplicate-vars 

###Check if there are any sites with name "."
cut -f2 entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2.map > id_snps_entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates3.txt 
grep -c "\." id_snps_entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates3.txt  
###There are some snps with that name, but should not come up in the next step 

###Find common snps
awk 'NR==FNR{a[$2]=$4;next}{if ($2 in a) print $1,a[$2],$4,$2}' OFS="\t" ../1KG_Affy6.0_training_samples.map entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2.map > overlap_snps.txt
cut -f4 overlap_snps.txt > id_overlap_snps.txt

###
plink2 --file ../1KG_Affy6.0_training_samples --extract id_overlap_snps.txt --recode --make-bed --out 1KG_Affy6.0_training_samples_overlap
plink2 --file entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2 --extract id_overlap_snps.txt --recode --make-bed --out entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap

##Merge
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.bed entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.bim entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.fam --recode --out merged_data
plink2 --bfile entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap --flip merged_data.missnp --make-bed --out source2

##Merge again
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data2

##Remove pb SNPs
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --exclude-snps rs3170863,rs3948464,rs5525 --make-bed --out temp
mv temp.bim 1KG_Affy6.0_training_samples_overlap.bim
mv temp.bed 1KG_Affy6.0_training_samples_overlap.bed
mv temp.fam 1KG_Affy6.0_training_samples_overlap.fam

plink2 --bfile source2 --exclude-snps rs3170863,rs3948464,rs5525 --make-bed --out temp
mv temp.bim entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.bim
mv temp.bed entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.bed
mv temp.fam entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.fam

##Merge again
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.bed entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.bim entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap.fam --recode --out merged_data3

###Remove unnecessary individuals
plink2 --bfile merged_data3 --remove remove_ind.txt --recode --make-bed --out merged_data3_ok

###Final number of SNPs in common between the two datasets: 127,318 SNPs 

####Run cal_concordance.pl (this script makes a barplot)
NA12878	0.919394743869681	0.0806052561303193
HG01938	0.917992742581567	0.0820072574184326
HG01941	0.918707488336292	0.0812925116637082
NA19088	0.918660362242574	0.0813396377574263

######################################################
###########CONCORDANCE BETWEEN THE MEGA TRAINING RUN AND HapMap (plink format)
#####################################################################################

###Find intersect with MEGA run 
awk 'NR==FNR{a[$2]=$4;next}{if ($2 in a) print $1,a[$2],$4,$2}' OFS="\t" ../hapmap3_b37_autosomes_training_samples.map entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2.map > overlap_snps_hapmap_mega.txt
cut -f4 overlap_snps_hapmap_mega.txt > id_overlap_snps_hapmap_mega.txt

plink2 --file ../hapmap3_b37_autosomes_training_samples --extract id_overlap_snps_hapmap_mega.txt --recode --make-bed --out hapmap3_b37_autosomes_training_samples_overlap
plink2 --file entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2 --extract id_overlap_snps_hapmap_mega.txt --recode --make-bed --out entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap

###Merge the datasets and flip
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.bed entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.bim entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.fam --recode --out merged_data_hapmap
plink2 --bfile entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap --flip merged_data_hapmap.missnp --make-bed --out source2

plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data_hapmap_2

###Remove again SNPs that have different genotypes
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --exclude merged_data_hapmap_2.missnp --make-bed --out temp
mv temp.bed hapmap3_b37_autosomes_training_samples_overlap.bed 
mv temp.fam hapmap3_b37_autosomes_training_samples_overlap.fam
mv temp.bim hapmap3_b37_autosomes_training_samples_overlap.bim

plink2 --bfile source2 --exclude merged_data_hapmap_2.missnp --make-bed --out entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap

##Merge again
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.bed entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.bim entrenamiento_clusterPAGE.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.fam --out merged_data_hapmap_3

##Remove inds
plink2 --bfile merged_data_hapmap_3 --remove remove_ind_hapmap.txt --make-bed --recode --out merged_data_hapmap_3_ok

###FINAL FAM FILE:
1 NA-21402 0 0 0 -9
2 NA-12878 0 0 0 -9
3 NA-21405 0 0 0 -9
6 NA-21685 0 0 0 -9
8 NA-19088 0 0 0 -9
1463 NA12878 NA12891 NA12892 2 -9
2571 NA21685 0 0 1 -9
2588 NA21402 0 0 1 -9
2589 NA21405 0 0 1 -9
NA19088 NA19088 0 0 1 -9

###Final number of SNPs in common between the two datasets: 279,879 SNPs

####Results
NA21402	0.940552881781055	0.0594471182189447
NA12878	0.940488568274147	0.0595114317258529
NA21405	0.94052608448651	0.0594739155134898
NA21685	0.940061598047728	0.0599384019522722
NA19088	0.940183079116332	0.0598169208836676



