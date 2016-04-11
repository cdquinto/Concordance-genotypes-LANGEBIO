###########CONCORDANCE BETWEEN THE MEGA TRAINING RUN WITH PAGE GLOBAL FILE AND 1KG GENOTYPE DATA (AFFY 6.0)
#####################################################################################

##In /home/cdquinto/test_langebio/analysis_entrenamiento_clusterGLOBAL, files: entrenamiento_clusterGLOBAL.95.snps.rs.ped/map

###Keep only autosomes and SNPs on the autosomes from the training sample
plink2 --file entrenamiento_clusterGLOBAL.95.snps.rs --geno 0 --not-chr 0,X,Y,XY,MT --recode --out entrenamiento_clusterGLOBAL.95.snps.rs_autosomes 

##Find duplicates in the IDs of the SNPs of the training sample
cut -f2 entrenamiento_clusterGLOBAL.95.snps.rs_autosomes.map > id_snps_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes.txt
sort id_snps_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes.txt | uniq -d > id_snps_id_snps_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_duplicates.txt

###Remove the duplicates
plink2 --file entrenamiento_clusterGLOBAL.95.snps.rs_autosomes --exclude id_snps_id_snps_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_duplicates.txt --recode --out entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates

###Find more duplicates with plink 
plink2 --file entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates --list-duplicate-vars

sed 1d plink.dupvar | cut -f4 > duplicate_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates.txt
cut -f1 -d" " duplicate_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates.txt > out1
cut -f2 -d" " duplicate_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates.txt > out2
cat out1 out2 > duplicate_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates.txt_2.txt
rm out1 out2
sort duplicate_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates.txt_2.txt | uniq -d > out

###Remove those SNPs from the MEGA genotypes
plink2 --file entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates --exclude duplicate_entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates.txt_2.txt --recode --out entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2
plink2 --file entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2 --list-duplicate-vars 

###Find common snps
awk 'NR==FNR{a[$2]=$4;next}{if ($2 in a) print $1,a[$2],$4,$2}' OFS="\t" ../1KG_Affy6.0_training_samples.map entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2.map > overlap_snps.txt
cut -f4 overlap_snps.txt > id_overlap_snps.txt

###
plink2 --file ../1KG_Affy6.0_training_samples --extract id_overlap_snps.txt --recode --make-bed --out 1KG_Affy6.0_training_samples_overlap
plink2 --file entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2 --extract id_overlap_snps.txt --recode --make-bed --out entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap

##Merge
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.bed entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.bim entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.fam --recode --out merged_data
plink2 --bfile entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap --flip merged_data.missnp --make-bed --out source2

##Merge again
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data2

##Remove pb SNPs
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --exclude-snps rs3170863,rs3948464,rs5525 --make-bed --out temp
mv temp.bim 1KG_Affy6.0_training_samples_overlap.bim
mv temp.bed 1KG_Affy6.0_training_samples_overlap.bed
mv temp.fam 1KG_Affy6.0_training_samples_overlap.fam

plink2 --bfile source2 --exclude-snps rs3170863,rs3948464,rs5525 --make-bed --out temp
mv temp.bim entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.bim
mv temp.bed entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.bed
mv temp.fam entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.fam

##Merge again
plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.bed entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.bim entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap.fam --recode --out merged_data3

###Remove unnecessary individuals
plink2 --bfile merged_data3 --remove ../analysis_entrenamiento_clusterPAGE/remove_ind.txt --recode --make-bed --out merged_data3_ok

###Final number of SNPs in common between the two datasets: 122,581 SNPs 

####Run cal_concordance.pl (this script makes a barplot)
NA12878	0.91833971006926	0.0816602899307397
HG01938	0.916908003687358	0.0830919963126423
HG01941	0.917605501668285	0.0823944983317153
NA19088	0.917556554441553	0.0824434455584471

##Moved files to /home/cdquinto/test_langebio/analysis_entrenamiento_clusterGLOBAL/concordance_1KG_MEGA

#####################################################################################
###########CONCORDANCE BETWEEN THE MEGA TRAINING RUN AND HapMap (plink format)
#####################################################################################

###Find intersect with MEGA run 
awk 'NR==FNR{a[$2]=$4;next}{if ($2 in a) print $1,a[$2],$4,$2}' OFS="\t" ../hapmap3_b37_autosomes_training_samples.map entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2.map > overlap_snps_hapmap_mega.txt
cut -f4 overlap_snps_hapmap_mega.txt > id_overlap_snps_hapmap_mega.txt

plink2 --file ../hapmap3_b37_autosomes_training_samples --extract id_overlap_snps_hapmap_mega.txt --recode --make-bed --out hapmap3_b37_autosomes_training_samples_overlap
plink2 --file entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2 --extract id_overlap_snps_hapmap_mega.txt --recode --make-bed --out entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap

###Merge the datasets and flip
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.bed entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.bim entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.fam --recode --out merged_data_hapmap
plink2 --bfile entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap --flip merged_data_hapmap.missnp --make-bed --out source2

plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data_hapmap_2

###Remove again SNPs that have different genotypes
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --exclude merged_data_hapmap_2.missnp --make-bed --out temp
mv temp.bed hapmap3_b37_autosomes_training_samples_overlap.bed 
mv temp.fam hapmap3_b37_autosomes_training_samples_overlap.fam
mv temp.bim hapmap3_b37_autosomes_training_samples_overlap.bim

plink2 --bfile source2 --exclude merged_data_hapmap_2.missnp --make-bed --out entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap

##Merge again
plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.bed entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.bim entrenamiento_clusterGLOBAL.95.snps.rs_autosomes_noduplicates2_overlap_hapmap.fam --out merged_data_hapmap_3

##Remove inds
plink2 --bfile merged_data_hapmap_3 --remove remove_ind_hapmap.txt --make-bed --recode --out merged_data_hapmap_3_ok

###Final number of SNPs in common between the two datasets: 268,907 SNPs

####Results
NA21402	0.939557170322825	0.0604428296771746
NA12878	0.939516263987178	0.0604837360128223
NA21405	0.939560889080612	0.0604391109193885
NA21685	0.939172278891959	0.060827721108041
NA19088	0.939207607090927	0.0607923929090727

