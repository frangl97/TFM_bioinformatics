#!/bin/bash

#### First establish the folder and the location of the PLINK program

FOLDER="/home/frangl97/Escritorio/Filtering_BCF_merge_file"
PLINK="/home/frangl97/Programas/PLINK"


#### Filter by DP 

# The minimum and maximum values for filpering by depth read was 495.85 and 1310.58 respectively.
# Therefore, the minimum value was set to 496 and the maximum value was set to 1311 respectively

bcftools filter -i 'DP>=496 && DP<= 1311' ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed.bcf -o ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_dp.bcf

#Obtain the statistics of the number of variants filtered on this step

bcftools stats ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_dp.bcf > ${FOLDER}/Statistics_BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_dp.txt


#### Filter by QUAL

# Retain all the variants that their quality is greater or equal to 60

bcftools filter -i 'QUAL>=60' ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_dp.bcf -o ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_dp_qual_60.bcf

# Obtain the statistics of the number of variants filtered on this step

bcftools stats ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_dp_qual_60.bcf > ${FOLDER}/Statistics_BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_dp_qual_60.txt


#### Create an id for each variant

# Create an id for each variant because the PLINK program needs them. The id is established combining the chromosome and position of the variant
# The parameter -O is for indicating the type of output file, and b is for obtain a compressed bcf file as output

bcftools annotate --set-id '%CHROM\_%POS' ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_dp_qual_60.bcf -O b -o BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60.bcf


#### Convert the bcf file into vcf file for PLINK

bcftools view ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60.bcf > ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60.vcf


#### Filter by MAF

# Use the parameter --maf 0.05 for filtering the file discarding all variants with a minor allele frequency less than 0.05

${PLINK}/./plink --vcf ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60.vcf --maf 0.05 --make-bed --out ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf --allow-extra-chr

# Once the filtering by maf is done, convert the result into vcf file for viewing the number of samples which have passed the filter

${PLINK}/./plink --bfile ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf --recode vcf --out ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_def --allow-extra-chr

#Now obtain the statistics of the number of variants filtered on this step

bcftools stats ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_def.vcf > ${FOLDER}/Statistics_BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_def.txt


#### Filter by no missing rate

# Use the parameter --geno 0 for filtering the file discarding all variants that are not identified in one of the samples

${PLINK}/./plink --vcf ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_def.vcf --geno 0 --make-bed --out ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0 --allow-extra-chr

# Once the filtering by no missing rate is done, convert again the result into a vcf for viewing the number of samples which have passed the filter

${PLINK}/./plink --bfile ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0 --recode vcf --out ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_def --allow-extra-chr

# Obtain the statistics of the number of variants filtered on this step

bcftools stats ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_def.vcf > ${FOLDER}/Statistics_BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_def.txt


#### Pruning variants in linkage disequilibrium

# For that is used the parameter --indep-parwise for obtaining a file with the ids of the variants that are (out) or not (in) in linkage disequilibrium

${PLINK}/./plink --bfile ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0  --indep-pairwise 50 5 0.2 --out ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_LD  --allow-extra-chr

# Now use bcftools for keeping only the variants that passed are not in linkage disequilibrium (prune.in)

cp ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_LD.prune.in ${FOLDER}/selected_variants.txt

bcftools filter -i "ID==@/home/frangl97/Escritorio/Filtering_BCF_merge_file/selected_variants.txt" ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_def.vcf -o ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld.vcf

# Now use bcftools in order to know the set of variants which have passed the linkage disequilibrium filter

bcftools stats ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld.vcf > ${FOLDER}/Statistics_BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld.txt

# Once we have the file, use bcftools query for renaming the samples, due to plink change them

bcftools query -l ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld.vcf > ${FOLDER}/samples_names.txt

# We have to change the name, and then, apply the following command: 

bcftools reheader -s ${FOLDER}/samples_names.txt ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld.vcf -o ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld_renamed.vcf

#### Discarding the variants with indeterminations (N's) 

# For identifying those variants, use the command 'grep' with the vcf file. Grep the N's, after eliminate the comment rows '#' and obtain only the third column

grep 'N' BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld_renamed.vcf | grep -v '^#' | cut -f3 > variants_to_discard.txt

#For discarding these identified variants, we use the program VCFtools

vcftools --vcf BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld_renamed.vcf --exclude variants_to_discard.txt --recode --out BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld_renamed_def.vcf
