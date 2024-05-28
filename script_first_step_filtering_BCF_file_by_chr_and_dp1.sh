#!/bin/bash

# Specify the folder where the bcf file is located

FOLDER="/home/frangl97/Escritorio/Filtering_BCF_merge_file"

#First of all, the bcf file has to be indexed. For that use bcftools index
# The parameter -f forces the indexation even if an index already exist 

bcftools index -f {FOLDER}/BCF_merge_Davis_Cebas_Inrae.bcf


#### Filter by chromosomes

# Specifically, almond has 8 chromosomes, and as the reference genome used was Nonpareil, we have to put the name or ids of the 8 chromosomes
# For that, use bcftools view -r, for selecting regions, -O for especify the output format (b is a bcf compressed file) and -o for specify the name of the output file

bcftools view -r CM037988.1,CM037989.1,CM037990.1,CM037991.1,CM037992.1,CM037993.1,CM037994.1,CM037995.1 ${FOLDER}/BCF_merge_Davis_Cebas_Inrae.bcf -O b -o ${FOLDER}/BCF_file_Davis_Cebas_Inrae_filtered_by_chr.bcf


#### Filter by dp > 1

#Eliminate variants that have a dept coverage = 1. For that use bcftools filter. The parameter -i indicates to include all the variants that satisfy the expression

bcftools filter -i 'DP>1' ${FOLDER}/BCF_file_Davis_Cebas_Inrae_filtered_by_chr.bcf -o ${FOLDER}/BCF_file_Davis_Cebas_Inrae_filtered_by_chr_dp1.bcf

# Obtain statistics of the filtering. For that use bcftools stats

bcftools stats ${FOLDER}/BCF_file_Davis_Cebas_Inrae_filtered_by_chr_dp1.bcf > ${FOLDER}/Statistics_BCF_file_Davis_Cebas_Inrae_filtered_by_chr_dp1.bcf

#### Obtain a file for analysing the DP. For that use bcftools query 

bcftools query -f "%CHROM\t%POS\t%DP\t%TYPE\n" ${FOLDER}/BCF_file_Davis_Cebas_Inrae_filtered_by_chr_dp1.bcf > ${FOLDER}/Information_DP_BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp1.txt


### Run the Rmd script for calculating the statistics of DP values



#### Rename the samples 

#Rename the samples in order to change the route of the file for the sample name. For that, use bcftools query -l to obtain a file with the name of each sample and the order. For chosing the new name, add a blank space into the line and write the name that do you want to use

bcftools query -l ${FOLDER}/BCF_file_Davis_Cebas_Inrae_filtered_by_chr_dp1.bcf > ${FOLDER}/sample_names_first_filtering.txt

#After, use the command bcftools reheader -s and the sample name file with the bcf file in order to change the name of the samples from the bcf file

bcftools reheader -s ${FOLDER}/sample_names_first_filtering.txt ${FOLDER}/BCF_file_Davis_Cebas_Inrae_filtered_by_chr_dp1.bcf -o ${FOLDER}/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed.bcf
