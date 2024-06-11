#!/bin/bash

# Command for creating the Nonpareil database after following the instructions of the Manual

java -jar snpEff.jar build -gtf22 -v prunus.dulcis.nonpareil -noCheckCds -noCheckProtein

# Command for performing the snpEff analysis 

java -Xmx8g -jar snpEff.jar data/prunus.dulcis.nonpareil/ /home/frangl97/Escritorio/Resultados_snpEff_Davis_Cebas_Inrae_geno_0/BCF_file_merge_Davis_Cebas_Inrae_filtered_by_chr_dp_1_renamed_set_ids_dp_qual_60_maf_geno_0_ld_renamed_def.vcf > /home/frangl97/Escritorio/Resultados_snpEff_Davis_Cebas_Inrae_geno_0/output_snpEff.vcf
