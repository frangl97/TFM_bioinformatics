#!/bin/bash
#
#SBATCH -p eck-q
#SBATCH --chdir=/home/alumno11/LABORATORIO/Mapeo_Genomas_Davis/
#SBATCH -J Mapping_Scaffoldss
#SBATCH --cpus-per-task=10  # Number of CPUs (max is 64)
#SBATCH --mail-type=NONE

FOLDER="/home/alumno11/LABORATORIO/Mapeo_Genomas_Davis/Genome_reads"
GENOME="/home/alumno11/LABORATORIO/Mapeo_Genomas_Davis/Reference_genome/GCA_021292205.2_OSU_Pdul_2.5_genomic.fna"
GENOME_FOLDER="/home/alumno11/LABORATORIO/Mapeo_Genomas_Davis/Reference_genome"

#Fistly the genome is indexed

#bwa index ${GENOME}

for read_folder in ${FOLDER}/* ; do

  # For each genome, we map the reads into the reference using the BWA aligner
  name=$(echo "$read_folder" | cut -d'/' -f 7)
  echo "$name"
  
  # Do the mapping of the reads
  echo "The mapping process for $name begin"
  bwa mem -v 1 -t 10 ${GENOME} ${read_folder}/*trimmed_1.fq.gz ${read_folder}/*trimmed_2.fq.gz | samtools sort -o ${read_folder}/${name}_sorted.bam
  echo "The alignment for $name has ended"
  #Create a summary for the alignment
  samtools flagstat ${read_folder}/${name}_sorted.bam > ${read_folder}/${name}_summary_alignment.txt
  echo "The summary of the alignment is done"
  # Filter the bam file retaining only the reads which mapped only in one position
  samtools view -b -q 40 ${read_folder}/${name}_sorted.bam -o ${read_folder}/${name}_sorted_by_mapq_filtered.bam
  #Count the number of reads in the bam file after the filtering
  samtools view -c ${read_folder}/${name}_sorted_by_mapq_filtered.bam -o ${read_folder}/number_of_counts_af_filtering_by_mapq_${name}.txt
  rm ${read_folder}/${name}_sorted.bam
  echo "The filtering by MAPQ and not primary and supplementaty alignment for $name has ended"
  #Filter the bam file discarding unmapped reads, duplicates reads, secondary and supplementary aligmnents using the flags (-F 3332)
  samtools view -b -F 3332 ${read_folder}/${name}_sorted_by_mapq_filtered.bam -o ${read_folder}/${name}_sorted_filtered.bam
  #Count the number of reads that have passed the filters
  samtools view -c ${read_folder}/${name}_sorted_filtered.bam -o ${read_folder}/number_of_counts_af_filtering_by_dup_unmapped_not_primary_${name}.txt

done

