#!bin/bash

mkdir output_ens 

input_path="../rawdata"
ref_path="/home/mira/Documents/programming/Bioinformatics/STAR-test/STAR-ensembl"
genome_fasta="hg19.fa"
GTF="hg19.ensGene.gtf"

cd rawdata

input_files=$(ls)

cd ../output_ens 



for read in $input_files 
do 
    STAR --runThreadN 8 --runMode alignReads --genomeDir ${ref_path}/genomegen --sjdbGTFfile ${ref_path}/${GTF} -sjdbOverhang 100 --readFilesIn ${input_path}/${read} --outFileNamePrefix ${read} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic --chimOutType Junctions

done

