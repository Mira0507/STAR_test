#!/bin/bash 

# Source of reference:
# FASTA: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ ---> hg19.fa.gz
# GTF: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/  ---> hg19.ensGene.gtf.gz 

# For more info, check this out: 
# https://genome.ucsc.edu/FAQ/FAQgenes.html#hg19


fasta="hg19.fa"

gtf="hg19.ensGene.gtf"

mkdir genomegen

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./genomegen --genomeFastaFiles $fasta --sjdbGTFfile $gtf 



