#!/bin/bash 


# Reference source:
#
# transcriptome=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz -> unzip
# genome=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.primary_assembly.genome.fa.gz -> 



cd STAR-gencode 

mkdir genomegen

gtf=*.gtf
fasta=*.fa



STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genomegen --genomeFastaFiles $fasta --sjdbGTFfile $gtf 

cd ..
