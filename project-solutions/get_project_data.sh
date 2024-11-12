#!/bin/bash

# Set the raw repo link
rawdatalink1=ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR208
rawdatalink2=ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR201

# Get lists of the samples to download
LIST1=$(cat << 'END_HEREDOC'
ERR2088830/DC01_R1.fastq.gz
ERR2088830/DC01_R2.fastq.gz
ERR2088831/DC02_R1.fastq.gz
ERR2088831/DC02_R2.fastq.gz
ERR2088832/DC03_R1.fastq.gz
ERR2088832/DC03_R2.fastq.gz
ERR2088833/DC04_R1.fastq.gz
ERR2088833/DC04_R2.fastq.gz
ERR2088834/DC05_R1.fastq.gz
ERR2088834/DC05_R2.fastq.gz
END_HEREDOC
) 

LIST2=$(cat << 'END_HEREDOC'
ERR2015047/SC01_R1.fastq.gz
ERR2015047/SC01_R2.fastq.gz
ERR2015048/SC02_R1.fastq.gz
ERR2015048/SC02_R2.fastq.gz
ERR2015049/SC03_R1.fastq.gz
ERR2015049/SC03_R2.fastq.gz
ERR2015050/SC04_R1.fastq.gz
ERR2015050/SC04_R2.fastq.gz
ERR2015051/SC05_R1.fastq.gz
ERR2015051/SC05_R2.fastq.gz
ERR2015052/SC06_R1.fastq.gz
ERR2015052/SC06_R2.fastq.gz
ERR2015053/SC07_R1.fastq.gz
ERR2015053/SC07_R2.fastq.gz
ERR2015054/SC08_R1.fastq.gz
ERR2015054/SC08_R2.fastq.gz
ERR2015055/SC09_R1.fastq.gz
ERR2015055/SC09_R2.fastq.gz
ERR2015056/SC10_R1.fastq.gz
ERR2015056/SC10_R2.fastq.gz
END_HEREDOC
) 


## Make directory if not exists
echo "Downloading datasets into a directory `project/data`"

mkdir -p data1
for sample in $LIST1; do
  echo "Downloading $sample"
  wget -P data1 -np ${rawdatalink1}/$sample
  gunzip data1/*.gz
done

mkdir -p data2
for sample in $LIST2; do
  echo "Downloading $sample"
  wget -P data2 -np ${rawdatalink2}/$sample
  gunzip data2/*.gz
done