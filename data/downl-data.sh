#!/bin/bash

# Download raw files
# $ chmod +x downl-data.sh
# $ ./downl-data.sh
# if there is an error "... /bin/bash^M: bad interpreter: No such file or directory"
# https://stackoverflow.com/questions/14219092/bash-script-and-bin-bashm-bad-interpreter-no-such-file-or-directory

# Reference: https://zenodo.org/record/1193466#.X5kngIhKiUk

# Link to where the files are stored 
LIST=$(cat << 'END_HEREDOC'
mut_lib1_R1.fq.gz
mut_lib1_R2.fq.gz
WT_lib1_R1.fq.gz
WT_lib1_R2.fq.gz
Drosophila_melanogaster.BDGP6.dna.fa.gz
Drosophila_melanogaster.BDGP6.85.sample.gtf
END_HEREDOC
) 

rawdatalink=https://zenodo.org/record/1193466/files

## Make directory if not exists
#echo "Creating a directory ./data/"
#mkdir -p ../data/



# Download each fasta read sequence file into the directory
for file in $LIST; do
	    echo "Downloading $file"
	        wget -P ../data -np ${rawdatalink}/$file
done

# If zipped
# echo "Unzipping files"
gunzip ../data/*.fa.gz
