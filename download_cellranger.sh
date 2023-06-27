#!/bin/bash
#SBATCH --mail-user=xxxx 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=11
#SBATCH --time=100:30:00
#SBATCH --mem=20g
#SBATCH --exclusive
#SBATCH -J download.sh

module load cellranger/7.0.0
module load sratoolkit

filename=/home/weir/beegfs/scRNA/'PRJNA823380_sralist.txt'
SraAccList_test=`cat $filename`
for line in $SraAccList_test
do
        echo “Downloading $line”
        
mkdir sample_$line
cd sample_$line
# download 
# --split-files will download full file, -split-3 will download truncated file
fastq-dump --split-files --gzip $line
# rename based on data access in SRA run selector
rename 1.fastq.gz S1_L001_I1_001.fastq.gz *1.fastq.gz
rename 2.fastq.gz S1_L001_I2_001.fastq.gz *2.fastq.gz
rename 3.fastq.gz S1_L001_R1_001.fastq.gz *3.fastq.gz
rename 4.fastq.gz S1_L001_R2_001.fastq.gz *4.fastq.gz
cd ../
cellranger count --id=cellranger_$line \
--fastqs sample_$line \
--sample $line \
--transcriptome /home/weir/beegfs/scRNA/refdata-gex-GRCh38-2020-A
done

# reference:
# publication: https://www.sciencedirect.com/science/article/pii/S0092867422014635?via%3Dihub#abs0010
# github: https://github.com/gatelabnw/csf_aging
# NCBI: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA823380



