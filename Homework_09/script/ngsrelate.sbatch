#!/bin/bash

#SBATCH --account=priority-bioe-591-genomics        
#SBATCH --job-name=ngsrelate                            
#SBATCH --partition=priority              
#SBATCH --nodes=1                       
#SBATCH --ntasks-per-node=1             
#SBATCH --cpus-per-task=1              
#SBATCH --time=0-10:00:00                 
#SBATCH --output=ngsrelate-%j.out
#SBATCH --error=ngsrelate-%j.err

module load Mamba/23.11.0-0;
eval "$(conda shell.bash hook)"

conda activate ngsrelate

# gzip the vcf
bcftools view -Oz -o Koala_MaxMissing10.recode.vcf sample.recode.vcf

# index the vcf
tabix -p vcf Koala_MaxMissing10.recode.vcf

# generate a sample ID file 
bcftools query -l Koala_MaxMissing10.recode.vcf > sample_ids.txt

# call ngsRelate
ngsRelate -h Koala_MaxMissing10.recode.vcf -T GT -c 1 -O sample.ngsrelate.out -z sample_ids.txt