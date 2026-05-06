#!/bin/bash

#!/bin/bash

#SBATCH --account=priority-bioe-591-genomics        
#SBATCH --job-name=admixture                            
#SBATCH --partition=priority              
#SBATCH --nodes=1                       
#SBATCH --ntasks-per-node=1             
#SBATCH --cpus-per-task=1              
#SBATCH --time=0-10:00:00                 
#SBATCH --output=admixture-%j.out
#SBATCH --error=admixture-%j.err


module load Mamba/23.11.0-0;
eval "$(conda shell.bash hook)"
conda activate admixture

plink2 \
  --vcf /home/g33g348/evo_genomics/students/osterhoudt/Homework_10/squirrels.vcf\
  --allow-extra-chr \
  --make-bed \
  --out species

awk 'BEGIN{OFS="\t"} {if(!($1 in a)) a[$1]=++i; $1=a[$1]; print}' species.bim > species.int.bim
cp species.bed species.int.bed
cp species.fam species.int.fam