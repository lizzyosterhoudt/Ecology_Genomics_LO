#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=vcf_hwk                          # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-00:30:00                 # Maximum job run time
#SBATCH --output=vcf_hwk-%j.out
#SBATCH --error=vcf_hwk-%j.err

module load Mamba/23.11.0-0
eval "$(conda shell.bash hook)"

conda activate vcftools

vcftools\
    --vcf Diglossa_cyanea_samples.vcf --remove-indels --recode --out step1
#step 1 pass


vcftools\
    --vcf step1.recode.vcf --minQ 40 --recode --out step2 

#step 2 pass


vcftools\
    --vcf step2.recode.vcf --max-missing 0.8 --recode --out step3
#step 3 pass 


vcftools\
    --vcf step3.recode.vcf --thin 50 --recode --out step4

#step 4 pass
    
vcftools\
    --vcf step4.recode.vcf --hwe 0.005 --recode --out step5
#step 5 pass

vcftools\
    --vcf step5.recode.vcf --geno-r2 --out step5_ld
    
vcftools\
    --vcf step5.recode.vcf --hardy --out step5_hardy
    
echo -e "file\tsnp_count" > snp_counts.tsv

for f in *.vcf; do
    count=$(bcftools view -H -v snps "$f" | wc -l)
    echo -e "${f}\t${count}" >> snp_counts.tsv
done
    