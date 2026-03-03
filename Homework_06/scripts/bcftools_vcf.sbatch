#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=vcf                             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=8                     # number of nodes to allocate
#SBATCH --ntasks-per-node=2            # number of descrete tasks - keep at one 
#SBATCH --cpus-per-task=2             # number of cores to allocate
#SBATCH --time=0-02:30:00                 # Maximum job run time
#SBATCH --output=bcfmpile_vcf-%j.out
#SBATCH --error=bcfmpile_vcf-%j.err

module load Mamba/23.11.0-0
eval "$(conda shell.bash hook)"

conda activate bcftools

bcftools stats Diglossa_cyanea_samples.vcf | grep ^SN > sample_stats.txt
bcftools view -H Diglossa_cyanea_samples.vcf | wc -l > variant_stats.txt
bcftools view -H -v snps Diglossa_cyanea_samples.vcf | wc -l > snps_stats.txt