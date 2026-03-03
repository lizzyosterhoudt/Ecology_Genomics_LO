#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=merge.bam                            # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=4            # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=2              # number of cores to allocate
#SBATCH --time=0-02:30:00                 # Maximum job run time
#SBATCH --output=merge.bam-%j.out
#SBATCH --error=merge.bam-%j.err

module load Mamba/23.11.0-0
eval "$(conda shell.bash hook)"
conda activate align

samtools merge -@ 16 -b sample_bam_list.txt -O BAM -o Diglossa_cyanea_merged.bam

