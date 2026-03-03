#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=samtools_stats                            # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-00:30:00                 # Maximum job run time
#SBATCH --output=samtoolsstats.out
#SBATCH --error=samtoolsstats.err

module load Mamba/23.11.0-0
eval "$(conda shell.bash hook)"

conda activate align

samtools flagstat Diglossa_cyanea_161093.sorted.bam > Diglossa_cyanea_161093_alignment_stats.txt
samtools depth Diglossa_cyanea_161093.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > Diglossa_cyanea_161093_depth.txt