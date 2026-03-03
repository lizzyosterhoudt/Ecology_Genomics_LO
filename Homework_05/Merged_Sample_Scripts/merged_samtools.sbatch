#!/bin/bash
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=merge.samtools                            # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=4                      # number of nodes to allocate
#SBATCH --ntasks-per-node=12           # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=4              # number of cores to allocate
#SBATCH --time=0-02:30:00                 # Maximum job run time
#SBATCH --output=merge.samtools-%j.out
#SBATCH --error=merge.samtools-%j.err

module load Mamba/23.11.0-0
eval "$(conda shell.bash hook)"
conda activate align

samtools view -b Diglossa_cyanea_merged.sam | samtools sort -o Diglossa_cyanea.merged.sorted.bam
samtools index Diglossa_cyanea.merged.sorted.bam

samtools flagstat Diglossa_cyanea.merged.sorted.bam > Diglossa_cyanea.merged.alignement_stats.txt
samtools depth Diglossa_cyanea.merged.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > Diglossa_cyanea.merged.sorted.depth.txt

