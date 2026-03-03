#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=fastp                             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-00:30:00                 # Maximum job run time
#SBATCH --output=fastp_Diglossa_c_demo-%j.out
#SBATCH --error=fastp_Diglossa_c_demo-%j.err

module load Mamba/23.11.0-0
eval "$(conda shell.bash hook)"

conda activate fastp

shopt -s nullglob

for r1 in /home/g33g348/evo_genomics/course-materials/data/raw_reads/Diglossa_cyanea/*_R1_001.fastq.gz; do
    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    sample=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')

    if [ ! -f "$r2" ]; then
        echo "Missing R2 for $sample"
        continue
    fi

    echo "Processing $sample"

    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$HOME/evo_genomics/students/osterhoudt/Homework_04/second_try/trimmed_reads/${sample}_R1.trimmed.fastq.gz" \
        -O "$HOME/evo_genomics/students/osterhoudt/Homework_04/second_try/trimmed_reads/${sample}_R2.trimmed.fastq.gz" \
        -h "$HOME/evo_genomics/students/osterhoudt/Homework_04/second_try/trimmed_reads/${sample}.fastp.html" \
        -j "$HOME/evo_genomics/students/osterhoudt/Homework_04/second_try/trimmed_reads/${sample}.fastp.json"
done 
        