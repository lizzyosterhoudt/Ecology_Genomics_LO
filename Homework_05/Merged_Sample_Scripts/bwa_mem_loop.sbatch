#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=align_loop                            # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=4            # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=2              # number of cores to allocate
#SBATCH --time=0-02:30:00                 # Maximum job run time
#SBATCH --output=align_loop-%j.out
#SBATCH --error=align_loop-%j.err

module load Mamba/23.11.0-0
eval "$(conda shell.bash hook)"
conda activate align


READDIR=/home/g33g348/evo_genomics/students/osterhoudt/Homework_05/trimmed_reads
OUTDIR=/home/g33g348/evo_genomics/students/osterhoudt/Homework_05
REF=/home/g33g348/evo_genomics/students/osterhoudt/Homework_05/hemoglobin_references.fasta
shopt -s nullglob

for r1 in "$READDIR"/*_R1_001_trimmed.fastq.gz; do

    r2="${r1/_R1_001_trimmed.fastq.gz/_R2_001_trimmed.fastq.gz}"
    sample=$(basename "$r1" | sed 's/_R1_001_trimmed.fastq.gz//')

    if [ ! -f "$r2" ]; then
        echo "Missing R2 for $sample"
        continue
    fi

    echo "Processing $sample"

    bwa mem -t 4 "$REF" "$r1" "$r2" | \
        samtools sort -o "$OUTDIR/${sample}.sorted.bam"

done