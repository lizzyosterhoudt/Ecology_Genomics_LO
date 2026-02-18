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
#SBATCH --output=/home/g33g348/evo_genomics/students/osterhoudt/04_homework/err_out_files/219530_fastp_Diglossa_c_demo-%j.out
#SBATCH --error=/home/g33g348/evo_genomics/students/osterhoudt/04_homework/err_out_files/219530_fastp_Diglossa_c_demo-%j.err

# load module and activate mamba
module load Mamba/23.11.0-0
eval "$(conda shell.bash hook)"

# activate your fastp environment
conda activate fastp

# run fastp
fastp \
-i ~/evo_genomics/students/osterhoudt/04_homework/Diglossa_cyanea/Diglossa_cyanea_219530_R1_001.fastq.gz \
-I ~/evo_genomics/students/osterhoudt/04_homework/Diglossa_cyanea/Diglossa_cyanea_219530_R2_001.fastq.gz \
-o ~/evo_genomics/students/osterhoudt/04_homework/trimmed_reads/Diglossa_cyanea_219530_R1_001_trimmed.fastq.gz \
-O ~/evo_genomics/students/osterhoudt/04_homework/trimmed_reads/Diglossa_cyanea_219530_R2_001_trimmed.fastq.gz \
-h ~/evo_genomics/students/osterhoudt/04_homework/trimmed_reads/219530_sample.fastp.html \
-j ~/evo_genomics/students/osterhoudt/04_homework/trimmed_reads/219530_sample.fastp.json

