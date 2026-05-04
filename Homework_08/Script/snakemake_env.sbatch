#!/bin/bash

#SBATCH --account=priority-bioe-591-genomics        
#SBATCH --job-name=snake_install                             
#SBATCH --partition=priority              
#SBATCH --nodes=1                       
#SBATCH --ntasks-per-node=1             
#SBATCH --cpus-per-task=1              
#SBATCH --time=0-10:00:00                 
#SBATCH --output=snake_install-%j.out
#SBATCH --error=snake_install-%j.err


module load Mamba/23.11.0-0;
eval "$(conda shell.bash hook)"
mamba env create -f snakemake.yaml