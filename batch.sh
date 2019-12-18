#!/bin/bash -l
#SBATCH --job-name=this_is_keyu
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kyuxu@ucdavis.edu

#install AlphaSeqOptLinux in your home directory yourself.
#copy it to your Farm home directory, and unzip the file to get to the binaries.
#move this shell script, pedigree files, AlphaSeqOpt.txt LP.jl and fix_budget.jl to the TestDataASO1 folder in AlphaSeqOptLinux.
#move Alphasim.R to the phase folder in AlphaSeqOptLinux.

module load R
module load julia/1.1.0

bash run.sh > output.txt
