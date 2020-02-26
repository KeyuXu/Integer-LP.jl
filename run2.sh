#!/bin/bash -l
#SBATCH --job-name=this_is_keyu
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kyuxu@ucdavis.edu

#1) This bash script was used to get median of frequent haplotypes from LP and AlphaSeOpt.
#2) this is the second script based on the output from first script. 

module load julia/1.1.0
module load R

RELPATH=/home/keyu1996/alpha/AlphaSeqOptLinux_2

#generation for each of five populations
gen=(5 10 15 30 50)
#individuals for each of five populations
ind=(6000 11000 16000 31000 51000)

for ((i = 0; i < 5; i++))
do
  Rscript --vanilla pedigree.R ${ind[i]} Pedigree.txt # change name of pedigree for each of populations
# locate TestDataASO1 folder and generate new AlphaSeqOptSpec.txt
  echo "NumberOfIndividualsToSequence ,${ind[i]}" > alpha.txt 
  cat alpha.txt AlphaSeqOpt.txt > AlphaSeqOptSpec.txt
     for a in 1 2 3 4 5 6 7 8 9 10 #replications for each of populations
    do
# locate Phase folder and rename each haplotypes alleles to run AlphaSeqOptLinux.
  cd $RELPATH/TestDataASO1/Phase
  cat $a,${gen[i]}_haps.txt  > Chromosome1.txt
# rename hapblock file as input file.
  cd $RELPATH/TestDataASO1
  cat $a,${gen[i]}_hapblock.txt > hapblock.txt  
  ./AlphaSeqOptLinux
  # remove the last two lines from results of TopIndividualsToSequence.txt.
  head -n -3 TopIndividualsToSequence.txt > TopIndividuals.txt
  # output results of AlphaSeqOpt.
  julia fix_budget_Alpha.jl > $a,${gen[i]}_Alpha.txt 
   done
done