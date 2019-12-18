#!/bin/bash -l

#1) download AlphaSeqOptLinux to your local computer
#from https://alphagenes.roslin.ed.ac.uk/wp/software-2/alphaseqopt/
#2) copy it to your Farm home directory, and unzip the file to get to the binaries.
# scp -P 2022 AlphaSeqOptLinux.zip yourusername@farm.cse.ucdavis.edu:.
# unzip AlphaSeqOptLinux.zip 
#move this shell script, pedigree files, AlphaSeqOpt.txt LP.jl and fix_budget.jl to the TestDataASO1 folder in AlphaSeqOptLinux.
#move Alphasim.R to the phase folder in AlphaSeqOptLinux.

# create absolute path for AlphaSeqOpt
RELPATH=/home/keyu1996/alpha/AlphaSeqOptLinux_2

#generation for each of five populations
gen=(5 10 15 30 50)

#individuals for each of five populations
ind=(6000 11000 16000 31000 51000)

#pedigree for each of five populations
ped=(1 2 3 4 5)

#the range of minimum number of animals in AlphaSeqOptSpec.txt for each populations
range1=(4900 9900 14900 29900 49900)
range2=(6000 11000 16000 31000 51000)

#combine simulation, AlphaSeqOpt and LP
for ((i = 0; i < 5; i++))
do
	cat Pedigree${ped[i]}.txt > Pedigree.txt # change name of pedigree for each of populations
     for a in 1 2 3 4 5  #replications for each of populations
    do

# locate Phase folder and get haplotypes alleles from alphasim.R.
  cd $RELPATH/TestDataASO1/Phase
  Rscript --vanilla Alphasim.R $a ${gen[i]} ${ind[i]} Chromosome1.txt

# save and rename haplotypes alleles files.
  cat Chromosome1.txt > $a,${gen[i]}_haps.txt

# locate TestDataASO1 folder and generate new AlphaSeqOptSpec.txt to run AlphaSeqOptLinux.
  cd $RELPATH/TestDataASO1
  for b in $(seq ${range1[i]} 100 ${range2[i]})
  do
  echo "NumberOfIndividualsToSequence ,$b" > alpha.txt
  cat alpha.txt AlphaSeqOpt.txt > AlphaSeqOptSpec.txt
  ./AlphaSeqOptLinux

# output the results of the last line "OverallPopFootprint" in TopIndividualsToSequence.txt.
# if it's 100.00000 which means we get minimum number of animals. save and rename each results.
  cat TopIndividualsToSequence.txt > $a,$b.txt
done

# output and save each unique haplotypes matrices files.
  cat HaplotypesIndividualsCarryPerCore.txt > $a_${gen[i]}_hapblock.txt

# output results from LP and fixbudget which locate in TestDataASO1 folder.
  julia LP.jl > $a,${gen[i]}_LP.txt
  julia fix_budget.jl > $a,${gen[i]}_fixbudget.txt
    done
done
