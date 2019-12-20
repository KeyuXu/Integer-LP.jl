#!/bin/bash -l
#SBATCH --job-name=this_is_keyu
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kyuxu@ucdavis.edu

#1)download AlphaSeqOptLinux to your local computer
#from https://alphagenes.roslin.ed.ac.uk/wp/software-2/alphaseqopt/
#2)copy it to your Farm home directory, and unzip the file to get to the binaries.
#scp -P 2022 AlphaSeqOptLinux.zip yourusername@farm.cse.ucdavis.edu:.
#unzip AlphaSeqOptLinux.zip.
#3)download Integer-LP.jl-master zip from Github: https://github.com/KeyuXu/Integer-LP.jl
#4)copy it to your Farm home directory, and unzip it.
#unzip Integer-LP.jl-master.zip.
#5)move run.sh, AlphaSeqOpt.txt, LP.jl and fix_budget.jl these four files to the TestDataASO1 folder 
#in AlphaSeqOptLinux_2 which released from AlphaSeqOptLinux.zip.
#cd AlphaSeqOptLinux_2
#touch run.sh 
#mv run.sh ~/AlphaSeqOptLinux_2/TestDataASO1
#move other files in the same way.
#6)move Alphasim.R to the phase folder in AlphaSeqOptLinux_2.
#touch Alphasim.R 
#mv Alphasim.R ~/AlphaSeqOptLinux_2/TestDataASO1/phase
#7)Then we can submit this job in the server.
#cd ~RELPATH/TestDataASO1
#sbatch run.sh

module load R
module load julia/1.1.0

#change obsolute path based on yourself.
RELPATH=/home/keyu1996/alpha/AlphaSeqOptLinux_2

#generation for each of five populations
gen=(5 10 15 30 50)
#individuals for each of five populations
ind=(6000 11000 16000 31000 51000)
#pedigree for each of five populations
ped=(1 2 3 4 5)
#the range of minimum number of animals in AlphaSeqOptSpec.txt for each populations
range1=(4900 8400 11350 17200 20300)
range2=(5200 8600 11600 17500 21000)

for ((i = 0; i < 5; i++))
do 
	Rscript --vanilla pedigree.R ${ind[i]} Pedigree.txt # get pedigrees for each populations 
     for a in 1 2 3  #replications for each of populations
    do
# locate Phase folder and get haplotypes alleles from alphasim.R.
  cd $RELPATH/TestDataASO1/Phase 
  Rscript --vanilla Alphasim.R $a ${gen[i]} ${ind[i]} Chromosome1.txt 
# save and rename haplotypes alleles files.  
  cat Chromosome1.txt > $a,${gen[i]}_haps.txt  
# locate TestDataASO1 folder and generate new AlphaSeqOptSpec.txt to run AlphaSeqOptLinux.
  cd $RELPATH/TestDataASO1 
  for b in $(seq ${range1[i]} 50 ${range2[i]}) 
  do
  echo "NumberOfIndividualsToSequence ,$b" > alpha.txt 
  cat alpha.txt AlphaSeqOpt.txt > AlphaSeqOptSpec.txt
  ./AlphaSeqOptLinux
# output the results of the last line "OverallPopFootprint" in TopIndividualsToSequence.txt.
# if it's 100.00000 which means we get minimum number of animals. save and rename each results.
  cat TopIndividualsToSequence.txt > $a,$b.txt   
done
# output and save each unique haplotypes matrices files.
  cat HaplotypesIndividualsCarryPerCore.txt > $a,${gen[i]}_hapblock.txt 
# output results from LP which locate in TestDataASO1 folder.
  julia LP.jl > $a,${gen[i]}_LP.txt 
    done
done
