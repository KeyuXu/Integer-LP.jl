#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
seed = args[1]
ngenerations = args[2]
ninds = args[3]
args[4] = "Chromosome1.txt"
#output_path_hap =paste(seed,"Hap_",ngenerations,".txt")
#output_path_geno = paste(seed,"Geno_",ngenerations,ninds,".txt")

library(AlphaSimR)
set.seed(seed)
# Creating 1000 Founder Haplotypes with 10 chr, chrmosome length equals to 0.015M
founderPop = runMacs(nInd=1000, nChr=10, segSites=1000,manualGenLen=0.015)
# Setting Simulation Parameters
SP = SimParam$new(founderPop)
#the trait was controlled by 15 QTL per chromosome
SP$addTraitA(nQtlPerChr=15)
#equal ratio of sex for each generation
SP$setGender("yes_sys")
#snp array loci
SP$addSnpChip(1000)
#generate the initial population from founderpop
pop = newPop(founderPop,simParam=SP)
popt=c(pop)
#modeling 5 generation of selection and mating
for(generation in 1:ngenerations){
  pop = selectCross(pop=pop, nFemale=500, nMale=25, use="bv", nCrosses=1000,selectTop = TRUE,simParam=SP)
  popt = c(popt,pop)
}
# output haplotypes
a<-pullSnpHaplo(popt, simParam=SP)
# add animals id
ind1<-rep(c(1:ninds), each = 2)
a3<-cbind(ind1,a)
write.table(a3,file = args[4], row.names=FALSE,col.names = FALSE)
# output genotypes
#a1<-pullSnpGeno(popt, simParam=SP)
# add animals id
#ind<-c(1:ninds)
#a2<-cbind(ind,a1)
#write.table(a2,output_path_geno,row.names=FALSE,col.names = FALSE)
