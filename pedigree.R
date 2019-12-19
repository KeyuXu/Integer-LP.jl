#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
ninds = args[1]
args[2] = "Pedigree.txt"

# generate pedigree file
b<-seq(1,ninds)
c<-rep(0,ninds)
d<-rep(0,ninds)
ped<-cbind(b,c,d)
write.table(ped, file = args[2], row.names=FALSE,col.names = FALSE)
