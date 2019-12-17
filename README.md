# Optimzing-sequence-resources-by-ILP

## Introduction
An integer linear programming problem is a mathematical optimization or feasibility program in which some or all of the variables are restricted to be integers. In many settings the term refers to integer linear programming (ILP), in which the objective function and the constraints (other than the integer constraints) are linear. Then we will apply this integer linear programming into optimizing sequencing animals.

In my project, the main objection is to select minimum animals covering all haplotypes in the population and to maximize untility for population given by fix budget of animals. And ILP is compared to an existing method AlphaSeqOpt.

## Requirements
Project is created with:
* R version: 3.6.0
* Julia version: 1.1.1
* AlphaSeqOpt: Linux, Mac OS, Windows

## Setup:
To run this project, there are four steps:
* simulate genotypes or haplotypes data by using Alphasim package in R, `install.packages("ALphaSim")`
* run AlphaSeqOpt to select minimum of animals based on input haplotypes data.
* select minimum of animals in Julia LP based on getting haplotypes blocks from AlphaSeqOpt or coding, `using Pkg; Pkg.add("JuMP"); Pkg.add("GLPK") `
* maximum largest number of haoplotypes given by certain animals by modifying Julia LP. 

## Code Example:
```
$ cd ../AlphaSeqOptLinux_2/TestDataASO1/Phase
$ Rscript --vanilla Alphasim.R
$ cd ../AlphaSeqOptLinux_2/TestDataASO1
$ julia LP.jl 
$ julia fix_budget.jl
```

## Help
For details, please check origin code and shell script.
