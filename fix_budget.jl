
input="HaplotypesIndividualsCarryPerCore.txt"
using DataFrames,LinearAlgebra,CSV,DelimitedFiles,SparseArrays,GLPK,JuMP #,GLPKMathProgInterface,MathProgBase

function fix_budget(hapblock)
    
    df = readdlm(input,Int64)#read hapblock file
    m,n = size(df)
    
# Column 1 Individual Name
# Column 2-3 Maternal and Paternal Haplotypes block 1
# Column 4-5 Maternal and Paternal Haplotypes block 2
# etc.
# However, 7 in column 3 is not the same haplotype as 7 in column 5
# Each block uses numerical labels to denote haplotypes
    num_cores=Int((n-1)/2)
    
    num_haps=zeros(Int64,num_cores)
    breaks=zeros(Int64,num_cores)
    
# number of haplotypes for each block
    for i in 1:num_cores
        num_haps[i]=length(union(df[:,2*i],df[:,2*i+1]))
    end

    for i in 1:1
        breaks[i]=0
    end

    for i in 2:num_cores
        breaks[i]=num_haps[i-1]+breaks[i-1]
    end

# Need to rewrite the data in a form of a sparse matrix XMat such that 
# Each row is one haplotype
# Each column is one individual
# sparse[i,j,1] if individual j has contains haplotypes i
# Column Sums of XMat are the number of distinct haplotypes that each
# animal possesses. Row Sums of XMat are the number of animals that each
# contain a given haplotype

    rowindex=[]
    colindex=[]
    value=[]
    for k in 1:2
        push!(rowindex,Int64[])
        push!(colindex,Int64[])
        push!(value,Int64[])
    end

    for j in 1:num_cores
        haps=union(df[:,2*j],df[:,2*j+1])
            for i in 1:m
                val1=df[i,2*j]
                val2=df[i,2*j+1]
                ind1=findall(x -> x==val1,haps)
                ind2=findall(x -> x==val2,haps)
                push!(rowindex[1],breaks[j]+ ind1[1])
                push!(rowindex[2],breaks[j]+ ind2[1])
                push!(colindex[1],i)
                push!(colindex[2],i)
                push!(value[1],1)
                push!(value[2],1)
            end
    end
    
    rowindex=[rowindex[1];rowindex[2]]
    colindex=[colindex[1];colindex[2]]
    value=[value[1];value[2]]
    XMat=sparse(rowindex,colindex,value)

#To make sparse matrix only has value "1",one unique haplotypes in one animal in one block was calculated twice. 
    mat=replace!(XMat, 2=>1)

#End of rewrite of data. mat is ready
    M,N =size(mat)
    
#get rowsums of each haplotypes
    rowsums=sum(mat,dims=2)
    rmat=hcat(mat,rowsums)

#remove 1% of frequnecy halotypes from matrix
    a= Int64(N * 0.01)
    freq=rmat[rmat[:,end].>a,:]

#get new matrix
    freq1=freq[:,1:(end-1)]

#End of rewrite of new matrix.
    m,n =size(freq1)
    
# get colsums of each animals which could be coefficients of animals in objective function
    colsums=sum(freq1,dims=1)

# bj = 1, assume all animals will be sequenced.
    c=transpose(ones(Int64,N))

# objective function is integer linear programming
# It is the sum of all the x variables (one x variable for each animal)
# (c1*x1 + c2*x2 + c3*x3 + etc. ) c1, c2 and c3 are colsums. 
# x =1 means animal choosen
# x = 0 means animal not choosen
#
# Part to implement the linear programming in GLPK function
# "max" means maxmize the objective function
# but subject to two constraint 
# each haplotypes will not be tagged at same animals, which set the rhs of first constraint equals to 4.
# given fix budget is 10, which set the rhs of second constraint equals to 10.
# Bin means all x variables constrained to be one or zero

   model = Model(with_optimizer(GLPK.Optimizer))
   @variable(model, x[1:n],Bin)
   @objective(model, Max, sum([colsums[i]*x[i] for i= 1:n]))
   @constraint(model,freq1 * x .<=35)
   @constraint(model,c * x .<= 100)

#Models are solved with the JuMP.optimize! function
   @time optimize!(model)

#check if the solver found an optimal solution
   termination_status(model)

# the primal only inform that the primal solutions are feasible
   primal_status(model)

#final objective solutions
   return objective_value(model)
end

   values = fix_budget(input)

#minimum of animals in the solution
println("Objective of fix_budget: ", values)


