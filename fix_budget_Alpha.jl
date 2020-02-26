input="hapblock.txt"
input1="TopIndividuals.txt"
using DataFrames,LinearAlgebra,CSV,DelimitedFiles,SparseArrays,GLPK,JuMP,Statistics #,GLPKMathProgInterface,MathProgBase   
function Alpha(hapblock,topindividuals)
    df = readdlm(input,Int64)#read hapblock file
    df1= readtable(input1,separator = ' ',header=true)# read top individuals sequence file
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

# sort 100 top individuals sequence file
top1=sort(df1[1:100,:ID])
num1=fill(1,100)
df3 = DataFrame(top1=top1,num1=num1)
top1 = collect(1:1:N)
num2=fill(0,N)
df4 =DataFrame(top1=top1,num2=num2)
df5=join(df3, df4, kind = :outer, on = [:top1])
df6=DataFrame(colwise(col -> recode(col, missing=>0), df5), names(df5))
df7=sort!(df6)
#solution of AlphaSeqOpt
num=hcat(df7[:num1]...)
num3=convert(Array{Int64}, num)
#median of haplotypes
colsum=Matrix(mat * num3')
com=hcat(colsum,rmat[:,end])
percenty=Matrix(com)
d=percenty[percenty[:,1].!=0,:]
e=convert(Array{Int64}, d)
j=size(e) 
median(e[:,2]) 
     return j,median(e[:,2])
end

total,med = Alpha(input,input1)

#total number of identified haplotypes selected by 100 animals
println("total number:", total)
#median of frequent haplotypes in the solution
println("median of haplotypes:", med)


