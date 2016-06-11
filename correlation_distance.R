#!/usr/bin/env Rscript

########################################################################
##
## File: correlation_distance.R
##
## Created by: Anja S.
## Created on: June 06, 2016
##
## Description:   Calculates correlation or distance between the 
##  (haar wavelet transformed) gene expression profiles for each gene pair.
##     
## Input:   File with (haar wavelet transformed) gene expression profiles.
##
## Output: File with GeneIDs and the their distance or 
##  correlation matrix (saved as lower triangle).
##          
########################################################################

#install CRAN packages
#install.packages("optparse")
#install.packages("infotheo")
#install.packages("schoolmath")

#load packages
library(optparse) #to parse command line options
library(infotheo) #to calculate the mutual information
library(schoolmath)

#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help=".dat file with (e.g. haar wavelet transformed) gene expression data from the Allen Brain Atlas"),
  make_option(c("-o","--output"), action="store", help="result folder (voluntary)"),
  
  #choose at least one measurement to compare gene expression profiles
  make_option(c("-e","--euclidean"), action="store_true", default=FALSE, help="Compare gene expression profiles with euclidean distance"),
  make_option(c("-s","--spearman"), action="store_true", default=FALSE, help="Compare gene expression profiles with spearman correlation coefficient"),
  make_option(c("-p","--pearson"), action="store_true", default=FALSE, help="Compare gene expression profiles with pearson correlation coefficient"),
  make_option(c("-r","--rio"), action="store_true", default=FALSE, help="Compare gene expression profiles with RIO (Relative intensity overlap)"),
  make_option(c("-m","--mutual_information"), action="store_true", default=FALSE, help="Compare gene expression profiles with mutual information")
  
)

#parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

#check if file with the gene expression data is given
if (is.null(opt$input)) {
  stop("No file with (haar wavelet transformed) gene expression data is given!")
}

#check if at least one measurement (to compare gene expression profiles) has been selected
if(!any(opt$euclidean, opt$spearman, opt$pearson, opt$rio, opt$mutual_information)){
  stop("Select a measurement to compare gene expression profiles!")
}

#user specified an output/result folder
if (!is.null(opt$output)) {
  output_folder<-opt$output
  
}else{ #standard output folder
  output_folder<-"output"
}

#create output folder
dir.create(path=output_folder, recursive=TRUE)


########################################Functions###############################################

##The functions comp_euclidean, comp_spearman, comp_pearson, comp_rio and comp_mutual_information 
##calculate the lower triangle of the correlation or distance matrix for the given genes.
##Due to the high memory consumption, only the gene expression profile of gene i and a limited 
##number of gene expression profiles of gene j are loaded per iteration.

##gene j    a b c ...               
##gene i  
##    a    1
##    b    2 3 
##    c    4 5 6                  
##   ...                 


#computes euclidean distance between gene expression profiles for each gene pair (gene_i, gene_j)
comp_euclidean<-function(inputfile, number_genes, gene_limit, outputfile){
  
  #first connection -> to read gene expression profile of the first gene of each gene pair
  connection1<-file(inputfile)
  open(connection1)
  #skip header + first gene
  gene1<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE,  skip=1, nlines=1)

  #select gene expression profile from first gene
  for (i in 2:number_genes){
    
    gene_i<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE, nlines=1)[-1] 
    
    #saves distance between gene i to gene_1, gene_2, ... gene_i-1
    distance_array<-numeric(length=i-1)
    
    #second connection -> to read gene expression profiles of the second genes of each gene pair
    #(only a limited number can be read because of memory consumption)
    connection2<-file(inputfile)
    open(connection2)
    #skip headline
    genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, skip=1, nlines=gene_limit)
    genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1] #remove column with the GeneIDs
    
    #start with the first gene in genes_j
    pos=1
    
    #select gene expression profile from second gene
    for(j in 1:(i-1)){
      
      #load new gene expression profiles (second genes)
      if(!is.decimal((j-1)/gene_limit) && j!=1){
        
        pos=1
        #a full set of genes can be loaded
        if(i-j>=gene_limit){
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=gene_limit)
          genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1]
        }
        #there are only a few genes to load
        else{
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=i-j)
          genes_j<-matrix(genes_j, nrow=i-j, byrow=TRUE)[,-1]
        }
      }
      
      #genes_j contains more than one gene
      if(is.matrix(genes_j)){
        distance_array[j]<-sqrt(sum((gene_i- genes_j[pos,])^2))
        
      } else{ #genes_j contains only one gene
        
        distance_array[j]<-sqrt(sum((gene_i- genes_j)^2))
      }

      pos=pos+1
    
    }
    
    #save distance between gene_i to all genes in genes_j
    write.table(t(round(distance_array, digits = 8)), sep="\t", append=TRUE, file=outputfile, row.names=FALSE, col.names=FALSE)
    close(connection2)
  }
  
  close(connection1)

}

#computes spearman correlation between gene expression profiles for each gene pair
comp_spearman<-function(inputfile, number_genes, gene_limit, outputfile){
  
  #first connection -> to read gene expression profile of the first gene of each gene pair
  connection1<-file(inputfile)
  open(connection1)
  #skip header + first gene
  gene1<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE,  skip=1, nlines=1)
  
  #select gene expression profile from first gene
  for (i in 2:number_genes){
    
    gene_i<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE, nlines=1)[-1] 
    
    #saves correlation between gene i to gene_1, gene_2, ... gene_i-1
    distance_array<-numeric(length=i-1)
    
    #second connection -> to read gene expression profiles of the second genes of each gene pair
    #(only a limited number can be read because of memory consumption)
    connection2<-file(inputfile)
    open(connection2)
    #skip headline
    genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, skip=1, nlines=gene_limit)
    genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1] #remove column with the GeneIDs
    
    #start with the first gene in genes_j
    pos=1
    
    #select gene expression profile from second gene
    for(j in 1:(i-1)){
      
      #load new gene expression profiles (second genes)
      if(!is.decimal((j-1)/gene_limit) && j!=1){
        
        pos=1
        #a full set of genes can be loaded
        if(i-j>=gene_limit){
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=gene_limit)
          genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1]
        }
        #there are only a few genes to load
        else{
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=i-j)
          genes_j<-matrix(genes_j, nrow=i-j, byrow=TRUE)[,-1]
        }
      }
      
      #genes_j contains more than one gene
      if(is.matrix(genes_j)){
        distance_array[pos]<-cor(gene_i, genes_j[pos,], method="spearman")
        
        
      } else{ #genes_j contains only one gene
        
        distance_array[pos]<-cor(gene_i, genes_j, method="spearman")
      }
      
      pos=pos+1
      
    }
    
    #save correlation between gene_i to all genes in genes_j
    write.table(t(round(distance_array, digits = 8)), sep="\t", append=TRUE, file=outputfile, row.names=FALSE, col.names=FALSE)
    close(connection2)
  }
  
  close(connection1)
  
}


#computes pearson correlation between gene expression profiles for each gene pair
comp_pearson<-function(inputfile, number_genes, gene_limit, outputfile){
  
  #first connection -> to read gene expression profile of the first gene of each gene pair
  connection1<-file(inputfile)
  open(connection1)
  #skip header + first gene
  gene1<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE,  skip=1, nlines=1)
  
  #select gene expression profile from first gene
  for (i in 2:number_genes){
    
    gene_i<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE, nlines=1)[-1] 
    
    #saves correlation between gene i to gene_1, gene_2, ... gene_i-1
    distance_array<-numeric(length=i-1)
    
    #second connection -> to read gene expression profiles of the second genes of each gene pair
    #(only a limited number can be read because of memory consumption)
    connection2<-file(inputfile)
    open(connection2)
    #skip headline
    genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, skip=1, nlines=gene_limit)
    genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1] #remove column with the GeneIDs
    
    #start with the first gene in genes_j
    pos=1
    
    #select gene expression profile from second gene
    for(j in 1:(i-1)){
      
      #load new gene expression profiles (second genes)
      if(!is.decimal((j-1)/gene_limit) && j!=1){
        
        pos=1
        #a full set of genes can be loaded
        if(i-j>=gene_limit){
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=gene_limit)
          genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1]
        }
        #there are only a few genes to load
        else{
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=i-j)
          genes_j<-matrix(genes_j, nrow=i-j, byrow=TRUE)[,-1]
        }
      }
      
      #genes_j contains more than one gene
      if(is.matrix(genes_j)){
        distance_array[pos]<-cor(gene_i, genes_j[pos,], method="pearson")
        
        
      } else{ #genes_j contains only one gene
        
        distance_array[pos]<-cor(gene_i, genes_j, method="pearson")
      }
      
      pos=pos+1
      
    }
    
    #save correlation between gene_i to all genes in genes_j
    write.table(t(round(distance_array, digits = 8)), sep="\t", append=TRUE, file=outputfile, row.names=FALSE, col.names=FALSE)
    close(connection2)
  }
  
  close(connection1)
  
}



#computes RIO (Relative intensity overlap) between gene expression profiles for eachgene pair
comp_rio<-function(inputfile, number_genes, gene_limit, outputfile){
  
  #first connection -> to read gene expression profile of the first gene of each gene pair
  connection1<-file(inputfile)
  open(connection1)
  #skip header + first gene
  gene1<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE,  skip=1, nlines=1)
  
  #select gene expression profile from first gene
  for (i in 2:number_genes){
    
    gene_i<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE, nlines=1)[-1] 
    
    #saves RIO between gene i to gene_1, gene_2, ... gene_i-1
    distance_array<-numeric(length=i-1)
    
    #second connection -> to read gene expression profiles of the second genes of each gene pair
    #(only a limited number can be read because of memory consumption)
    connection2<-file(inputfile)
    open(connection2)
    #skip headline
    genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, skip=1, nlines=gene_limit)
    genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1] #remove column with the GeneIDs
    
    #start with the first gene in genes_j
    pos=1
    
    #select gene expression profile from second gene
    for(j in 1:(i-1)){
      
      #load new gene expression profiles (second genes)
      if(!is.decimal((j-1)/gene_limit) && j!=1){
        
        pos=1
        #a full set of genes can be loaded
        if(i-j>=gene_limit){
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=gene_limit)
          genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1]
        }
        #there are only a few genes to load
        else{
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=i-j)
          genes_j<-matrix(genes_j, nrow=i-j, byrow=TRUE)[,-1]
        }
      }
      
      #genes_j contains more than one gene
      if(is.matrix(genes_j)){
        distance_array[j]<-rio(gene_i, genes_j[pos,])
        
      } else{ #genes_j contains only one gene
        
        distance_array[j]<-rio(gene_i, genes_j)
      }
      
      pos=pos+1
      
    }
    
    #save RIO between gene_i to all genes in genes_j
    write.table(t(round(distance_array, digits = 8)), sep="\t", append=TRUE, file=outputfile, row.names=FALSE, col.names=FALSE)
    close(connection2)
  }
  
  close(connection1)
  
}

#calculates RIO
rio<-function(a, b){
  rio_ab<-sum(a*b)/sum( ( max(abs(a),abs(b)) )^2 )
    return(rio_ab)
}

#computes mutual information between gene expression profiles for each gene pair
comp_mutual_information<-function(inputfile, number_genes, gene_limit, outputfile){
  
  #first connection -> to read gene expression profile of the first gene of each gene pair
  connection1<-file(inputfile)
  open(connection1)
  #skip header + first gene
  gene1<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE,  skip=1, nlines=1)
  
  #select gene expression profile from first gene
  for (i in 2:number_genes){
    
    gene_i<-scan(file=connection1, sep="\t", what=numeric(), quiet=TRUE, nlines=1)[-1] 
    
    #saves mutual information between gene i to gene_1, gene_2, ... gene_i-1
    distance_array<-numeric(length=i-1)
    
    #second connection -> to read gene expression profiles of the second genes of each gene pair
    #(only a limited number can be read because of memory consumption)
    connection2<-file(inputfile)
    open(connection2)
    #skip headline
    genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, skip=1, nlines=gene_limit)
    genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1] #remove column with the GeneIDs
    
    #start with the first gene in genes_j
    pos=1
    
    #select gene expression profile from second gene
    for(j in 1:(i-1)){
      
      #load new gene expression profiles (second genes)
      if(!is.decimal((j-1)/gene_limit) && j!=1){
        
        pos=1
        #a full set of genes can be loaded
        if(i-j>=gene_limit){
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=gene_limit)
          genes_j<-matrix(genes_j, nrow=gene_limit, byrow=TRUE)[,-1]
        }
        #there are only a few genes to load
        else{
          genes_j<-scan(file=connection2, sep="\t", what=numeric(), quiet=TRUE, nlines=i-j)
          genes_j<-matrix(genes_j, nrow=i-j, byrow=TRUE)[,-1]
        }
      }
      
      #genes_j contains more than one gene
      if(is.matrix(genes_j)){
        mutinformation(gene_i, genes_j[pos,], method="emp")
        
      } else{ #genes_j contains only one gene
        
        mutinformation(gene_i, genes_j, method="emp")
      }
      
      pos=pos+1
      
    }
    
    #save mutual information between gene_i to all genes in genes_j
    write.table(t(round(distance_array, digits = 8)), sep="\t", append=TRUE, file=outputfile, row.names=FALSE, col.names=FALSE)
    close(connection2)
  }
  
  close(connection1)
  
}


################################################################################################


#read gene expression data
inputfile<-opt$input

#select GeneIDs (first column of the input file)
com<-paste0("cut -f1 ", inputfile)
gene_ids<-as.numeric(system(command=com, intern=TRUE)[-1])

#save the number of genes
com<-paste0("wc -l ", inputfile, " | awk '{ print $1 }'")
number_genes<-as.numeric(system(command=com, intern=TRUE))-1 #minus headline

##create partial name of the output file
#remove path from the input filename
raw_outputfile<-sub(".*\\/", "", inputfile)
#remove extension from the input filename
raw_outputfile<-sub("(\\.[[:alpha:]]+$)", "", raw_outputfile)

#remark for the output file
remark<-"#File contains lower triangle of the distance or correlation matrix for the following GENEIDs:"

#start stopwatch
ptm <- proc.time()

#the input file is too large to load all gene profiles at once 
#-> load only a group of genes, do the calculations and overwrite this gene group, repeat... 
gene_limit<-100

#method(s) to compare gene expression profiles
distance_measures<-names(which(opt==TRUE))


#calculate distance or correlation between each gene pair
#save the lower triangle of the distance or correlation matrix in the output file
for (method in distance_measures){
  
  #create full output filename
  outputfile<-paste0(output_folder, "/", raw_outputfile, "_distance_matrix_", method,".txt")
  #add remark to the output file
  write.table(t(remark), sep="\t", append=FALSE, quote=FALSE, file=outputfile, row.names=FALSE, col.names=FALSE)
  #save GeneIDs
  write.table(t(gene_ids), sep="\t", append=TRUE, quote=FALSE, file=outputfile, row.names=FALSE, col.names=FALSE)
  
  switch(method, 
         euclidean={
            comp_euclidean(inputfile, number_genes, gene_limit, outputfile)
        },
         spearman={
           comp_spearman(inputfile, number_genes, gene_limit, outputfile)
           
        },
         pearson={
           comp_pearson(inputfile, number_genes, gene_limit, outputfile)    
        },
        rio={
           comp_rio(inputfile, number_genes, gene_limit, outputfile)   
        },
        mutual_information={
          comp_mutual_information(round(data[,-1]), inputfile, number_genes, gene_limit, outputfile) #Werte mussten gerundet werden?   
        }
  )
  
 #stop stopwatch
  print(proc.time() - ptm)
}
