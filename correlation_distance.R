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
## Output: File with the GeneIDs and the lower triangle of their 
##  distance/correlation matrix.
##          
########################################################################

#install CRAN packages
#install.packages("optparse")
#install.packages("DescTools")
#install.packages("data.table")

#load packages
library(optparse) #to parse command line options
library(DescTools) #to calculate the mutual information
library(data.table) #provides fast way to load the data

#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help=".dat file with (e.g. haar wavelet transformed) gene expression profiles from the Allen Brain Atlas."),
  make_option(c("-o","--output"), action="store", help="Result folder (optional)."),
  
  #choose a measurement to compare (haar wavelet transformed) gene expression profiles
  make_option(c("-e","--euclidean"), action="store_true", default=FALSE, help="Compare gene expression profiles with euclidean distance."),
  make_option(c("-s","--spearman"), action="store_true", default=FALSE, help="Compare gene expression profiles with spearman correlation coefficient."),
  make_option(c("-p","--pearson"), action="store_true", default=FALSE, help="Compare gene expression profiles with pearson correlation coefficient."),
  make_option(c("-r","--rio"), action="store_true", default=FALSE, help="Compare gene expression profiles with RIO (Relative intensity overlap.)"),
  make_option(c("-m","--mutual_information"), action="store_true", default=FALSE, help="Compare gene expression profiles with mutual information.")
  
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

#computes euclidean distance between two gene expression profiles
comp_euclidean<-function(gene_i, gene_j){

  dist<-sqrt(sum((gene_i- gene_j)^2))
  return(dist)
}

#computes spearman correlation between two gene expression profiles
comp_spearman<-function(gene_i, gene_j){
  
  correlation<-cor(gene_i, gene_j, method="spearman")
  return(correlation)
}

#computes pearson correlation between two gene expression profiles
comp_pearson<-function(gene_i, gene_j){
  
  correlation<-cor(gene_i, gene_j, method="pearson")
  return(correlation)
}

#calculates RIO between two gene expression profiles
comp_rio<-function(gene_i, gene_j){
  
  rio_ij<-sum(gene_i*gene_j)/sum( apply(rbind(abs(gene_i), abs(gene_j)), MARGIN=2, max)^2 )
  return(rio_ij)
}

#computes mutual information between two gene expression profiles
comp_mutual_information<-function(gene_i, gene_j){
  
      mut<-MutInf(gene_i, gene_j)
      return(mut)
}


################################################################################################


#gene expression data
inputfile<-opt$input

#load gene expression profiles
data<-fread(inputfile, header=TRUE, sep="\t")
data<-as.matrix(data)

#method(s) to compare gene expression profiles
distance_measures<-names(which(opt==TRUE))

#save GENEIDs
gene_ids<-data[,1]

#save the number of genes
number_genes<-nrow(data)

##create partial name of the output file
#remove path from the input filename
raw_outputfile<-sub(".*\\/", "", inputfile)
#remove extension from the input filename
raw_outputfile<-sub("(\\.[[:alpha:]]+$)", "", raw_outputfile)

#remark for the output file
remark<-"#File contains lower triangle of the distance or correlation matrix for the following GENEIDs:"

#calculate distance or correlation matrix for each distance measure
for (method in distance_measures){
  
  switch(method, 
         euclidean={
           distance_function<-comp_euclidean
         },
         spearman={
           distance_function<-comp_spearman
           
         },
         pearson={
           distance_function<-comp_pearson 
         },
         rio={
           distance_function<-comp_rio
         },
         mutual_information={
           distance_function<-comp_mutual_information   
         }
  )
  
  #create full output filename
  outputfile<-paste0(output_folder, "/", raw_outputfile, "_distance_matrix_", method,".txt")
  #add remark to the output file
  write.table(t(remark), sep="\t", append=FALSE, quote=FALSE, file=outputfile, row.names=FALSE, col.names=FALSE)
  #save GeneIDs
  write.table(t(gene_ids), sep="\t", append=TRUE, quote=FALSE, file=outputfile, row.names=FALSE, col.names=FALSE)
  
  #calculate distance or correlation between each gene pair (gene i, gene j)
  #save the lower triangle of the distance or correlation matrix in the output file
  for (i in 2:number_genes){

    #load gene expression profile from gene i
    gene_i<-data[i,][-1] #remove column with the GeneIDs
      
    #saves distance between gene i to gene_1, gene_2, ... gene_i-1
    distance_array<-numeric(length=i-1)
      
    #select gene expression profile from second gene
    for(j in 1:(i-1)){
      
      gene_j<-data[j,][-1]  #remove column with the GeneIDs
      distance_array[j]<-distance_function(gene_i, gene_j)
        
    }
      
    #save all distances to gene_i 
    write.table(t(round(distance_array, digits = 8)), sep="\t", append=TRUE, file=outputfile, row.names=FALSE, col.names=FALSE)
  }
    
}
