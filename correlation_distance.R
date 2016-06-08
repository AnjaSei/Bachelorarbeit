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
## Output: File with GeneIDs (first line) and the upper triangle of the 
##  distance/correlation matrix (second line).
##          
########################################################################

#install CRAN packages
#install.packages("optparse")
#install.packages("infotheo")
#install.packages("data.table")

#load packages
library(optparse) #to parse command line options
library(infotheo) #to calculate the mutual information
library(data.table) #fast way to load the data

#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help=".dat file with (e.g. haar wavelet transformed) gene expression data from the Allen Brain Atlas"),
  make_option(c("-o","--output"), action="store", help="result folder (voluntary)"),
  
  #choose a measurement to compare (haar wavelet transformed) gene expression profiles
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

#computes euclidean distance between gene expression profiles for each gene pair
comp_euclidean<-function(data, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-sqrt(sum( (data[reference_gene,][-1] - data[gene,][-1])^2 ))
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}

#computes spearman correlation between gene expression profiles for each gene pair
comp_spearman<-function(data, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-cor(data[reference_gene,][-1], data[gene,][-1], method="spearman")
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}

#computes pearson correlation between gene expression profiles for each gene pair
comp_pearson<-function(data, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-cor(data[reference_gene,][-1], data[gene,][-1], method="pearson")
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}

#computes RIO (Relative intensity overlap) between gene expression profiles for eachgene pair
comp_rio<-function(data, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-rio(data[reference_gene,][-1], data[gene,][-1])
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}

#calculates RIO
rio<-function(a, b){
  rio_ab<-sum(a*b)/sum( ( max(abs(a),abs(b)) )^2 )
    return(rio_ab)
}

#computes mutual information between gene expression profiles for each gene pair
comp_mutual_information<-function(data, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-mutinformation(data[reference_gene,][-1], data[gene,][-1], method="emp")
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}


################################################################################################


#read gene expression data
inputfile<-opt$input
data<-fread(inputfile, header=TRUE, sep="\t", data.table=FALSE)
data<-as.matrix(data)

#method(s) to compare gene expression profiles
distance_measures<-names(which(opt==TRUE))
#save GENEIDs
gene_ids<-data[,1]
#save the number of genes
number_genes<-nrow(data)

##create partial name of the output file
#remove path from the inputfilename
outputfile<-sub(".*\\/", "", inputfile)
#remove extension from the inputfilename
outputfile<-sub("(\\.[[:alnum:]]+$)", "", outputfile)

#remark for the outputfile
remark<-"#The following line contains the distance or correlation between the gene expression profiles for each gene pair. Only the upper triangle of the distance/correlation matrix is saved as a vector.\n#Formula to access the distance between the gene expression profile of the i-th and the j-th gene : number_genes*(i-1)+j-i*(i+1)/2"


#calculate distance or correlation between each gene pair
#save the upper triangle of the distance or correlation matrix in a vector
for (method in distance_measures){
  
  switch(method, 
         euclidean={
           distance_array<-comp_euclidean(data, number_genes)
        },
         spearman={
           distance_array<-comp_spearman(data, number_genes)
           
        },
         pearson={
           distance_array<-comp_pearson(data, number_genes)    
        },
        rio={
           distance_array<-comp_rio(data, number_genes)   
         },
         mutual_information={
           distance_array<-comp_mutual_information(round(data), number_genes) #Werte mussten gerundet werden?   
         }
  )
  
  #save results
  write.table(t(c("GENEID", gene_ids)), sep="\t", append=FALSE, file=paste0(output_folder, "/", outputfile, "_distance_array_",method,".txt"), row.names=FALSE, col.names=FALSE)
  write.table(remark, append=TRUE, file=paste0(output_folder, "/", outputfile, "_distance_array_",method,".txt"), row.names=FALSE, col.names=FALSE)
  write.table(t(round(distance_array, digits = 10)), sep="\t", append=TRUE, file=paste0(output_folder, "/", outputfile, "_distance_array_",method,".txt"), row.names=FALSE, col.names=FALSE)
  
}
