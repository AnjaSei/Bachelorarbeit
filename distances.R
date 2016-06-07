#!/usr/bin/env Rscript

########################################################################
##
## File: distances.R
##
## Created by: Anja S.
## Created on: June 06, 2016
##
## Description: 
##     
## Input: 
##
## Output: 
##          
########################################################################

#install CRAN packages
#install.packages("optparse")
#install.packages("infotheo")
#install.packages("bigmemory")

#load packages
library(optparse) #to parse command line options
library(infotheo) #to calculate the mutual information
library(bigmemory) 


#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help=".dat file with haar wavelet transformed gene expression data from the Allen Brain Atlas"),
  make_option(c("-o","--output"), action="store", help="result folder (voluntary)"),
  
  #choose a measurement to compare haar wavelet transformed gene expression profiles
  make_option(c("-e","--euclidean"), action="store_true", default=FALSE, help="Compare haar wavelet transformed gene expression profiles with euclidean distance"),
  make_option(c("-s","--spearman"), action="store_true", default=FALSE, help="Compare haar wavelet transformed gene expression profiles with spearman correlation coefficient"),
  make_option(c("-p","--pearson"), action="store_true", default=FALSE, help="Compare haar wavelet transformed gene expression profiles with pearson correlation coefficient"),
  make_option(c("-r","--rio"), action="store_true", default=FALSE, help="Compare haar wavelet transformed gene expression profiles with RIO (Relative intensity overlap)"),
  make_option(c("-m","--mutual_information"), action="store_true", default=FALSE, help="Compare haar wavelet transformed gene expression profiles with mutual information")
  
)

#parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

#check if file with gene expression data is given
if (is.null(opt$input)) {
  stop("No file with haar wavelet transformed gene expression data is given!")
}

#check if a measurement (to compare haar wavelet transformed gene expression profiles) has been selected
if(!any(opt$euclidean, opt$spearman, opt$pearson, opt$rio, opt$mutual_information)){
  stop("Select a measurement to compare haar wavelet transformed gene expression profiles!")
}
#allow only one measurement method
if(sum(c(opt$euclidean, opt$spearman, opt$pearson, opt$rio, opt$mutual_information))>1){
  stop("Choose only ONE measurement to compare haar wavelet transformed gene expression profiles!")
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

comp_euclidean<-function(data_transformed, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-sqrt(sum( (data_transformed[reference_gene,][-1] - data_transformed[gene,][-1])^2 ))
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}


comp_spearman<-function(data_transformed, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-cor(data_transformed[reference_gene,][-1], data_transformed[gene,][-1], method="spearman")
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}

comp_pearson<-function(data_transformed, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-cor(data_transformed[reference_gene,][-1], data_transformed[gene,][-1], method="pearson")
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}


comp_rio<-function(data_transformed, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-rio(data_transformed[reference_gene,][-1], data_transformed[gene,][-1])
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}


rio<-function(a, b){
  rio_ab<-sum(a*b)/sum( ( max(abs(a),abs(b)) )^2 )
    return(rio_ab)
}

comp_mutual_information<-function(data_transformed, number_genes){
  
  distance_array<-numeric(length=number_genes*(number_genes-1)/2)
  
  pos=1
  for (reference_gene in 1:(number_genes-1)){
    for(gene in (reference_gene+1):number_genes){
      distance_array[pos]<-mutinformation(data_transformed[reference_gene,][-1], data_transformed[gene,][-1], method="emp")
      pos=pos+1
      
    }
  }
  
  return(distance_array)
  
}


################################################################################################


##read haar wavelet transformed data
#file with the gene expression data
inputfile=opt$input
#inputfile<-"/home/anja/Schreibtisch/Daten/final/output/haar_wavelet_transformed_data_cube_edge_length_16.txt"
data_transformed<-read.big.matrix(inputfile, header=TRUE, sep="\t", type="double")

#how to compare gene expression profiles
distance_measures<-names(which(opt==TRUE))
number_genes<-nrow(data_transformed)
#number_genes<-10


for (method in distance_measures){
  
  switch(method, 
         euclidean={
           distance_array<-comp_euclidean(data_transformed, number_genes)
           write.table(t(round(distance_array, digits = 8)), col.names=FALSE, sep="\t", append=FALSE, row.names=FALSE, file=paste0(output_folder,strsplit(inputfile, ".txt"),"distance_array_",method,".txt"))
           
         },
         spearman={
           distance_array<-comp_spearman(data_transformed, number_genes)   
           write.table(t(round(distance_array, digits = 8)), col.names=FALSE, sep="\t", append=TRUE, row.names=FALSE, file=paste0(output_folder,strsplit(inputfile, ".txt"),"distance_array_",method,".txt"))
           
           
         },
         pearson={
           distance_array<-comp_pearson(data_transformed, number_genes)    
           write.table(t(round(distance_array, digits = 8)), col.names=FALSE, sep="\t", append=TRUE, row.names=FALSE, file=paste0(output_folder,strsplit(inputfile, ".txt"),"distance_array_", method, ".txt"))
         },
         rio={
           distance_array<-comp_rio(data_transformed, number_genes)   
           write.table(t(round(distance_array, digits = 8)), col.names=FALSE, sep="\t", append=TRUE, row.names=FALSE, file=paste0(output_folder,strsplit(inputfile, ".txt"),"distance_array_", method, ".txt"))
           
         },
         mutual_information={
           distance_array<-comp_mutual_information(round(data_transformed), number_genes) #Werte mussten gerundet werden?   
           write.table(t(round(distance_array, digits = 8)), col.names=FALSE, sep="\t", append=TRUE, row.names=FALSE, file=paste0(output_folder,strsplit(inputfile, ".txt"),"distance_array_", method, ".txt"))
         }
  )
}
