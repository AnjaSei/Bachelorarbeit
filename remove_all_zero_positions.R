#!/usr/bin/env Rscript

########################################################################
##
## File: remove_all_zero_positions.R
##
## Created by: Anja S.
## Created on: June 13, 2016
##
## Description: Creates reduced gene expression profiles by removing
##  the positions which were zero in all original (transformed) gene 
##  expression profiles.
##     
## Input: File with the (haar wavelet transformed) gene expression profiles.
##
## Output: File with the reduced gene expression profiles.
##          
########################################################################

#install CRAN packages
#install.packages("optparse")

#load packages
library(optparse) #to parse command line options

#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help="File with the (transformed) gene expression data."),
  make_option(c("-o","--output"), action="store", help="Result folder (optional).")
)

#parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

#check if the file with the gene expression data is given
if (is.null(opt$input)) {
  stop("No file with the gene expression data is given!")
}

#user specified an output/result folder
if (!is.null(opt$output)) {
  output_folder<-opt$output
  
}else{ #standard output folder
  output_folder<-"output"
}

#create output folder
dir.create(path=output_folder, recursive=TRUE)

#save input file
inputfile=opt$input

#count number of genes
com<-paste0("wc -l ", inputfile, " | awk '{ print $1 }'")
number_genes<-as.numeric(system(command=com, intern=TRUE))-1 #minus headline

#count data points per gene expression profile
com2<-paste0('awk "{print NF; exit}" ', inputfile)
columns<-as.numeric(system(command=com2, intern=TRUE))-1 #minus GeneID

#open a connection to the gene expression file
connection<-file(inputfile)
open(connection)
#read the headline of the gene expression file
headline<-scan(file=connection, nlines=1,  what=character(), quiet=TRUE, sep="\t")


##save for each position of the gene expression profiles whether it is zero for all profiles(-> 0 ) or not (1) 
output_array<-rep(0, columns)

#for each gene do:
for(gene in 1:number_genes) {
  
  #read (haar wavelet transformed) gene expression data
  data <- scan(file=connection, nlines=1,  what=numeric(), quiet=TRUE, sep="\t")
  
  #remove gene id
  input_array<-data[-1] 
  
  #if a position in the input_array is unequal zero then this position in the output_array will be changed to 1/TRUE 
  output_array<- input_array | output_array
  
}
close(connection)

#check which positions in the gene expression profiles are zero for all genes
zero_positions<-which(!output_array)


##create the reduced gene expression profiles

#open a new connection to the original (transformed) gene expression profiles
connection<-file(inputfile)
open(connection)

#update the headline
headline<-scan(file=connection, nlines=1,  what=character(), quiet=TRUE, sep="\t")
headline_reduced<-headline[-(zero_positions+1)]

##create name of the output file
#remove path from the inputfilename
raw_outputfile<-sub(".*\\/", "", inputfile)
#remove extension from the inputfilename
raw_outputfile<-sub("(\\.[[:alpha:]]+$)", "", raw_outputfile)
outputfile<-paste0(output_folder, "/", raw_outputfile, "_without_zero_data_points.txt")

write.table(t(headline_reduced), col.names=FALSE, sep="\t", quote=FALSE, row.names=FALSE, file=outputfile, append=FALSE) 

#read again the gene expression profile for each gene and remove the all-zero positions
for(gene in 1:number_genes) {
  
  #read old (transformed) gene expression profile
  data_old<-scan(file=connection, nlines=1,  what=numeric(), quiet=TRUE, sep="\t") #inclusive GeneID
  
  #remove positions which were zero for all gene expression profiles
  data_reduced<-data_old[-(zero_positions+1)] #zero_positions+1 -> start to count after the GeneID
  
  #save reduced gene expression profile
  write.table(t(data_reduced), col.names=FALSE, sep="\t", append=TRUE, row.names=FALSE, file=outputfile)
  
}
close(connection)
