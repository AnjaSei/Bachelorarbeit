#!/usr/bin/env Rscript

########################################################################
##
## File: ranked_co_expression_groups.R
##
## Created by: Anja S.
## Created on: June 08, 2016
##
## Description: Calculates a ranked co-expression group for each gene.  
##     
## Input:   File with the GeneIDs and their distance/correlation matrix 
##  (lower triangle).
##
## Output: File with a ranked co-expression group for each gene.
##          
########################################################################

#install CRAN packages
#install.packages("optparse")

#load packages
library(optparse) #to parse command line options

#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help="File with the GeneIDs and their distance/correlation matrix (lower triangle)."),
  make_option(c("-o","--output"), action="store", help="Result folder (voluntary)"),
  
  make_option(c("-k", "--rcg"), action="store", type="integer", help="Number of genes per ranked co-expression group (RCG); without the gene that defines the RCG"),
  make_option(c("-d","--dist"), action="store", help="Enter the distance or correlation measurement which was used to calculate the distance/correlation matrix. Available options are: euclidean, spearman, pearson, rio and mutual_information.")
)


################################Functions#########################################

#returns position of gene i and gene j in the distance array (condensed distance matrix)
position<-function(i,j,n){return(n*(i-1)+j-(n*(n+1)/2-(n-i+1)*(n-i+2)/2))} 
#position2<-function(i,j,n){return(n*(i-1)+j-( (i-1)*i/2 + (n-i+1)*(i-1) ))}

#calculates RCG for gene i
rcg_increasing_abs<-function(rcg_i, size_rcg){
  return(sort(abs(rcg_i), decreasing=FALSE)[1:size_rcg])
}

rcg_decreasing<-function(rcg_i, size_rcg){
  return(sort(rcg_i, decreasing=TRUE)[1:size_rcg])
}

##################################################################################

#parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

#check if file with the distance or correlation matrix is given
if (is.null(opt$input)) {
  stop("No file with the distance or correlation matrix is given!")
}

#check if the size of the ranked co-expression group is given
if(is.null(opt$rcg)){
  stop("Enter the size of the ranked co-expression group!")
}

#check if a correlation or distance measure is given
if(is.null(opt$dist)){
  stop("Choose the measurement which was used to calculate the distance or correlation matrix from the input file!")
}

#check if the user selected measurement is available
if(!(opt$dist %in% c("euclidean", "spearman", "pearson", "rio", "mutual_information"))){
  stop("Only euclidean, spearman, pearson, rio and mutual_information are available measurements!")
}

#user specified an output/result folder
if (!is.null(opt$output)) {
  output_folder<-opt$output
  
}else{ #standard output folder
  output_folder<-"output"
}

#create output folder
dir.create(path=output_folder, recursive=TRUE)

#read file with the distance or correlation matrix
inputfile=opt$input
connection<-file(inputfile)
open(connection)
gene_ids<-scan(file=connection,  what=numeric(), quiet=TRUE, sep="\t",  nlines=1, skip=1) 
distance_array<-scan(file=connection, what=numeric(), quiet=TRUE, sep="\t")
close(connection)


###distance or correlation matrix from the input file:     
##Gene j    a b c ... 
##Gene i  
##    a    1
##    b    2 3 
##    c    4 5 6                  
##   ...                          -> distance_array=c(1,2,3,4,5,6 ...)

#euclidean, spearman, pearson, rio or mutual_information
method<-opt$dist

#size of the ranked co-expression group
size_rcg<-opt$rcg

number_genes<-length(gene_ids) 

##create name of the output file
#remove path from the inputfilename
raw_outputfile<-sub(".*\\/", "", inputfile)
#remove extension from the inputfilename
raw_outputfile<-sub("(\\.[[:alpha:]]+$)", "", raw_outputfile)
outputfile<-paste0(output_folder, "/", raw_outputfile, "_ranked_co_expression_groups_k_", size_rcg, ".txt")

#write headline to the output file
headline<-c("REFERENCE_GENE")
headline<-append(headline, paste0("GENE", 1:size_rcg))
write.table(t(headline), sep="\t", file=outputfile, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)

#decide how to calculate RCGs
switch(method, 
       euclidean={  
         calc_rcg<-rcg_increasing_abs 
       },
       spearman={   
         calc_rcg<-rcg_decreasing
       },
       pearson={    
         calc_rcg<-rcg_decreasing
       },
       rio={
         calc_rcg<-rcg_decreasing
       },
       mutual_information={
         calc_rcg<-rcg_decreasing   
       }
)

#calculate RCG for each gene
for(i in 1:number_genes){
  
  group_i<-numeric(length=number_genes-1)
  names(group_i)<-gene_ids[-i]
  
  pos=1
  
  #select distance or correlation between gene i and gene j
  for(j in 1:number_genes){
    
    if(i!=j){
      if(i>j){
        dij<-position(i, j, number_genes)
      }
      else{
        dij<-position(j, i, number_genes)
      }
      group_i[pos]<-distance_array[dij]
      pos=pos+1
    }

  }
  
  rcg_i<-calc_rcg(group_i, size_rcg)
  
  write.table(t(c(gene_ids[i], names(rcg_i))), sep="\t", file=outputfile, col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)
  #write.table(t(round(rcg_i, digits=10)), sep="\t", file=outputfile, col.names=TRUE, row.names=FALSE, append=TRUE)
}



#spearman
#moeglichst nah an Eins

#pearson
#moeglichst nah an Eins

#RIO
#moeglichst nah an Eins

#mutual information
#kann nicht negativ sein -> Maximum gesucht



