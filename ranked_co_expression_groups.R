#!/usr/bin/env Rscript

########################################################################
##
## File: ranked_co_expression_groups.R
##
## Created by: Anja S.
## Created on: June 08, 2016
##
## Description: Calculates ranked co-expression group for each gene.  
##     
## Input:   File with the GeneIDs and their distance/correlation matrix 
##  (lower triangle).
##
## Output: File with ranked co-expression group for each gene.
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
  
  #which measurement was used to calculate the distance or correlation matrix from the input file?
  make_option(c("-e","--euclidean"), action="store_true", default=FALSE, help="Choose this option if the distace matrix was calculated with the euclidean distance."),
  make_option(c("-s","--spearman"), action="store_true", default=FALSE, help="Choose this option if the distace matrix was calculated with the Spearman correlation coefficient."),
  make_option(c("-p","--pearson"), action="store_true", default=FALSE, help="Choose this option if the distace matrix was calculated with the Pearson correlation coefficient."),
  make_option(c("-r","--rio"), action="store_true", default=FALSE, help="Choose this option if the distace matrix was calculated with the RIO (Relative intensity overlap)."),
  make_option(c("-m","--mutual_information"), action="store_true", default=FALSE, help="Choose this option if the distace matrix was calculated with the mutual information.")
  
)


################################Functions#########################################

#Returns position of gene i and gene j in the distance array (condensed distance matrix)
#grosses Dreieck minus kleines Dreieck
position<-function(i,j,n){return(n*(i-1)+j-(n*(n+1)/2-(n-i+1)*(n-i+2)/2))} 

#anderes kleines Dreieck + Rechteck
position2<-function(i,j,n){return(n*(i-1)+j-( (i-1)*i/2 + (n-i+1)*(i-1) ))}


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

#check if exactly one correlation or distance measure is given
if(sum(opt==TRUE)!=1){
  stop("Choose the measurement which was used to calculate the distance or correlation matrix from the input file!")
}

#user specified an output/result folder
if (!is.null(opt$output)) {
  output_folder<-opt$output
  
}else{ #standard output folder
  output_folder<-"output"
}

#create output folder
dir.create(path=output_folder, recursive=TRUE)


#file with the distance or correlation matrix
inputfile=opt$input

#open a connection
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

#euclidean, spearman, pearson, RIO or mutual information
method<-names(which(opt==TRUE))
#size of the ranked co-expression group
size_rcg<-opt$rcg
print(size_rcg)
number_genes<-length(gene_ids) 
#number_genes<-10


##create name of the output file
#remove path from the inputfilename
raw_outputfile<-sub(".*\\/", "", inputfile)
#remove extension from the inputfilename
raw_outputfile<-sub("(\\.[[:alpha:]]+$)", "", raw_outputfile)
outputfile<-paste0(output_folder, "/", raw_outputfile, "_ranked_co_expression_groups_k_", size_rcg, ".txt")

#write headline to the output file
headline<-c("REFERENCE_GENE")
headline<-append(headline, paste0("GENE",1:size_rcg))
write.table(t(headline), sep="\t", file=outputfile, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)


for(i in 1:number_genes){
  
  rcg_i<-numeric(length=number_genes-1)
  names(rcg_i)<-gene_ids[1:number_genes][-i] #change to gene_ids[-i]
  
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
      rcg_i[pos]<-distance_array[dij]
      pos=pos+1
    }

  }
  
  #sort rcg for gene i and select top "similar" genes
  switch(method, 
    euclidean={  
      rcg_i<-sort(abs(rcg_i), decreasing=FALSE)[1:size_rcg] 
      print("e")
    },
    spearman={   
      rcg_i<-sort(abs(rcg_i), decreasing=TRUE)[1:size_rcg] 
      print("s")
    },
    pearson={    
      rcg_i<-sort(abs(rcg_i), decreasing=TRUE)[1:size_rcg]
      print("p")
    },
    rio={
      rcg_i<-sort(abs(rcg_i), decreasing=TRUE)[1:size_rcg] 
      print("r")
    },
    mutual_information={
      rcg_i<-sort(abs(rcg_i), decreasing=TRUE)[1:size_rcg]    
      print("m")
    }
  )
  
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



