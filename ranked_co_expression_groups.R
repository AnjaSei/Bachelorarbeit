#!/usr/bin/env Rscript

########################################################################
##
## File: ranked_co_expression_groups.R
##
## Created by: Anja S.
## Created on: June 08, 2016
##
## Description:   
##     
## Input:   File with the GeneIDs and their distance/correlation matrix 
##  (lower triangle)
##
## Output: 
##          
########################################################################

#install CRAN packages
#install.packages("optparse")

#load packages
library(optparse) #to parse command line options

#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help=""),
  make_option(c("-o","--output"), action="store", help="result folder (voluntary)"),
  
  #number of genes per ranked co-expression group (RCG); without the gene that defines the RCG
  make_option(c("-k", "--rcg"), action="store_true", default=FALSE, help="Compare gene expression profiles with euclidean distance")
)

#parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

#check if file with the distance or correlation array is given
if (is.null(opt$input)) {
  stop("No file with the distance or correlation array is given!")
}

#check if the size of the ranked co-expression group is given
if(is.null(opt$rcg)){
  stop("Enter the size of the ranked co-expression group!")
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

###distance or correlation matrix from the input file:     
##Gene j    a b c ... 
##Gene i  
##    a    1
##    b    2 3 
##    c    4 5 6                  
##   ...                          -> distance_array=c(1,2,3,4,5,6 ...)

#open a connection
connection<-file(inputfile)
open(connection)
gene_ids<-scan(file=connection,  what=numeric(), quiet=TRUE, sep="\t",  nlines=1, skip=1)
distance_array<-scan(file=connection, what=numeric(), quiet=TRUE, sep="\t")
close(connection)

#search for a measurement
method<-sub(".*(euclidean|spearman|pearson|rio|mutual_information).*", "\\1", outputfile)

##create partial name of the output file
#remove path from the inputfilename
outputfile<-sub(".*\\/", "", inputfile)
#remove extension from the inputfilename
outputfile<-sub("(\\.[[:alpha:]]+$)", "", outputfile)


#grosses Dreieck minus kleines Dreieck
position<-function(i,j,n){return(n*(i-1)+j-(n*(n+1)/2-(n-i+1)*(n-i+2)/2))} 
#anderes kleines Dreieck + Rechteck
position2<-function(i,j,n){return(n*(i-1)+j-( (i-1)*i/2 + (n-i+1)*(i-1) ))}
  
#euclidean
#moeglichst nah an Null


number_genes<-length(gene_ids) #evtl. "GENEIDs" abziehen
size_rcg<-opt$k

outputfile="/home/anja/Schreibtisch/bla.txt"

headline<-c("REFERENCE_GENE")
headline<-append(headline, paste0("GENE",1:size_rcg))
write.table(t(headline), sep="\t", file=outputfile, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)


for(i in 1:number_genes){
  
  rcg_i<-numeric(length=number_genes-1)
  names(rcg_i)<-gene_ids[1:number_genes][-i] #change to gene_ids[-i]
  pos=1
  
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
    },
    spearman={   
      rcg_i<-sort(abs(rcg_i), decreasing=TRUE)[1:size_rcg] 
    },
    pearson={    
      rcg_i<-sort(abs(rcg_i), decreasing=TRUE)[1:size_rcg]
    },
    rio={
      rcg_i<-sort(abs(rcg_i), decreasing=TRUE)[1:size_rcg] 
    },
    mutual_information={
      rcg_i<-sort(abs(rcg_i), decreasing=TRUE)[1:size_rcg]    
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

