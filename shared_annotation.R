#!/usr/bin/env Rscript

#cut -f1 c2.cp.kegg.v5.0.entrez.gene2anno-nodisease.tsv | sort -n | uniq | wc -l
#5034
# cut -f1 c5.cc.v5.0.entrez.gene2anno.tsv | sort -n | uniq | wc -l
#5270

#awk '{s+=$2} END {print s}' test.txt

########################################################################
##
## File: shared_annotation.R
##
## Created by: Anja S.
## Created on: June 12, 2016
##
## Description: 
##     
## Input:   File with the ranked co-expression groups.
##
## Output: 
##          
########################################################################

#install CRAN packages
#install.packages("optparse")
#install.packages("data.table")

#load packages
library(optparse) #to parse command line options
library(data.table)

#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help="File with ranked co-expression groups."),
  make_option(c("-o","--output"), action="store", help="result folder (voluntary)"),
  
  make_option(c("-a", "--annotation"), action="store", help="File with GeneIDs (first column) and annotation (second column).")
  
)

#parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

#check if file with the ranked co-expression groups is given
if (is.null(opt$input)) {
  stop("No file with the ranked co-expression groups is given!")
}

#check if file with the annotation is given
if (is.null(opt$annotation)) {
  stop("No annotation file is given!")
}

#user specified an output/result folder
if (!is.null(opt$output)) {
  output_folder<-opt$output
  
}else{ #standard output folder
  output_folder<-"output"
}

#create output folder
dir.create(path=output_folder, recursive=TRUE)



#load the annotation file
annotationfile<-opt$annotation
annotation<-fread(annotationfile, header=FALSE, sep="\t", data.table=FALSE)
colnames(annotation)<-c("GENEID", "ANNOTATION")

#load file with the ranked co-expression groups
inputfile=opt$input
#open a connection
connection<-file(inputfile)
open(connection)

#count number of ranked co-expression groups
com<-paste0("wc -l ", inputfile, " | awk '{ print $1 }'")
number_rcg<-as.numeric(system(command=com, intern=TRUE))-1 #minus headline

headline_inputfile <- scan(file=connection, nlines=1,  what=character(), quiet=TRUE, sep="\t")
size_rcg<-length(headline_inputfile)

##create name of the output file
#remove path from the inputfilename
rcg_filename<-sub(".*\\/", "", inputfile)
#remove extension from the inputfilename
rcg_filename<-sub("(\\.[[:alpha:]]+$)", "", rcg_filename)
#remove path from the annotation filename
annotationfilename<-sub(".*\\/", "", annotationfile)
#remove extension from the inputfilename
annotationfilename<-sub("(\\.[[:alpha:]]+$)", "", annotationfilename)
outputfile<-paste0(output_folder, "/", rcg_filename, "_annotation_", annotationfilename, ".txt")


#write headline to the resultfile
headline_resultfile<-c("#ranked_co_expression_group", "number_shared_annotation_pairs")
write.table(t(headline_resultfile), outputfile, row.names=FALSE, col.names=FALSE, sep="\t", append=FALSE, quote=FALSE)

#for each ranked co-expression group do:
for(i in 1:number_rcg) {
  
  #load GeneIDs for each ranked co-expression group
  rcg <- scan(file=connection, nlines=1,  what=numeric(), quiet=TRUE, sep="\t")
  
  #count gene pairs with shared annotation per rcg
  counter=0
  for(gene1 in 1:(size_rcg-1)){
    for(gene2 in (gene1+1):size_rcg){
      geneID1<-rcg[gene1]
      geneID2<-rcg[gene2]
      #print(paste0(geneID1, "|", geneID2))
      #select annotations for both GeneIDs
      anno_gene1<-subset(annotation$ANNOTATION, annotation$GENEID==geneID1) #or annotation[annotation$GENEID==gene1,]$ANNOTATION
      anno_gene2<-subset(annotation$ANNOTATION, annotation$GENEID==geneID2) #or annotation[annotation$GENEID==gene2,]$ANNOTATION

      #gene pair share at least one annotation
      if(sum(anno_gene1 %in% anno_gene2)>=1){
       counter=counter+1
        #print(paste0(anno_gene1, "|", anno_gene2))
      }
    }
  }
  
  #save results
  write.table(t(c(rcg[1], counter)), outputfile, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE)

}


