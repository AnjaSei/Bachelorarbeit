#!/usr/bin/env Rscript

########################################################################
##
## File: shared_annotation.R
##
## Created by: Anja S.
## Created on: June 12, 2016
##
## Description: Determines the number of gene pairs which share an annotation 
##  in each ranked co-expression group (RCG).
##     
## Input:   File with the ranked co-expression groups.
##
## Output: File which contains number of gene pairs with shared 
##  annotation for each ranked co-expression group.
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
  make_option(c("-i", "--input"), action="store", help="File with the ranked co-expression groups."),
  make_option(c("-o","--output"), action="store", help="Result folder (optional)."),
  
  make_option(c("-a", "--annotation"), action="store", help="File with the GeneIDs (first column) and their annotation (second column).")
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
inputfile<-opt$input
rcg<-fread(inputfile, header=TRUE, sep="\t", data.table=FALSE)
rcg<-as.matrix(rcg)
ref_genes<-rcg[,1]

#only around 5000 GeneIDs have an annotation -> set GeneIDs without known annotation to NA
#to save running time (later in the code)
rcg[which(!(rcg %in% annotation$GENEID))]<-NA

#number genes per RCG
size_rcg<-ncol(rcg)
#number of RCGs
number_rcg<-nrow(rcg)

##create name of the output file
#remove path from the RCG-filename
rcg_filename<-sub(".*\\/", "", inputfile)
#remove extension from the RCG-filename
rcg_filename<-sub("(\\.[[:alpha:]]+$)", "", rcg_filename)
#remove path from the annotation filename
annotationfilename<-sub(".*\\/", "", annotationfile)
#remove extension from theannotation filename
annotationfilename<-sub("(\\.[[:alpha:]]+$)", "", annotationfilename)
outputfile<-paste0(output_folder, "/", rcg_filename, "_annotation_", annotationfilename, ".txt")

#saves the number of gene pairs with shared annotation for each RCG
shared_anno<-numeric(length=number_rcg)

#for each ranked co-expression group do:
for(i in 1:number_rcg) {
  
  #load GeneIDs for each ranked co-expression group
  rcg_i <- rcg[i,]
  
  #count gene pairs with shared annotation per RCG
  counter=0
  for(gene1 in 1:(size_rcg-1)){
    for(gene2 in (gene1+1):size_rcg){
      
      #select GeneIDs
      geneID1<-rcg_i[gene1]
      geneID2<-rcg_i[gene2]
      
      #both GeneIDs have an annotation 
      if(!is.na(geneID1) && !is.na(geneID2)){
        
        #select annotations
        anno_gene1<-annotation[which(annotation$GENEID==geneID1),]$ANNOTATION 
        anno_gene2<-annotation[which(annotation$GENEID==geneID2),]$ANNOTATION  
        
        #gene pair share at least one annotation
        if(any(anno_gene1 %in% anno_gene2)){
          counter=counter+1
        }
        
      }
    }
    
  }
  #save the number of gene pairs with shared annotation
  shared_anno[i]<-counter

}

#save results in the outputfile
result_matrix<-cbind(ref_genes, shared_anno)
colnames(result_matrix)<-c("ranked_co_expression_group", "number_shared_annotation_pairs")
write.table(result_matrix, outputfile, row.names=FALSE, col.names=TRUE, sep="\t", append=FALSE, quote=FALSE)

