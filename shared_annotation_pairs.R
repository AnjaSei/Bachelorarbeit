#!/usr/bin/env Rscript

########################################################################
##
## File: shared_annotation_pairs.R
##
## Created by: Anja S.
## Created on: June 22, 2016
##
## Description: Determines the number of gene pairs with shared annotation 
##  per ranked co-expression group (RCG).
##     
## Input:   File with the ranked co-expression groups.
##
## Output: File with number of gene pairs with shared annotation for 
##  different sizes per RCG.
##          
########################################################################

#install CRAN packages
#install.packages("optparse")
#install.packages("data.table")
#install.packages("schoolmath")

#load packages
library(optparse) #to parse command line options
library(data.table)
library(schoolmath)

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
rcgs<-fread(inputfile, data.table=FALSE, header=TRUE)
rcgs<-as.matrix(rcgs)

size_rcg<-max(rcgs[,"rank"])
geneIDs<-unique(rcgs[,"RCG"])

#build a new matrix: each line represents a RCG with its GeneIDs
if(is.decimal(nrow(rcgs)/size_rcg)){
  stop("RCGs have different size!")
}
rcgs<-matrix(rcgs[,"gene"], ncol=size_rcg, byrow=TRUE)
#number of RCGs
number_rcg<-nrow(rcgs)

#only some GeneIDs have an annotation -> set GeneIDs with unknown annotation to NA to save running time
rcgs[which(!(rcgs %in% annotation$GENEID))]<-NA


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

#k is the number of genes with the most correlated 3D gene expression profile to the reference gene
all_k<-c(seq_len(5),seq(10, size_rcg-1, 10))

#saves the number of gene pairs with shared annotation for each k per RCG
shared_anno<-matrix(data=NA, ncol=length(all_k)+1, nrow=number_rcg)
colnames(shared_anno)<-c("RCG", paste0("k=",all_k))
shared_anno[,1]<-geneIDs

#for each ranked co-expression group do:
for(rcg in 1:number_rcg) {

  #load GeneIDs
  rcg_i <- rcgs[rcg,]

  #count gene pairs with shared annotation for each k per RCG 
  counter=0
  for(j in 2:size_rcg){
    for(i in 1:(j-1)){
      
      #select GeneIDs
      geneID1<-rcg_i[i]
      geneID2<-rcg_i[j]
      
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
      
      #save number of shared gene pairs for a specific k
      if( (i==j-1) &&  sum(all_k==j-1) ){
        shared_anno[rcg, which(all_k==j-1)+1]<-counter    
      }
    }
  }

}

#save the results in the outputfile
comment<-paste0("#File with the ranked co-expression groups: ", inputfile, 
                "\n#File with the annotation: ", annotationfile,
                "\n#k=number of most correlated genes to each gene's RCG")
write.table(comment, outputfile, row.names=FALSE, col.names=FALSE, sep="\t", append=FALSE, quote=FALSE )
write.table(shared_anno, outputfile, row.names=FALSE, col.names=TRUE, sep="\t", append=TRUE, quote=FALSE)

