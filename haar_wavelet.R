#!/usr/bin/env Rscript

########################################################################
##
## File: haar_wavelet.R
##
## Created by: Anja S.
## Created on: May 30, 2016
##
## Description: Performes Haar wavelet transformation on gene 
##  expression data from the Allen Brain Atlas
##     
## Input: .dat with gene expression data
##  Path to output folder (optional)
##  Cube edge length (optional)
##
## Output: .txt with Haar wavelet transformed gene expression data
##          
########################################################################

#install CRAN packages
#install.packages("waveslim")
#install.packages("optparse")
#install.packages("schoolmath")

#load packages
library(waveslim) #haar wavelet transformation
library(optparse) #to parse command line options
library(schoolmath)


#command line options
option_list <- list(
  make_option(c("-i", "--input"), action="store", help=".dat file with gene expression data from the Allen Brain Atlas"),
  make_option(c("-o","--output"), action="store", help="result folder (voluntary)"),
  make_option(c("-c","--cube_edge_length"), action="store", type="integer", help="user specified cube  edge length; must be a multiple of 2 (voluntary)")
)

#parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

#check if file with gene expression data is given
if (is.null(opt$input)) {
  stop("No file with gene expression data is given!")
}

#if cube edge length was specified by the user: check if it is a multiple of 2
if(!is.null(opt$cube_edge_length)){
  if(is.decimal(log2(opt$cube_edge_length))){
    stop("Cube edge length have to be a multiple of 2!")
  }
}

#user specified an output/result folder
if (!is.null(opt$output)) {
  output_folder<-opt$output
  
}else{ #standard output folder
    output_folder<-"output"
}

#create output folder
dir.create(path=output_folder, recursive=TRUE)



##########Functions##########

#fill given gene expression data (saved in intensities) from one gene centrally into a 3D array (brain)
fill_brain<-function(intensities, coordinates, cube_edge_length, x_shift, y_shift, z_shift, brain){
  
  #position in intensities
  pos=1
  
  #Insert intensities of the "old" cuboid as centrally as possible in the "new" cube, move coordinates accordingly
  for (coord in 1:ncol(coordinates)){
    
    #move data points to the center of the "new cube"
    x<-coordinates[1,coord]+x_shift
    y<-coordinates[2,coord]+y_shift
    z<-coordinates[3,coord]+z_shift
    
    #coordinates are not cut off (can happen when the user selects the cube edge length )
    if(all(c(x,y,z)>=1) && all(c(x,y,z)<=cube_edge_length) ){
      
      #save intensities from 1D array in 3D array (cube)
      brain[x,y,z]<-intensities[pos]
      pos<-pos+1
    }
  }
  
  return(brain)
  
}



#############################

#file with the gene expression data
inputfile=opt$input

#open a connection
connection<-file(inputfile)
open(connection)

#count number of lines of the gene expression file
com<-paste0("wc -l ", inputfile, " | awk '{ print $1 }'")
row_number<-as.numeric(system(command=com, intern=TRUE))

##save coordinates  of the gene expression data in a matrix 
#(first row -> x-coordinates, second row -> y-coordinates, third row -> z-coordinates)
headline <- scan(file=connection, nlines=1,  what=character(), quiet=TRUE, sep="\t")
#split each coordinate to x-, y- and z-value
coordinates<-strsplit(headline, split=",")[-1]
#convert from character to numeric
coordinates<-sapply(coordinates, as.numeric)


#user specified cube edge length
if(!is.null(opt$cube_edge_length)){
  cube_edge_length<-opt$cube_edge_length
  
} else { #calculate cube edge length
  
  #calculate range of the x-, y- resp. z- coordinates
  #it is  okay to consider only the first and the last coordinate because only whole data planes 
  #(with only zero-intensities) and no single coordinates are missing in the gene expression data set
  
  range_x<-coordinates[1, ncol(coordinates)]-coordinates[1,1]+1 #==range_y<-max(coordinates[1,])-min(coordinates[1,])+1
  range_y<-coordinates[2, ncol(coordinates)]-coordinates[2,1]+1 #==range_y<-max(coordinates[2,])-min(coordinates[2,])+1
  range_z<-coordinates[3, ncol(coordinates)]-coordinates[3,1]+1 #==range_z<-max(coordinates[3,])-min(coordinates[3,])+1
  
  cube_edge_length<-max(range_x,range_y, range_z)
  
  #cube_edge_length is no multiple of 2
  if(is.decimal(log2(cube_edge_length))){
    
    #widen cube_edge_length to a multiple of 2
    cube_edge_length<-2^ceiling(log2(cube_edge_length))
    
  }
}

#number of decompositions of the haar wavelet transformation
decomposition_depth=log2(cube_edge_length)

#center of the given gene expression coordinates, considered as a cuboid
#(min_x+max_x)/2;  m_x_old=(min(coordinates[1,])+max(coordinates[1,]))/2
m_x_old<-(coordinates[1,1]+coordinates[1,ncol(coordinates)])/2
#(min_y+max_y)/2;  m_y_old=(min(coordinates[2,])+max(coordinates[2,]))/2
m_y_old<-(coordinates[2,1]+coordinates[2,ncol(coordinates)])/2
#(min_z+max_z)/2;  m_z_old=(min(coordinates[3,])+max(coordinates[3,]))/2
m_z_old<-(coordinates[3,1]+coordinates[3,ncol(coordinates)])/2

#center of the "new" cube (->brain)
m_x=m_y=m_z=cube_edge_length/2

#shift data points of the "old" cuboid along the x-, y- or z-axis
x_shift<-ceiling(m_x-m_x_old)
y_shift<-ceiling(m_y-m_y_old)
z_shift<-ceiling(m_z-m_z_old)

#"empty" 3D brain (filled with nulls)
brain<-array(data=0, dim=c(cube_edge_length, cube_edge_length, cube_edge_length))


#name of the result file (stores haar wavelet transformed expression data)
resultfile=paste0(output_folder,"/haar_wavelet_transformed_data_cube_edge_length_",cube_edge_length,".txt")

#build a line for the column names of the resultfile
column_names<-c("GENEID")
filter<-c("HLL","LHL","LLH","HHL","HLH","LHH","HHH")
data_points<-cube_edge_length
decomposition_step=1

repeat{
  data_points<-data_points/2
  for (filtertyp in filter){
    column_names<-c(column_names, rep(paste0(filtertyp, as.character(decomposition_step)), data_points^3))
  }
  if(data_points==1){
    column_names<-c(column_names, paste0("LLL",decomposition_step))
    break
  }
  decomposition_step<-decomposition_step+1
  
}
write.table(t(column_names), col.names=FALSE, sep="\t", row.names=FALSE, file=resultfile)


#for each gene do:
for(gene in 1:(row_number-1)) {

  #read gene_id and intensities for each coordinate of the gene
  data <- scan(file=connection, nlines=1,  what=numeric(), quiet=TRUE, sep="\t")
  
  #select geneID
  gene_id<-data[1]
  
  #intensities for one gene at each coordinate (without gen_id), planes with only zero-intensities are missing 
  intensities<-data[-1]
  
  #create 3D brain (cube) with the intensities of the gene at the specific coordinates
  brain<-fill_brain(intensities, coordinates, cube_edge_length, x_shift, y_shift, z_shift, brain)
  
  #apply haar wavelet transformation
  data_transform<-dwt.3d(brain, wf="haar", J=decomposition_depth)
  
  #concatenate the data_transform-list to a vector
  data_transform_vector<-Reduce(c,data_transform)
  
  #write transformed gene intensity values to the result file
  write.table(t(c(gene_id,round(data_transform_vector, digits = 8))), col.names=FALSE, sep="\t", append=TRUE, row.names=FALSE, file=resultfile)

}
close(connection)
