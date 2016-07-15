#!/bin/bash

#calculate for all given RCG files with all given annotation files the number 
#of gene pairs with shared annotation for different k's for each geneID (RCG)

#contains RCG files calculated with different distance/correlation measurements
rcgs_folder=$1

#contains annotation files
annotation_folder=$2

#outputfolder to save the results
output_folder=$3

mkdir ${output_folder}/

g++ shared_annotation_pairs.cpp -o shared_annotation_pairs

for anno_file in ${annotation_folder}/*
do
	tempname=$(basename ${anno_file})
	suboutdir=${tempname%.tsv}

	mkdir ${output_folder}/${suboutdir} #create a suboutputdir for each annotation
	
	for rcg_file in ${rcgs_folder}/*
	do
			
			#determines the number of gene pairs with shared annotation for different k per RCG
			tempname=$(basename ${rcg_file})
			./shared_annotation_pairs -r ${rcg_file} -a ${anno_file} > ${output_folder}/${suboutdir}/shared_annotation_pairs_${tempname%.tsv}.txt

			#Rscript shared_annotation_pairs.R -i ${rcg_file} -o ${output_folder} -a ${anno_file}
	
	done

done








