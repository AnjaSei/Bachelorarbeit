#!/bin/bash

#calculates ranked co-expression groups with k between 1 and 100 for all given distance/correlation matrices
r_script=$1  		#ranked_co_expression_groups.R
input_folder=$2		#contains files each with a distance/correlation matrix 
output_folder=$3	#stores the generated files (RCGs)


#call ranked_co_expression_groups.R to calculate RCGs with k=100 for all given distance/correlation matrices
#for distance/correlation matrix do:
for dist_file in ${input_folder}/*.{txt,gz}
do 
	#echo $dist_file	
	fname=$(basename "$dist_file")
	dist="unknown"

	#extract the distance/correlation measure
	case $fname in
		*euclidean*) 
			dist="euclidean"
		;;
		*spearman*) 
			dist="spearman"
		;;
		*pearson*) 
			dist="pearson"
		;;
		*rio*) 
			dist="rio"
		;;
		*mutual_information*) 
			dist="mutual_information"
		;;
	esac
	

	if [ "$dist" != "unknown" ]
	    then
		echo $dist
		outputdir=${output_folder}/${dist}/
		#Rscript $r_script -i $dist_file -o $outputdir -k 100 -d $dist

		#after generating RCGs with k=100 -> determine remaining RCGs between 1 and 90
		all_k=(1 2 3 5 10 20 30 40 50 60 70 80 90)
		for k in "${all_k[@]}"
		do
			#name of the old RCG file with k=100
		 	rcg_k_100=$( basename ${outputdir}/*k_100.txt)
			
		 	#name of the new RCG file with smaller k
		 	rcg_new=${rcg_k_100%k*}k_$k.txt
			
		 	#select first k+1 columns from the old RCG file (k=100) to build smaller RCGs
		 	cut -f$(seq -s, 1 1 $((k+1))) ${outputdir}/${rcg_k_100} > ${outputdir}/${rcg_new}
		done
	    else
		echo "Unknown distance/correlation measure in "$fname
	fi

done



