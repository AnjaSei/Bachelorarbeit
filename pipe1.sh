#!/bin/bash

#create result file which contains for a specific annotation the number of 
#shared gene pairs per RCG size and distance/correlation measure
#e.g. 
#distance 	k=1 	k=2 	... 	k=100
#euclidean	3	5		700
#pearson	4	7		800

#input files contain the number of gene pairs with shared annotation per gene 
#and RCG size between 1 and 100 for a specific measurement
#inputfolder contains such inputfiles with different measurement but the same annotation
input_folder=$1

#build headline: distance k=1 k=2 ... k=100
printf "%s\t" "distance"
for i in `seq 1 5`
do
	printf "%s\t" k=${i}
done
for i in `seq 10 10 100`
do
	if [ $i != 100 ]
		then
			printf "%s\t" k=${i}
		else
			printf "%s\n" k=${i}
	fi
done


for shared_anno_file in ${input_folder}/*.txt
do
	#echo $shared_anno_file
	fname=$(basename $shared_anno_file)
	dist="unknown"

	#extract the distance/correlation measure
	case $fname in
		*euclidean*) 
			dist="euclidean"
		;;
		*spearman*) 
			dist="spearman"
		;;
		*quadraticpearson*)
			dist="quadraticpearson"
		;;
		*pearson*) 
			dist="pearson"
		;;
		*relintoverlap*) 
			dist="rio"
		;;
		*mutual_information*) 
			dist="mutual_information"
		;;
		
		*cosine*)
			dist="cosine"
		;;
	esac
	
	if [ "$dist" != "unknown" ]
	    then
		
		#colums=awk '{print NF; exit}' $shared_anno_file
		#printf "%s\t%s" ${dist}
		#awk '{OFS="\t"}{for (i=2; i<=NF; i++) a[i-1]+=$i } END {for (i=2; i<=NF; i++) printf a[i-1] OFS; printf "\n"}' $shared_anno_file

		printf "%s\t" ${dist}
		awk '{OFS="\t"}{for (i=2; i<=NF; i++) sum[i-1]+=$i } END {for (i=2; i<=NF; i++) printf "%s%s", sum[i-1], (i==NF?"\n":"\t")}' $shared_anno_file


	    else
		echo "Unknown distance/correlation measure in "$shared_anno_file
	fi



done
