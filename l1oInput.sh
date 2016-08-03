#!/bin/bash
#############################################################################
# File: l1oInput.sh
#
# Created by: Anja S.
# Created on: August 03, 2016
#
# Description: Creates Inputfiles for the Leave-One-Out-Validation
#
# Input: Results fromm OMIM-search "brain"
#  Map OMIMGeneIDs to OMIMPhenoIDs
#  Map OMIMIDs to EntrezGeneIDs
#	For details see below.
#
# Output: Map EntrezGeneID to brain-associated OMIMPhenoIDs
#  Map EntrezGeneID to similarOMIMID and relatedOMIMID
#
############################################################################


thisDir=$(pwd)

######################################### Necessary Data sets and Preprocessing #######################################

#1. Data set -> brain-disease associated OMIMPhenoIDs
#Results from OMIM search "brain":
#Go to http://omim.org/search/advanced/entry
#Search for "brain"
#Search in: Text OR Clinical Synopsis
#Only Records with: Gene Map Locus
#MIM Number Prefix: #phenotype description, molecular basis known
#->Search
#Download As Tab-delimited File
#Attention: This are only the top 200 results!
#Solution: Download step by step 200 results and concatenate this files
#Files are saved in  OMIM_brain_search/

cd OMIM_brain_search/
#Concatenate result files
cat 1-200.txt 201-400.txt 401-600.txt 601-800.txt 801-1000.txt 1001-1200.txt 1201-1400.txt 1401-1504.txt > search-brain-text-clinical-synopsis-known-molecular-basis.txt

#Remove header and description (OMIMPhenoIDs with known molecular basis start with #), select column with the OMIMPhenoID and remove '#' from each OMIMPhenoID
grep "^#" search-brain-text-clinical-synopsis-known-molecular-basis.txt | cut -f1 | sed 's/#//g' > search-brain-text-clinical-synopsis-known-molecular-basis-wo-header-PhenoOMIMIDs.txt

cd $thisDir

#2. Data set -> Map OMIMGeneIDs to OMIMPhenoIDs
#Results from OMIM Gene Map "brain"
#Go to http://omim.org/search/advanced/entry
#Search for "brain"
#Search in: Text OR Clinical Synopsis
#Only Records with: Gene Map Locus
#MIM Number Prefix: #phenotype description, molecular basis known
#->Search
#Retrieve Corresponding Gene Map (returns all genes which are associated with the brain-associated phenotypes plus ALL other with the gene associated phenotypes; it is possible that a phenotype with known molecular basis has also associated genes which were determined by e.g. statistical methods -> remove this entries later)
#Phenotype Only Entries 
#Download As Tab-delimited File
#Attention: This are only the top 200 results!
#Solution: Download step by step 200 results and concatenate this files
#Files are saved in  OMIM_gene_map_brain/

cd OMIM_gene_map_brain/
#Concatenate result files
cat 1-200.txt 201-400.txt 401-600.txt 601-800.txt 801-1000.txt 1001-1200.txt 1201-1400.txt 1401-1466.txt > gene-map-search-brain-text-clinical-synopsis-known-molecular-basis.txt

#Remove header + description for the Phenotype Mapping Key and take only entries with phenomapkey=3, select OMIMGeneID (column 5) and OMIMPhenoID (column 8)
grep -E "^Gene Map Search|^Downloaded|^Cytogenetic location|^Copyright|^$|^Phenotype Mapping Key|^1 - the disorder is placed on the map|^2 - the disorder was placed on the map|^3 - the molecular basis of the disorder|^4 - a contiguous gene duplication" -v gene-map-search-brain-text-clinical-synopsis-known-molecular-basis.txt | awk -F'\t' '{OFS="\t"}{ if( $10 == 3) print $5,$8 }'  > gene-map-search-brain-text-clinical-synopsis-known-molecular-basis-wo-header-phenomapkey-3-OMIMGeneID-OMIMPhenoID.txt 

cd $thisDir

#3. Data set -> Map OMIMIDs to EntrezGeneIDs
#Download mim2gene.txt from https://omim.org/static/omim/data/mim2gene.txt
#Remove the header
grep "^#" -v mim2gene.txt > mim2gene_wo_header.txt


#4. Data set -> Similarity Matrix for OMIMPhenoID pairs
#Download the SoftPanel Similarity Matrix from http://www.isb.pku.edu.cn/SoftPanel/html/download.html

##Take only the first three columns from the original SoftPanel Similarity Matrix (OMIMID, similarOMIDID, similarityScore) and remove quotes around each OMIMID and similarityScore
#cut -f1,2,3 SimilarityMatrix_MTH02_SoftPanel_original | sed 's/\"//g' > SimilarityMatrix_MTH02_SoftPanel

##Take only OMIMID phenotype pairs with similarity score greater or equal 0.5
#awk '{ if( $3 >= 0.5) print $0 }' SimilarityMatrix_MTH02_SoftPanel > SimilarityMatrix_MTH02_SoftPanel_threshold_0.5.tsv




############################# Produce Files for the Leave-One-Out-Validation ###########################


#Join over OMIMPhenoID to get associated OMIMGeneIDs
join -t $'\t' -1 1 -2 2 <(sort -t $'\t' -k1n OMIM_brain_search/search-brain-text-clinical-synopsis-known-molecular-basis-wo-header-PhenoOMIMIDs.txt) <( sort -t $'\t' -k2n OMIM_gene_map_brain/gene-map-search-brain-text-clinical-synopsis-known-molecular-basis-wo-header-phenomapkey-3-OMIMGeneID-OMIMPhenoID.txt) > PhenoOMIMID2GeneOMIMID.txt

#Create first inputfile for Leave-One-Out: EntrezGeneID to OMIMPhenoID
#Join over OMIMGeneID to get associated EntrezGeneIDs; then select OMIMPhenoID (column 2) and EntrezGeneID (column 4)
join -t $'\t' -1 2 -2 1 <(sort -t $'\t' -k2n PhenoOMIMID2GeneOMIMID.txt) <( sort -t $'\t' -k1n mim2gene_wo_header.txt) | awk  -F"\t" '{OFS="\t"}{print $4,$2}'> EntrezGeneID2OMIMPhenoID-OMIM-search-brain-in-clinical-synopsis-and-text_known-molecular-basis-with-gene-locus.txt


#Create second inputfile for Leave-One-Out: EntrezGeneID to similarOMIMID and relatedOMIMID
#Join over OMIMPhenoID to get similar OMIMPhenoIDs, then select EntrezGeneID (column 2), similar OMIMPhenoID (column 3), related OMIMPhenoID (column 1), similarityScore (column 4)
join -t $'\t' -1 2 -2 1 <(sort -t $'\t' -k2n EntrezGeneID2OMIMPhenoID-OMIM-search-brain-in-clinical-synopsis-and-text_known-molecular-basis-with-gene-locus.txt) <( sort -t $'\t' -k1n SimilarityMatrix_MTH02_SoftPanel_threshold_0.5.tsv) | awk  -F"\t" '{OFS="\t"}{print $2,$3,$1,$4}' > EntrezGeneID2simPheno2relPheno_OMIM-search-brain-in-clinical-synopsis-and-text_known-molecular-basis-with-gene-locus-SilarityMatrix_MTH02_SoftPanel_threshold_0.5.tsv



#Remove temporary files
rm OMIM_brain_search/search-brain-text-clinical-synopsis-known-molecular-basis.txt
rm OMIM_gene_map_brain/gene-map-search-brain-text-clinical-synopsis-known-molecular-basis.txt
rm mim2gene_wo_header.txt
rm PhenoOMIMID2GeneOMIMID.txt


