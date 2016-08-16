/**********************************************************************
 File: l10_summary.cpp

 Created by: Anja S.
 Created on: August 10, 2016

 Description: Calculates the Leave-One-Out summary for a user specified
  locus size from a file containing the prioritization results.

 Input: File containing the prioritization results
  Locus size (contains at most 2* Locus_Size + 1 candidates)

 Output: File containing the Leave-One-Out summary (number of occurences of
  the true disease gene among the top k candidates, the expected value and the
  p-value)

************************************************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <map>
#include <math.h>   //pow, ceil
#include <cmath>    //abs
#include <stdlib.h> //atoi
#include <algorithm>    //sort
#include<float.h>   //DBL_MAX
#include <set>
#include <sstream>  //to concatenate strings and ints
#include <stdexcept>

using namespace std;

#define UNSCORED DBL_MAX

//split the string tab by tab
vector<string> split(string line, char delimiter) {

    vector<string> container;
    size_t startpos = 0;
    size_t endpos = line.find(delimiter);

    while (endpos != string::npos) {

        container.push_back(line.substr(startpos, endpos-startpos));
        startpos = endpos + 1;
        endpos = line.find(delimiter, startpos);

        if (endpos == string::npos){

            container.push_back(line.substr(startpos));
        }
    }

    return(container);

}


bool sortByDistance ( vector<double> distance1, vector<double> distance2) {
    return (abs(distance1[3]) < abs(distance2[3]));
}
//smaller score is better than higher score
bool sortByL1OScore ( vector<double> score1, vector<double> score2) {
    return (score1[2] < score2[2]);
}

//call probHyperGeometric to calculate the pValue with the one-tailed Fisher’s exact test
double exec(const char* command) {
    char buffer[128];
    string pValue = "";
    FILE* probHyperGeometric = popen(command, "r");
    if (probHyperGeometric == NULL){
        cerr<<"popen failed!"<<endl;
    }
	while(fgets(buffer, sizeof(buffer), probHyperGeometric) != NULL){
		pValue += buffer;
	}

    pclose(probHyperGeometric);

    return atof(pValue.c_str());
}

//print this help message for the user
void printUsageInfo(const char *prgname) {
  cerr  << "Usage: " << prgname << " [options]" << endl
        << "options:" << endl
        << " -l     File containing the prioritization results" << endl
        << "        (output of mbc-evalCandGenesByCoExpr-leaveOneOut.pl) (mandatory)"<< endl
        << " -s     Locus size (contains at most 2* Locus_Size + 1 candidates)" << endl
        << " -h     Print this help message and exit." << endl;

}

int main(int argc, char* argv[])
{

    //File containing the prioritization results (output of mbc-evalCandGenesByCoExpr-leaveOneOut.pl)
	string l1oFilename = "";
	unsigned locusSize = 0;

	int index = 1;

	//parse command line arguments
	while(index < argc) {

            if(strcmp(argv[index], "-h") == 0){
                    printUsageInfo(argv[0]);
                	return(0);
        	}
        	//flag for the l1o file
	        else if(strcmp(argv[index], "-l") == 0){

        	    if (index + 1 < argc) {
                	l1oFilename = argv[index+1];
                }
            	else {
                		// filename missing after -r
                		cerr<<"File containing the prioritization results (output of mbc-evalCandGenesByCoExpr-leaveOneOut.pl) is missing after -l!"<<endl;
                		printUsageInfo(argv[0]);
                		return(1);
            	}
        	}
        	//flag for locus size
        	else if(strcmp(argv[index], "-s") == 0){
                if (index + 1 < argc) {
                	locusSize = atoi(argv[index+1]);
                }
            	else {
                    cerr<<"The locus size is missing after -s!"<<endl;
                    printUsageInfo(argv[0]);
                    return(2);
                }
        	}

		index++;
	}

	//check if mandatory arguments are given
    if(l1oFilename.empty()){

        cerr<<"File containing the prioritization results (output of mbc-evalCandGenesByCoExpr-leaveOneOut.pl) is missing!"<<endl;
        printUsageInfo(argv[0]);
        return(3);
    }
    if (locusSize <= 0){
        cerr<<"The locus size should be greater than zero!"<<endl;
        printUsageInfo(argv[0]);
        return(4);
    }

	//read file with the  prioritization results
    ifstream l1oFile;
    l1oFile.open(l1oFilename.c_str());

   	if( !l1oFile.is_open() ){

       		cerr <<"Could not open the file containing the prioritization results!"<< endl;
       		return(1);
   	}

    map<pair<int, double>, vector<vector<double> > > OMIMGeneMap; //save the candidates for each OMIM-diseaseGene-Pair
    map<pair<int, double>, vector<vector<double> > >::iterator OMIMGeneMapit;
    bool invalidOMIMGenePair = false; //disease gene has no expression profile
    pair<int,double> thisOMIMGenePair;

	//read RCG file line per line and fill OMIMGeneMap
	string line;
	while( getline(l1oFile, line) ){

            //new OMIM-diseaseGene-pair -> select the OMIMID and the geneID
        	if( line.substr(0, 3) == "+++" ){
                    //example:
                    //+++ OMIM ID=126200 gene=100133941 +++

                    invalidOMIMGenePair = false;

                    //get start position of the OMIM ID
                    string OMIM = "OMIM ID=";
                    size_t found_OMIM_start = line.find(OMIM);
                    if (found_OMIM_start == string::npos){
                        cerr<<"Could not find the OMIM ID in "<<line<<endl;
                        continue;
                    }
                    found_OMIM_start += OMIM.size();

                    //get start position of the gene
                    string gene = " gene=";
                    size_t found_geneID_start = line.find(gene);
                    if (found_geneID_start == string::npos){
                        cerr<<"Could not find the gene in "<<line<<endl;
                        continue;
                    }
                    found_geneID_start += gene.size();

                    //get the position after the gene
                    size_t found_geneID_end = line.find(" +++");
                    if (found_geneID_end == string::npos){
                        cerr<<"Could not find the gene in "<<line<<endl;
                        continue;
                    }

                    //each OMIM entry has an unique six-digit number
                    int OMIMID = atoi(line.substr(found_OMIM_start, 6).c_str());

                    //get the disease gene ID
                    double diseaseGene = atof(line.substr(found_geneID_start, found_geneID_end - found_geneID_start).c_str());

                    thisOMIMGenePair.first = OMIMID;
                    thisOMIMGenePair.second = diseaseGene;


        	}
        	//no expression data available for the diseaseGene
        	//-> skip the candidate of the actual OMIM-diseaseGene-Pair
        	else if(invalidOMIMGenePair){
                continue;
        	}

        	//line is comment or headline
            else if (line.substr(0,3) == "---" || line.substr(0,1) == "#" || line.substr(0,4) == "rank"){
                continue;
            }
            //candidate of the actual OMIM-diseaseGene-Pair
            else{
                //split line tab by tab and fill container
                vector<string> container = split(line, '\t');

                if( container.size() < 9){
                    cerr<<"line should have the form <rank>\\t<gene>\\t<chromosome>\\t<geneStart_bp>\\t<geneEnd_bp>\\t<score_mult>\\t<score_av>\\t<relRank_scored>\\t<relRank_all>"<<endl;
                    cerr<<"and maybe additional columns (extended output)!"<<endl;
                }

                //no expression profile available for the diseaseGene -> skip OMIM-diseaseGene-pair
                else if (atof(container[1].c_str()) == thisOMIMGenePair.second && container[5] == "unscored"){

                    OMIMGeneMap.erase(thisOMIMGenePair);
                    invalidOMIMGenePair = true;
                }

                //save candidate gene for the OMIM-diseaseGene-Pair
                else{
                        double geneID = atof(container[1].c_str());
                        double geneStart = atof(container[3].c_str());
                        double score;

                        //score_mult
                        if( container[5] == "unscored"){
                            score = UNSCORED;
                        }
                        else{
                            score = atof(container[5].c_str());
                        }

                        vector<double> thisCandidate;
                        thisCandidate.push_back(geneID);
                        thisCandidate.push_back(geneStart);
                        thisCandidate.push_back(score);

                        OMIMGeneMapit = OMIMGeneMap.find(thisOMIMGenePair);

                        //first candidate of the actual OMIM-diseaseGene-Pair
                        if ( OMIMGeneMapit == OMIMGeneMap.end() ){

                            vector<vector<double> > candidates;
                            candidates.push_back(thisCandidate);
                            OMIMGeneMap.insert(make_pair(thisOMIMGenePair, candidates));

                        }
                        //further candidate
                        else{
                            (OMIMGeneMapit -> second).push_back(thisCandidate);

                        }

                }

            }

    }

    l1oFile.close();

    //sort candidates for each OMIM-diseaseGene-pair by their close to the disease gene (absolute value)
    //to get the a list of candidate genes for a user specified locus size
    for (OMIMGeneMapit = OMIMGeneMap.begin(); OMIMGeneMapit != OMIMGeneMap.end(); OMIMGeneMapit ++){

        double diseaseGene = (OMIMGeneMapit -> first).second;
        vector<vector<double> > candidateList = OMIMGeneMapit -> second;

        //search the diseaseGene
        double diseaseGeneStart = 0;
        for (unsigned i = 0; i < candidateList.size(); i++){

            double candidateGene = candidateList[i][0];

            if( candidateGene == diseaseGene ){

                diseaseGeneStart = candidateList[i][1];
                break;
            }

        }

        //save distance of all candidateGenes to the diseaseGene
        for (unsigned i = 0; i < candidateList.size(); i++){

            double candidateGeneStart = candidateList[i][1];
            double difference = diseaseGeneStart - candidateGeneStart;
            (OMIMGeneMapit -> second)[i].push_back(difference);

        }

        //sort candidates by their absolute distance to the disease gene
        sort((OMIMGeneMapit -> second).begin(), (OMIMGeneMapit -> second).end(), sortByDistance);

    }

    map<pair<int, double>, vector<vector<double> > > OMIMGeneMapNew;
    map<pair<int, double>, vector<vector<double> > >::iterator OMIMGeneMapNewIt;

    //save the effective candidates (have an expression profile)
    for (OMIMGeneMapit = OMIMGeneMap.begin(); OMIMGeneMapit != OMIMGeneMap.end(); OMIMGeneMapit ++){

        unsigned thisLocusSizeT = 0; //count candidates in direction to the telomere
        unsigned thisLocusSizeC = 0; //count candidates in direction to the centromere

        vector<vector<double> > candidateList = OMIMGeneMapit -> second;
        vector<vector<double> > effectiveCandidateList;

        for (unsigned i = 0; i < candidateList.size(); i++){

            double dist = candidateList[i][3]; //distance on the chromsome between candidate gene and disease gene

            if (thisLocusSizeC < locusSize && dist > 0){

                thisLocusSizeC++;

                //take this candidate -> expression profile available
                if (candidateList[i][2] != UNSCORED){
                    effectiveCandidateList.push_back(candidateList[i]);
                }
            }
            if (thisLocusSizeT < locusSize && dist < 0){

                thisLocusSizeT++;
                //take this candidate -> expression profile available

                if (candidateList[i][2] != UNSCORED){
                    effectiveCandidateList.push_back(candidateList[i]);
                }
            }
            //candidateGene == diseaseGene
            if ((thisLocusSizeC <= locusSize || thisLocusSizeT <= locusSize) && dist == 0){
                //checked before if the disease gene has an expression profile
                effectiveCandidateList.push_back(candidateList[i]);
            }


        }

        OMIMGeneMapNew.insert(make_pair(OMIMGeneMapit->first, effectiveCandidateList));

    }

    //sort the effective candidates by their Leave-One-Out-Validation score
    for (OMIMGeneMapNewIt = OMIMGeneMapNew.begin(); OMIMGeneMapNewIt != OMIMGeneMapNew.end(); OMIMGeneMapNewIt++){

        sort((OMIMGeneMapNewIt -> second).begin(), (OMIMGeneMapNewIt -> second).end(), sortByL1OScore);
    }

    //build header for the result file
    cout<<"N"<<"\t"<<"effectiveCandidates"<<"\t"<<"OMIM-Gene-Pairs"<<"\t"
        <<"RankedFirstObserved"<<"\t"<<"RankedFirstExpected"<<"\t"<<"pValueRankedFirst"<<"\t"
        <<"RankedFirstToThirdObserved"<<"\t"<<"RankedFirstToThirdExpected"<<"\t"<<"pValueRankedFirstToThird"<<"\t"
        <<"RankedFirstToTenthObserved"<<"\t"<<"RankedFirstToTenthExpected"<<"\t"<<"pValueRankedFirstToTenth"<<"\t"
        <<"RankedSmallerTenPercentObserved"<<"\t"<<"RankedSmallerTenPercentExpected"<<"\t"<<"pValueRankedSmallerTenPercent"
        <<endl;

    unsigned N = locusSize;
    unsigned averageNumberEffectiveCandidates = 0;
    unsigned numberOMIMGenePairs = OMIMGeneMapNew.size();
    unsigned rankedFirstObserved = 0;
    unsigned rankedFirstToThirdObserved = 0;
    unsigned rankedFirstToTenthObserved = 0;
    unsigned rankedSmallerTenPercentObserved = 0;

    //count how often the true disease gene has position k in the candidate lists
    for (OMIMGeneMapNewIt = OMIMGeneMapNew.begin(); OMIMGeneMapNewIt != OMIMGeneMapNew.end(); OMIMGeneMapNewIt++){

        double diseaseGene = (OMIMGeneMapNewIt -> first).second;
        vector<vector<double> > candidateList = OMIMGeneMapNewIt -> second;

        //take only OMIM-gene pairs with at least 50 effective candidates
        if (candidateList.size() >= 50){
            averageNumberEffectiveCandidates += candidateList.size();

            for (unsigned i = 0; i < candidateList.size(); i++){

                double candidateGene = candidateList[i][0];

                if (diseaseGene == candidateGene){
                    if (i == 0) rankedFirstObserved++;
                    if (i < 3) rankedFirstToThirdObserved++;
                    if (i < 10) rankedFirstToTenthObserved++;
                    if (i < candidateList.size()*0.1) rankedSmallerTenPercentObserved++;

                }

            }
        }
    }
    averageNumberEffectiveCandidates /= OMIMGeneMapNew.size();

    //expected value = number_gene_phenotype_pairs * ( k/ average_number_candidates)
    double rankedFirstExpected = ceil( numberOMIMGenePairs * (1./averageNumberEffectiveCandidates));
    double rankedFirstToThirdExpected = ceil( numberOMIMGenePairs * (3./averageNumberEffectiveCandidates));
    double rankedFirstToTenthExpected = ceil( numberOMIMGenePairs * (10./averageNumberEffectiveCandidates));
    double rankedSmallerTenPercentExpected = ceil (numberOMIMGenePairs * 0.1); // == numberOMIMGenePairs * ( (0.1*averageNumberEffectiveCandidates) / averageNumberEffectiveCandidates )

    stringstream strstream;
    string command;
    strstream<<"./probHyperGeometric -at_least -a "<<rankedFirstObserved<<" -a_tot "<<numberOMIMGenePairs<<" -n "<<numberOMIMGenePairs<<" -n_tot "<<ceil(numberOMIMGenePairs*averageNumberEffectiveCandidates);
    command = strstream.str();
    double pValueRankedFirst = exec(command.c_str());
    strstream.str(""); strstream.clear();

    strstream<<"./probHyperGeometric -at_least -a "<<rankedFirstToThirdObserved<<" -a_tot "<<numberOMIMGenePairs * 3<<" -n "<<numberOMIMGenePairs <<" -n_tot "<<ceil(numberOMIMGenePairs*averageNumberEffectiveCandidates);
    command = strstream.str();
    double pValueRankedFirstToThird = exec(command.c_str());
    strstream.str(""); strstream.clear();

    strstream<<"./probHyperGeometric -at_least -a "<<rankedFirstToTenthObserved<<" -a_tot "<<numberOMIMGenePairs * 10<<" -n "<<numberOMIMGenePairs <<" -n_tot "<<ceil(numberOMIMGenePairs*averageNumberEffectiveCandidates);
    command = strstream.str();
    double pValueRankedFirstToTenth = exec(command.c_str());
    strstream.str(""); strstream.clear();

    strstream<<"./probHyperGeometric -at_least -a "<<rankedSmallerTenPercentObserved<<" -a_tot "<<ceil(ceil(numberOMIMGenePairs*averageNumberEffectiveCandidates)/10)<<" -n "<<numberOMIMGenePairs <<" -n_tot "<<ceil(numberOMIMGenePairs*averageNumberEffectiveCandidates);
    command = strstream.str();
    double pValueRankedSmallerTenPercent = exec(command.c_str());
    strstream.str(""); strstream.clear();

    cout<<N<<"\t"<<averageNumberEffectiveCandidates<<"\t"<<numberOMIMGenePairs<<"\t"
        <<rankedFirstObserved<<"\t"<<rankedFirstExpected<<"\t"<<pValueRankedFirst<<"\t"
        <<rankedFirstToThirdObserved<<"\t"<<rankedFirstToThirdExpected<<"\t"<<pValueRankedFirstToThird<<"\t"
        <<rankedFirstToTenthObserved<<"\t"<<rankedFirstToTenthExpected<<"\t"<<pValueRankedFirstToTenth<<"\t"
        <<rankedSmallerTenPercentObserved<<"\t"<<rankedSmallerTenPercentExpected<<"\t"<<pValueRankedSmallerTenPercent
        <<endl;

    return 0;
}
