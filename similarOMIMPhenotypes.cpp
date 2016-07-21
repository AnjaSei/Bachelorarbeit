/**********************************************************************
 File: similarOMIMPhenotypes.cpp

 Created by: Anja S.
 Created on: July 21, 2016

 Description: Extract similar OMIM Phenotypes.

 Input: File with the MimMiner similarity matrix
  OMIM Phenotype Similarity Threshold

 Output: File with all OMIM Phenotype pairs with a similarity score
  greater or equal a user specified threshold

************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <stdlib.h>  //atof
#include <string.h> //strcmp

using namespace std;

//split line tab by tab
vector<string> split(string line, char delimiter) {

    vector<string> container;
    size_t startpos = 0;
    size_t endpos = line.find(delimiter);

    while (endpos!=string::npos) {
        container.push_back(line.substr(startpos, endpos-startpos));
        startpos = endpos+1;
        endpos= line.find(delimiter, startpos);


        if (endpos==string::npos){
            container.push_back(line.substr(startpos));

        }
    }

    return(container);

}

//print help message for the user
void printUsageInfo(const char *prgname) {
  cerr  << "Usage: " << prgname << " [options]" << endl
        << "options:" << endl
        << " -m     File with the MimMiner similarity matrix (mandatory)" << endl
        << "        Each line has the form <OMIMID>\\t<similarityValue1>\\t<similarityValue2>..."<< endl
        << " -t     Threshold for phenotype similarity between 0 and 1 (mandatory) " << endl
        << " -h     Print this help message and exit." << endl;

}


int main(int argc, char* argv[])
{

    string input_file_name; //File with the MimMiner similarity matrix
	double threshold=-1;    //Phenotyp Similarity Threshold

	unsigned index=1;

	//parse command line arguments
	while(index<argc) {

		//flag for file with the MimMiner similarity matrix
        if(strcmp(argv[index], "-m")==0){
            if (index + 1 < argc) {
            	input_file_name = argv[index+1];
            }
            else { // filename missing after -m
                cerr<<"File with the MimMiner similarity matrix is missing after -m!"<<endl;
                printUsageInfo(argv[0]);
                return(1);
            }
        }

        //flag for the Phenotyp Similarity Matrix
		if(strcmp(argv[index], "-t")==0){
		    if (index + 1 < argc) {
		        threshold = atof(argv[index+1]);
		    }
		    else { // threshold is missing after -t
                cerr<<"Enter threshold after -t!"<<endl;
                printUsageInfo(argv[0]);
		        return(2);
		    }
		}

		index++;
	}

    if(input_file_name.empty()){
        cerr<<"File with the MimMiner similarity matrix is missing!"<<endl;
        printUsageInfo(argv[0]);
        return(1);
    }

    if(threshold<0 or threshold>1){
        cerr<<"The Phenotype simithreshold is missing or illegal (have to be between 0 and 1)!"<<endl;
        printUsageInfo(argv[0]);
        return(2);
    }


    ifstream input_file;
    input_file.open(input_file_name.c_str());

    if(!input_file.is_open()){
       	cerr <<"Could not open input file!"<< endl;
       	return(1);
   	}

    //save for each OMIM phenotype the similarity to all other phenotypes
	map<string, vector<double> > similarityMatrix;
	map<string, vector<double> >:: iterator similarityMatrix_it;

	vector<string> OMIMIDs;

	vector<string> container;
    string line;

    //for each line of the input file do:
	while(getline(input_file, line)){

		//split line tab by tab
        container=split(line, '\t');

        //save the OMIM ID per line
        OMIMIDs.push_back(container[0]);

        //save the similarity from all OMIM IDs to the current OMIM ID
        vector<double> similarityValues;
		for (unsigned pheno=1; pheno<container.size(); pheno++){

            similarityValues.push_back(atof(container[pheno].c_str()));
        }
        similarityMatrix.insert(make_pair(container[0], similarityValues));

	}



	vector<double> similarityValues;
	string OMIMID="";

	//print header
    cout    <<"#"<<argv[0]<<endl
            <<"#File with the MimMiner similarity matrix: "<<input_file_name<<endl
            <<"#Phenotyp Similarity Threshold: "<<threshold<<endl
            <<"#OMIM Phenotype ID"<<"\t"<<"similar OMIM Phenotyp ID"<<"\t"<<"similarity score"<<endl;

    //for each OMIM Phenotype ID: print OMIM Phenotype IDs with similarity value over the threshold
    for(similarityMatrix_it=similarityMatrix.begin(); similarityMatrix_it!=similarityMatrix.end(); similarityMatrix_it++){

        OMIMID=similarityMatrix_it->first;
        similarityValues=similarityMatrix_it->second;

        for (unsigned value=0; value<similarityValues.size(); value++){

            //the OMIM Phenotypes are similar enough
            if(similarityValues[value]>=threshold){

                cout<<OMIMID<<"\t"<<OMIMIDs[value]<<"\t"<<similarityValues[value]<<endl;
            }
        }
    }



    return 0;
}
