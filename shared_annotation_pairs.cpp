/**********************************************************************
 File: shared_annotation_pairs.cpp

 Created by: Anja S.
 Created on: June 28, 2016

 Description: Count the number of genes per RCG (ranked co-expression group)
   and per RCG size -1 which share at least one annotation.

 Input: File with the RCGs
  File with the annotations

 Output: File with the number of shared annotation pairs per RCG and
  per RCG size -1 (between 1 to 100)

************************************************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <map>
#include <set>


using namespace std;


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

//check if the vectors (with annotations) share at least one element
bool contains(vector<string>& annoGene1, vector<string>& annoGene2){

    for (unsigned a1 = 0; a1 < annoGene1.size(); a1++){
        for (unsigned a2 = 0; a2 < annoGene2.size(); a2++){

            if( annoGene1[a1] == annoGene2[a2] ){
                return(true);
            }
        }
    }

    return(false);

}

//print this help message for the user
void printUsageInfo(const char *prgname) {
  cerr  << "Usage: " << prgname << " [options]" << endl
        << "options:" << endl
        << " -r     File with the ranked co-expression groups (RCGs) (mandatory)" << endl
        << "        Each line has the form <RCG>\\t<GeneID>\\t<rank>\\t<distanceValue>"<< endl
        << "        Each RCG must be sorted by the rank."<< endl
        << " -a     File with the annotation (mandatory) " << endl
        << "        Each line has the form <GeneID>\\t<annotation>"<< endl
        << " -h     Print this help message and exit." << endl;

}

int main(int argc, char* argv[])
{

	string rcgFileName;
	string annoFileName;
	unsigned index = 1;

	//parse command line arguments
	while(index < argc) {
		//flag for the rcg file
	        if(strcmp(argv[index], "-r") == 0){

        	    if (index + 1 < argc) {
                	rcgFileName = argv[index+1];
            		}
            		else {
                		// filename missing after -r
                		cerr<<"File with the RCGs is missing after -r!"<<endl;
                		printUsageInfo(argv[0]);
                		return(1);
            		}
        	}

        	//flag for the annoation file
		if(strcmp(argv[index], "-a") == 0){

            		if (index + 1 < argc) {
		        	annoFileName = argv[index+1];
	    		}
	    		else {
                		// filename missing after -a
                		cerr<<"File with the annoation is missing after -a!"<<endl;
                		printUsageInfo(argv[0]);
	        		return(2);
	    		}
		}

		index++;
	}

	//check if mandatory arguments are given
    	if(rcgFileName.empty()){

       		cerr<<"File with the RCGs is missing!"<<endl;
        	printUsageInfo(argv[0]);
        	return(1);
    	}

    	if(annoFileName.empty()){

       		cerr<<"File with the annotation is missing!"<<endl;
        	printUsageInfo(argv[0]);
       		return(2);
    	}



	//read file with the RCGs
    	ifstream rcgFile;
    	rcgFile.open(rcgFileName.c_str());

   	if( !rcgFile.is_open() ){

       		cerr <<"Could not open the RCG file!"<< endl;
       		return(1);
   	}

   	vector<string> container; 		    //saves values of each line from the RCG file
   	map<string, vector<string> > RCG;	//key==RCG, values==GENEIDs within a RCG
   	map<string, vector<string> >::iterator RCGit;
	string line;

	//read RCG file line per line
	while( getline(rcgFile, line) ){

        	//line is a comment
        	if( line.substr(0, 1) == "#" ){
            		continue;
        	}

        	//line is the headline
        	if( line.substr(0,3) == "RCG" ){
            		continue;
        	}

		//split line tab by tab and fill container
		container = split(line, '\t');

		//each line should contain: RCG, GeneID, rank and distance value
		if( container.size() != 4 ){
	    		cerr<<"Line does not have the four entries RCG, GeneID, rank and distance value!"<<endl;
		}
		else{
			//search GeneID which defines this RCG
	    		RCGit = RCG.find(container[0]);

			//GeneID was not found -> insert it as a key and as the first value of the RCG
	    		if ( RCGit == RCG.end() ){
            			RCG.insert( make_pair(container[0], vector<string> (1,container[1])) );
	    		}
			//GeneID was found -> extend this RCG
			else{
			   	(RCGit->second).push_back(container[1]);
			}

        	}

    	}

	rcgFile.close();


	//read the annotation file
    	ifstream annotationFile;
    	annotationFile.open(annoFileName.c_str());

    	if( !annotationFile.is_open() ){

       		cerr <<"Could not open the annotation file!"<< endl;
       		return(1);
    	}

    	map<string, vector<string> > annotations; 		//save the annotations for each GeneID (==key)
    	map<string, vector<string> >::iterator annoIt;

	//read the annotation file line per line
    	while( getline(annotationFile, line) ){

        	//line is a comment
        	if( line.substr(0, 1) == "#" ){
            		continue;
        	}

		//split line tab by tab and fill the container
        	container = split(line, '\t');

        	if( container.size() < 2 ){
            		cerr<<"Line does not have the two entries GeneID and annotation!"<<endl;
       		}
        	else{
            		//search GeneID
            		annoIt = annotations.find(container[0]);

			//GeneID was not found -> insert it as a key
            		if ( annoIt == annotations.end() ){
                		annotations.insert( make_pair(container[0], vector<string> (1, container[1])) );
            		}
			//GeneID was found -> add additional annotation
            		else{
             			(annoIt->second).push_back(container[1]);
             		}
        	}
    	}
    	annotationFile.close();



	//saves the number of shared annotation pairs for each RCG for each RCG size -1 (==k)
	map<string, vector<unsigned> > sharedAnnoPairs;
	map<string, vector<unsigned> >::iterator sharedAnnoPairsIt;

    	//k=number of most correlated genes to each gene's RCG
	set<unsigned> allK;
	set<unsigned>::iterator allKit;
	for ( unsigned i = 1; i <= 5; i++ ) allK.insert(i);
	for ( unsigned i = 10; i <= 100; i+=10) allK.insert(i);


	//count the number of gene pairs with shared annotation for each RCG for each k
	for( RCGit = RCG.begin(); RCGit != RCG.end(); RCGit++ ) {

  		//load GeneIDs
  		vector<string> thisRCG = RCGit->second;

  		//count shared annotation pairs
  		vector<unsigned> thisSharedAnnoPairs;

  		//count gene pairs with shared annotation
  		unsigned counter=0;

		//for each gene pair do:
	  	for( unsigned j = 1; j < thisRCG.size(); j++ ){
            		for( unsigned i = 0 ; i < j; i++ ){

                		//select GeneIDs
		      		string geneID1=thisRCG[i];
		      		string geneID2=thisRCG[j];

				//select annotations
				annoIt=annotations.find(geneID1);

				vector<string> annoGene1;
				if( annoIt != annotations.end() ){
                    			annoGene1 = annoIt->second;
				}

				annoIt = annotations.find(geneID2);
				vector<string> annoGene2;
				if( annoIt != annotations.end() ){
                    			annoGene2 = annoIt->second;
				}

				//GeneIDs share at least one annotation
				if( contains(annoGene1, annoGene2) ){
                    			counter += 1;
				}

            		}

			allKit = allK.find(j);

	      		//save number of shared gene pairs for a specific k
	      		if( allKit != allK.end() ){
				thisSharedAnnoPairs.push_back(counter);
	      		}
        	}

        	sharedAnnoPairs.insert( make_pair(RCGit->first,thisSharedAnnoPairs) );

	}



	//save the results
	cout<<"#"<<argv[0]<<endl;
	cout<<"#File with the ranked co-expression groups: "<<rcgFileName<<endl;
	cout<<"#File with the annotation: "<<annoFileName<<endl;
	cout<<"#k=number of most correlated genes to each gene's RCG"<<endl;

    	//headline
	cout<<"RCG";
	for ( allKit = allK.begin(); allKit != allK.end(); allKit++ ){
        	cout<<'\t'<<"k="<<*allKit;
	}
	cout<<""<<endl;

	vector<unsigned> thisSharedAnnoPairs;
	for(sharedAnnoPairsIt = sharedAnnoPairs.begin(); sharedAnnoPairsIt != sharedAnnoPairs.end(); sharedAnnoPairsIt++){

        	cout<<sharedAnnoPairsIt->first;
        	thisSharedAnnoPairs = sharedAnnoPairsIt->second;
        	for (unsigned i = 0; i < thisSharedAnnoPairs.size(); i++ ){

                	cout<<'\t'<<thisSharedAnnoPairs[i];
        	}
        	cout<<""<<endl;
	}



    return 0;
}
