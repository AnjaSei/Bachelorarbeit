#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <map>
#include <stdlib.h>  //atoi
#include <set>


using namespace std;


//splits the string tab by tab
vector<string> split(string line, char delimiter) {

    vector<string> container;
    size_t startpos = 0;
    size_t endpos = line.find(delimiter);

   while (endpos!=string::npos) {
        container.push_back(line.substr(startpos, endpos-startpos));
        startpos = endpos+1;
        endpos= line.find(delimiter, startpos);
        //cout<<"start "<<start<<" end "<<end<<" npos "<<string::npos<<endl;

        if (endpos==string::npos){
            container.push_back(line.substr(startpos));
            //cout<<"ENDE "<<end<<"LALALALALALALALALA"<<line.length()<<endl;
        }
    }

    return(container);

}

//check if the vectors share at least one element
bool contains(vector<string>& anno_gene1, vector<string>& anno_gene2){

    for (unsigned a1=0; a1<anno_gene1.size(); a1++){
            for (unsigned a2=0; a2<anno_gene2.size(); a2++){
                if(anno_gene1[a1]==anno_gene2[a2]){
                    return(true);
                }


            }
        }

    return(false);

}


int main(int argc, char* argv[])
{

	string rcg_file_name;
	string anno_file_name;
	unsigned index=1;

	//parse command line arguments
	while(index<argc) {
		//flag for rcg file
        	if(strcmp(argv[index], "-r")==0){
            		if (index + 1 < argc) {
                		rcg_file_name = argv[index+1];
            		}
            		else { // filename missing after -r
				cerr<<"File with the RCGs is missingafter -r!"<<endl;
                		return(1);
            		}
        	}
		if(strcmp(argv[index], "-a")==0){
		    if (index + 1 < argc) {
		        anno_file_name = argv[index+1];
		    }
		    else { // filename missing after -a
			cerr<<"File with the annoation is missing after -a!"<<endl;
		        return(2);
		    }
		}

		index++;
	}

    	if(rcg_file_name.empty()){
        	cerr<<"File with the RCGs is missing!"<<endl;
        	return(1);
    	}

    	if(anno_file_name.empty()){
        	cerr<<"File with the annoation is missing!"<<endl;
        	return(2);
    	}



	//read file with RCGs
    	ifstream rcg_file;
    	rcg_file.open(rcg_file_name.c_str());

    	if(!rcg_file.is_open()){
        	cerr <<"Coud not open the RCG file!"<< endl;
        	return(1);
    	}
	
	
    	vector<string> container; 				//saves values of each line from the RCG file
    	map<unsigned, vector<unsigned> > rcgs;			//key==RCG, values==GENEIDs within a RCG
    	map<unsigned, vector<unsigned> >::iterator it_rcgs;
	string line;
	
	//read RCG file line per line
	while(getline(rcg_file, line)){
		//line is no comment
        	if(!(line.substr(0, 1) == "#")){

		 	//line is not the headline
		    	if(line.substr(0,3)!="RCG"){
		    
				//split line tab by tab and fill container
				container=split(line, '\t');
				
				//each line should contain: RCG, GeneID, rank and distance value
				if(container.size()!=4){
				    	cerr<<"Line does not have the four entries RCG, GeneID, rank and distance value!"<<endl;

				}
				else{
					//search GeneID which defines a RCG
				    	it_rcgs=rcgs.find(atoi(container[0].c_str()));

					//GeneID was not found -> insert it as a key
				    	if (it_rcgs == rcgs.end()){
				            	rcgs.insert( make_pair(atoi(container[0].c_str()), vector<unsigned> {1,atoi(container[1].c_str())} ));

				    	}
					//GeneID was found -> extend this RCG
				    	else{
				        	rcgs[atoi(container[0].c_str())].push_back(atoi(container[1].c_str()));
				    	}


				}

            		}

        	}

	}
    	rcg_file.close();


	//read the annotation file
    	ifstream annotation_file;
    	annotation_file.open(anno_file_name.c_str());
    	map<unsigned, vector<string> > annotations; 		//save the annotations for each GeneID (==key)
    	map<unsigned, vector<string> >::iterator it_anno;

    	if(!annotation_file.is_open()){
        	cerr <<"Coud not open the annotation file!"<< endl;
        	return(1);
    	}

	//read the annotation file line per line
    	while(getline(annotation_file, line)){

		//line is no comment
        	if(!(line.substr(0, 1) == "#")){
			
			//split line tab by tab and fill the container
               		container=split(line, '\t');

               		if(container.size()<2){
                    	cerr<<"Line does not have the two entries GeneID and annotation!"<<endl;

                	}
                	else{
                    		//search GeneID
                    		it_anno=annotations.find(atoi(container[0].c_str()));
				
				//GeneID was not found -> insert it as a key
                    		if (it_anno == annotations.end()){
                            		annotations.insert( make_pair(atoi(container[0].c_str()), vector<string> {1,container[1].c_str()} ));

                    		}
				//GeneId was found -> add additional annotation	
                   		else{
                        
                        		annotations[atoi(container[0].c_str())].push_back(container[1].c_str());
                    		}


                	}

        	}

    	}
    	annotation_file.close();

	//k=number of most correlated genes to each gene's RCG
	set<unsigned> all_k;
	for (unsigned i=1; i<=5; i++) all_k.insert(i);
	for (unsigned i=10; i<=100; i+=10) all_k.insert(i);

	set<unsigned>::iterator it_k;

	//saves the number of shared annotation pairs for each RCG for each k
	vector<unsigned> shared_anno;

	//for each ranked co-expression group do:
	for(it_rcgs = rcgs.begin(); it_rcgs != rcgs.end(); it_rcgs++) {

  		//load GeneIDs
  		vector<unsigned> rcg_i = rcgs[it_rcgs->first];
  		shared_anno.push_back(rcg_i[0]); //save reference gene

  		//count gene pairs with shared annotation for each k per RCG
  		unsigned counter=0;
  		
		unsigned i;
		//for each gene pair do:
	  	for( unsigned j=1; j<rcg_i.size(); j++){
	    		for(i=0; i<j; i++){

	      			//select GeneIDs
		      		unsigned geneID1=rcg_i[i];
		      		unsigned geneID2=rcg_i[j];

				//select annotations
				it_anno=annotations.find(geneID1);
				vector<string> anno_gene1;
				if(it_anno!=annotations.end()){
			    		anno_gene1=annotations[it_anno->first];
				}

				it_anno=annotations.find(geneID2);
				vector<string> anno_gene2;
				if(it_anno!=annotations.end()){
			    		anno_gene2=annotations[it_anno->first];
				}

				//GeneIDs share at least one annotation
				if(contains(anno_gene1, anno_gene2)){
			    		counter=counter+1;

				}

	      		}

			it_k=all_k.find(j);
	      		//save number of shared gene pairs for a specific k
	      		if( it_k!=all_k.end() ){
				shared_anno.push_back(counter);
	      		}
	    	}


	}



	//save the results
	cout<<"#"<<argv[0]<<endl;
	cout<<"#File with the ranked co-expression groups: "<<rcg_file_name<<endl;
	cout<<"#File with the annotation: "<<anno_file_name<<endl;
	cout<<"#k=number of most correlated genes to each gene's RCG"<<endl;

	unsigned pos=0;
	cout<<"RCG"<<'\t';
	for (it_k=all_k.begin(); it_k!=all_k.end(); it_k++) cout<<"k="<<*it_k<<'\t';
	cout<<""<<endl;
	for(unsigned i=0; i<rcgs.size(); i++){
	    for (unsigned j=0; j<=all_k.size(); j++){
		cout<<shared_anno[pos]<<'\t';
		pos++;
	    }
	    cout<<""<<endl;
	}


    return 0;
}
