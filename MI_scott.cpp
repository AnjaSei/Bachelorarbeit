#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <map>
#include <algorithm> //min_element, max_element
#include <cmath> //pow
#include <utility> //pair
#include <iostream>
#include <stdlib.h>  //atof

using namespace std;




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









int main(int argc, char* argv[])
{


	string input_file_name;
	input_file_name = argv[1];
    

    	ifstream input_file;
    	input_file.open(input_file_name.c_str());

    	string line;
    	vector<string> headline;

    	if(!input_file.is_open()){
        	cerr <<"Could not open input file!"<< endl;
        	return(1);
   	}

	int lines=0;
	vector<string> container;

	//two gene expression profiles
	vector<double> exprProfile1;
	vector<double> exprProfile2;
	while(getline(input_file, line)){

		if(lines==3){
			break;
		}
		if(line.substr(0,3)!="GEN"){
		        //fill container
		        container=split(line, '\t');
			for (unsigned l=1; l<container.size(); l++){
				if(lines==1){
					exprProfile1.push_back(atof(container[l].c_str()));
				}
				if(lines==2){
					exprProfile2.push_back(atof(container[l].c_str()));
				}
			}
            }

	lines++;
	}


	unsigned vecLen=exprProfile1.size();
	
	double min_x=*min_element(exprProfile1.begin(), exprProfile1.end());	

	double min_y=*min_element(exprProfile2.begin(), exprProfile2.end());

	double minimum=min(min_x, min_y);



	//calulate bin's width with Scott's rule
	//bin_size=(3.5*mean(sd(exprProfile1), sd(exprProfile2))) / vecLen^(1/3)
	
	//calculate average mean for the two expression profiles
	double mean_x=0;
	double mean_y=0;
	for (unsigned eV=0; eV<vecLen; eV++){
		mean_x+=exprProfile1[eV];
		mean_y+=exprProfile2[eV];
	}
	mean_x/=vecLen;
	mean_y/=vecLen;

	//calculate average standard deviation for the two expression profiles
	double sd=0;
	double sd_x=0;
	double sd_y=0;
	for (unsigned eV=0; eV<vecLen; eV++){
		sd_x+=pow(exprProfile1[eV]-mean_x, 2);
		sd_y+=pow(exprProfile2[eV]-mean_y, 2);
	}
	sd_x=sqrt((1./vecLen)*sd_x);
	sd_y=sqrt((1./vecLen)*sd_y);
	sd=(sd_x + sd_y)/2.;

	//Scott's Rule
	double bin_size=(3.5*sd)/pow(vecLen,1/3.);

	//formula mutual information: MI=(sum over x element X)(sum over y element Y) p(x,y)*log2(p(x,y)/(p(x)*p(y)))  
	//fill bins	
	double mutInfo=0.0;
	map<int, double> p_bin_x; //bin's probability for exprProfile1
	map<int, double> p_bin_y; //bin's probability for exprProfile2
	map<pair<int, int>, double> p_bin_xy; //bin pairs probability
	map<int, double>::iterator it_x;
	map<int, double>::iterator it_y;
	map<pair<int, int>, double>::iterator it_xy;


	int bin_x=0;
	int bin_y=0;
    	for (unsigned int eV = 0; eV < vecLen; eV++) {
		
		//select bins
		if(exprProfile1[eV]==minimum){
			bin_x=0;
		}
		else{	
			bin_x=ceil((exprProfile1[eV]-minimum)/bin_size)-1;
		}
		if(exprProfile2[eV]==minimum){
			bin_y=0;
		}
		else{	
			bin_y=ceil((exprProfile2[eV]-minimum)/bin_size)-1;
			
		}

		it_x = p_bin_x.find(bin_x);
      		it_y = p_bin_y.find(bin_y);
      		it_xy = p_bin_xy.find(make_pair(bin_x,bin_y));
		

		//check if bin exists
      		//no -> insert bin as a key
      		if (it_x == p_bin_x.end()){
       	 		p_bin_x.insert ( make_pair(bin_x, 1.0) );
        
      		}
      		//yes -> increase occurence of gene expression values in this bin
      		else{
        		p_bin_x[bin_x]++;
      		}

		//check if bin exists
      		//no -> insert bin as a key
      		if (it_y == p_bin_y.end()){
       	 		p_bin_y.insert ( make_pair(bin_y, 1.0) );
        
      		}
      		//yes -> increase occurence of gene expression values in this bin
      		else{
        		p_bin_y[bin_y]++;
      		}

		//check if bin pair exists
      		//no -> insert bin pair as a key
      		if (it_xy == p_bin_xy.end()){
        
        		p_bin_xy.insert( make_pair(make_pair(bin_x, bin_y), 1.0));
        
	      	}
	      	//yes -> increase occurence of gene expression value pairs in this bin
	      	else{
			p_bin_xy[make_pair(bin_x, bin_y)]++;
	      	}
      
	}

    	//determine the probability of each bin from exprProfile1
    	for(it_x = p_bin_x.begin(); it_x != p_bin_x.end(); it_x++) {
		p_bin_x[it_x->first] /= vecLen;
      
    	}
	
    	//determine the probability of each bin from exprProfile2
    	for(it_y = p_bin_y.begin(); it_y != p_bin_y.end(); it_y++) {
      		p_bin_y[it_y->first] /= vecLen;
		 
   	}
	
	
	//determine the probability of each bin pair from exprProfile1 and exprProfile2
    	for(it_xy = p_bin_xy.begin(); it_xy != p_bin_xy.end(); it_xy++) {
      		p_bin_xy[it_xy->first] /= vecLen;
      
    	}

	//calculate the mutual information ->(sum over X)(sum over Y) p(x,y)*log2(p(x,y)/(p(x)*p(y)))
    	for(it_xy = p_bin_xy.begin(); it_xy != p_bin_xy.end(); it_xy++) {
		mutInfo+=(it_xy->second)*log2((it_xy->second)/( p_bin_x[(it_xy->first).first] * p_bin_y[(it_xy->first).second]) );
		
     	}
	cout<<mutInfo<<endl;
	

  
  }
