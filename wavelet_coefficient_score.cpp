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
#include <stdio.h>
#include <string>


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




int main(){



    //transformed data
	string input_file_name;
	//input_file_name = "/home/anja/Dokumente/Studium/6.Semester/Bachelorarbeit/upload/data_sets/haar_wavelet_transformed_data_cube_edge_length_16_two_genes_700_1300.txt";
    //input_file_name= "/home/anja/Dokumente/Studium/6.Semester/Bachelorarbeit/upload/data_sets/haar_wavelet_transformed_data_cube_edge_length_64_10_genes.txt";
    input_file_name="/home/anja/Dokumente/Studium/6.Semester/Bachelorarbeit/upload/data_sets/temp_d_3_c_64/haar_wavelet_transformed_data_cube_edge_length_64_decomposition_level_3.txt";

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


	//save the start and end column for each subband
	map<string, pair<unsigned, unsigned> > subband;
    map<string, pair<unsigned, unsigned> >:: iterator it_subband;

	while(getline(input_file, line)){

		if(lines==3){
			break;
		}
        if( (line.substr(0,6)!="GENEID") && (line.substr(0,8)!="\"GENEID\"") ){

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
        //headline
        else{
            //fill subband map
            headline=split(line, '\t');
            for (unsigned column=1; column<headline.size(); column++){
                //search subband
                it_subband=subband.find(headline[column]);

                //found start column of a subband
                if(it_subband==subband.end()){
                    subband.insert(make_pair(headline[column], make_pair(column-1, column-1)));
                }
                //update end column of a subband
                else{
                    subband[headline[column]].second=column-1;
                }

            }



        }

	lines++;
	}


    //key=decomposition level, value=start and end column of all subbimages with the same decomposion level
    //special case: approximation -> key=0
    map<unsigned, vector<pair<unsigned, unsigned> > > decomposition;
    map<unsigned, vector<pair<unsigned, unsigned> > >::iterator it_decomposition;
    unsigned decompositionLevel;
    string subbandName;

    //save the start and end column for each suband at each decomposition level
    for (it_subband=subband.begin(); it_subband!=subband.end(); it_subband++){

        vector<pair<unsigned, unsigned> > startEndSubband;
        subbandName=it_subband->first;

        //subband is the approximation
        if( subbandName.substr(0,3)=="LLL" ){

            startEndSubband.push_back(it_subband->second);
            decomposition.insert(make_pair(0, startEndSubband));
            continue;
        }

        //erase filter
        decompositionLevel=atoi(subbandName.erase(0,3).c_str());
        it_decomposition=decomposition.find(decompositionLevel);

        //found the first subimage with this decomposition level
        if(it_decomposition==decomposition.end()){

            startEndSubband.push_back(it_subband->second);
            decomposition.insert(make_pair(decompositionLevel, startEndSubband));
        }
        //found another subimage with this decomposition level
        else{

            startEndSubband=decomposition[decompositionLevel];
            startEndSubband.push_back(it_subband->second);
            decomposition[decompositionLevel]=startEndSubband;
        }

    }

    /*for (it_decomposition=decomposition.begin(); it_decomposition!=decomposition.end(); it_decomposition++){
        vector< pair<unsigned, unsigned> > temp=it_decomposition->second;
        unsigned vecLen=temp.size();
        cout<<it_decomposition->first<<endl;
        cout<<"veclen "<<vecLen<<endl;
        for(unsigned v=0; v<vecLen; v++){
            cout<<(it_decomposition->second)[v].first<<'\t'<<(it_decomposition->second)[v].second<<endl;
        }
    }*/


    /*Calcuate the wavelet coefficient score from "A wavelet-coefficient score for comparison of
    two-dimensional climatic-data fields" by Livina et al.

    S0 is obtained from the approximation sub-image coefficients at the max. decomposition level.
    Sj is obtained from the detail sub-image coefficients of decomposition level j.
    E is the wavelet coefficient score and range between -infinity and 1 (E=1 for two identical Wavelet Profiles).


                                    2*( (sum over x element approxSubimage) (x_Profile1-x_Profile2)^2 )
    S0= 1- –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
           (sum over x element approxSubimage) [ (x_Profile1-mean(approxSubimageProfile1))^2 + (x_Profile2-mean(approxSubimageProfile2))^2 ]


                                    2*( (sum over each detailSubimage level j) (sum over x element detailSubimage)(x_Profile1-x_Profile2)^2 )
    Sj= 1- ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
           (sum over each detailSubimage level j) (sum over x element detailSubiage)[ (x_Profile1-mean(detailSubimageProfile1))^2 + (x_Profile2-mean(detailSubimageProfile2))^2 ]


                                                                    Sj + j*S0
    E= (sum over j=0,1,.., max. decomposition level)    ––––––––––––––––––––––––––––––––––––
                                                        (j+1)*(max. decomposition level +1)

    */


    map<unsigned, double> S; //key=decomposition level j, value=Sj with j=0,1,..,max. decomposition level
    map<unsigned, double>::iterator it_S;

    double Sj;
    unsigned startSubimage;
    unsigned endSubimage;
    double meanSubimageProfile1=0.0;
    double meanSubimageProfile2=0.0;
    double numerator=0.0;       //numerator of formular of Sj
    double denominator=0.0;     //denominator of formular of Sj
    vector<pair<unsigned, unsigned> > startEndSubimage; //saves start and end colum for each subimage of a specific decomposition level


    //E: -0.159857

    //calculate Sj
    //for each decomposition level
    for (it_decomposition=decomposition.begin(); it_decomposition!=decomposition.end(); it_decomposition++){

        startEndSubimage=it_decomposition->second;

        //for each subimage
        for (unsigned subimage=0; subimage<startEndSubimage.size(); subimage++){

            startSubimage=startEndSubimage[subimage].first;
            endSubimage=startEndSubimage[subimage].second;

            //for each exprValue in the subimage
            for (unsigned eV=startSubimage; eV<=endSubimage; eV++){

                meanSubimageProfile1+=exprProfile1[eV];
                meanSubimageProfile2+=exprProfile2[eV];
            }
            meanSubimageProfile1/=(endSubimage-startSubimage+1.);
            meanSubimageProfile2/=(endSubimage-startSubimage+1.);

            //for each exprValue in the subimage
            for (unsigned eV=startSubimage; eV<=endSubimage; eV++){

                numerator+=pow(exprProfile1[eV]-exprProfile2[eV], 2);
                denominator+=( pow(exprProfile1[eV]-meanSubimageProfile1, 2) + pow(exprProfile2[eV]-meanSubimageProfile2, 2) );
            }

            //reset
            meanSubimageProfile1=0.0;
            meanSubimageProfile2=0.0;
        }

        numerator*=2;
        Sj=1-numerator/denominator;
        S.insert(make_pair(it_decomposition->first, Sj));

        //reset
        numerator=0.0;
        denominator=0.0;

    }


    //calculate wavelet coefficient score
    double E=0.0;
    double maxDecLevel=S.size()-1;
    for(unsigned j=0; j<=maxDecLevel; j++){
        E+=(S[j]+j*S[0])/((j+1)*(maxDecLevel+1));

    }

    cout<<"E: "<<E<<endl;


}
