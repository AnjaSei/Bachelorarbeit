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
#include <string> //erase

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

    //0.0877507

    //transformed data
	string input_file_name;
	input_file_name = "/home/anja/Dokumente/Studium/6.Semester/Bachelorarbeit/upload/data_sets/haar_wavelet_transformed_data_cube_edge_length_16_two_genes_700_1300.txt";
    //input_file_name= "/home/anja/Dokumente/Studium/6.Semester/Bachelorarbeit/upload/data_sets/haar_wavelet_transformed_data_cube_edge_length_64_10_genes.txt";

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


    unsigned subband_size=0;

    vector<double> feature_exprProfile1;
    vector<double> feature_exprProfile2;

    //energy of a subband
    double energy_exprProfile1=0.0;
    double energy_exprProfile2=0.0;

    //standard deviation of a subband
    double sd_exprProfile1=0.0;
    double sd_exprProfile2=0.0;

    //mean of a subband
    double mean_exprProfile1=0.0;
    double mean_exprProfile2=0.0;

    //calculate Wavelet signature for each subband (energy, mean, standard deviation)
    for (it_subband=subband.begin(); it_subband!=subband.end(); it_subband++){

        subband_size=( (it_subband->second).second - (it_subband->second).first +1 );

        //for each gene expression value of a subband
        for (unsigned value=(it_subband->second).first; value<=(it_subband->second).second; value++){

            energy_exprProfile1+=pow(exprProfile1[value], 2);
            energy_exprProfile2+=pow(exprProfile2[value], 2);

            mean_exprProfile1+=exprProfile1[value];
            mean_exprProfile2+=exprProfile2[value];

        }

        energy_exprProfile1/=subband_size;
        energy_exprProfile2/=subband_size;

        mean_exprProfile1/=subband_size;
        mean_exprProfile2/=subband_size;

        //calculate standard deviation
        for (unsigned value=(it_subband->second).first; value<=(it_subband->second).second; value++){
                sd_exprProfile1+=pow(exprProfile1[value]-mean_exprProfile1, 2);
                sd_exprProfile2+=pow(exprProfile2[value]-mean_exprProfile2, 2);
        }

        sd_exprProfile1=sqrt((1./subband_size)*sd_exprProfile1);
        sd_exprProfile2=sqrt((1./subband_size)*sd_exprProfile2);


        feature_exprProfile1.push_back(energy_exprProfile1);
        feature_exprProfile1.push_back(mean_exprProfile1);
        feature_exprProfile1.push_back(sd_exprProfile1);

        feature_exprProfile2.push_back(energy_exprProfile2);
        feature_exprProfile2.push_back(mean_exprProfile2);
        feature_exprProfile2.push_back(sd_exprProfile2);

        //cout<<energy_exprProfile1<<'\t'<<mean_exprProfile1<<'\t'<<sd_exprProfile1<<endl;
        //cout<<energy_exprProfile2<<'\t'<<mean_exprProfile2<<'\t'<<sd_exprProfile2<<endl;

        //reset
        energy_exprProfile1=0.0;
        energy_exprProfile2=0.0;

        mean_exprProfile1=0.0;
        mean_exprProfile2=0.0;

        sd_exprProfile1=0.0;
        sd_exprProfile2=0.0;




    }


    //compare Wavelet signatures of the two expression profiles with the Euclidean distance
    double distance=0.0;


    //for each feature
    for (unsigned feature=0; feature<feature_exprProfile1.size(); feature++){

        distance+= pow(feature_exprProfile1[feature]- feature_exprProfile2[feature], 2);
    }

    distance=sqrt(distance);
    cout<<distance<<endl;

}
