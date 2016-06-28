#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <map>
#include <stdlib.h>  //atoi
#include <set>


using namespace std;



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

    while(index<argc) {
        if(strcmp(argv[index], "-r")==0){
            if (index + 1 < argc) {
                rcg_file_name = argv[index+1];
            }
            else { // filename missing after -i

                return(1);
            }
        }
        if(strcmp(argv[index], "-a")==0){
            if (index + 1 < argc) {
                anno_file_name = argv[index+1];
            }
            else { // filename missing after -i

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




    ifstream rcg_file;
    rcg_file.open(rcg_file_name.c_str());

    string line;
    vector<string> headline;

    if(!rcg_file.is_open()){
        cerr <<"Coud not open the RCG file!"<< endl;
        return(1);
    }

    vector<string> container;
    map<unsigned, vector<unsigned> > rcgs;
    map<unsigned, vector<unsigned> >::iterator it_rcgs;

    while(getline(rcg_file, line)){

        if(!(line.substr(0, 1) == "#")){

         //headline
            if(line.substr(0,3)=="RCG"){
                //cout<<"Found headline"<<endl;
                //split(line, '\t', headline);
            }
            else{
                //fill container
                container=split(line, '\t');
                if(container.size()!=4){
                    cerr<<"Line does not have  four entries!"<<container.size()<<" "<<container[0]<<endl;

                }
                else{
                    //cout<<container[0]<<" "<<container[1]<<" "<<container[2]<<" "<<container[3]<<endl;
                    it_rcgs=rcgs.find(atoi(container[0].c_str()));

                    if (it_rcgs == rcgs.end()){
                            rcgs.insert( make_pair(atoi(container[0].c_str()), vector<unsigned> {1,atoi(container[1].c_str())} ));

                    }
                    else{
                        //cout<<atoi((container[1]).c_str())<<endl;
                        rcgs[atoi(container[0].c_str())].push_back(atoi(container[1].c_str()));
                    }


                }

            }

        }

    }
    rcg_file.close();



    ifstream annotation_file;
    annotation_file.open(anno_file_name.c_str());
    map<unsigned, vector<string> > annotations;
    map<unsigned, vector<string> >::iterator it_anno;

    if(!annotation_file.is_open()){
        cerr <<"Coud not open the annotation file!"<< endl;
        return(1);
    }

    while(getline(annotation_file, line)){

        if(!(line.substr(0, 1) == "#")){

               container=split(line, '\t');
               if(container.size()!=2){
                    cerr<<"Line does not have  two entries!"<<container.size()<<" "<<container[0]<<endl;

                }
                else{
                    //cout<<container[0]<<" "<<container[1]<<" "<<container[2]<<" "<<container[3]<<endl;
                    it_anno=annotations.find(atoi(container[0].c_str()));

                    if (it_anno == annotations.end()){
                            annotations.insert( make_pair(atoi(container[0].c_str()), vector<string> {1,container[1].c_str()} ));

                    }
                    else{
                        //cout<<atoi((container[1]).c_str())<<endl;
                        annotations[atoi(container[0].c_str())].push_back(container[1].c_str());
                    }


                }

        }

    }
    annotation_file.close();

    /*for(it_rcgs = rcgs.begin(); it_rcgs != rcgs.end(); it_rcgs++) {
        //cout <<it_rcgs->first<<endl;
        vector<int> temp=rcgs[it_rcgs->first];
        for(unsigned i=0; i<temp.size(); i++){
            //cout<< temp[i] << endl;
        }
   	}*/

    //for each rcg do
    /*cout<<annotations.size()<<endl;
    for(it_anno = annotations.begin(); it_anno != annotations.end(); it_anno++) {
       cout <<it_anno->first<<" ";
        vector<string> temp=annotations[it_anno->first];
        for(unsigned i=0; i<temp.size(); i++){
           cout<< temp[i] << " ";
        }
       cout<<""<< endl;
   	}*/






//unsigned all_k[] = {1,2,3,4,5,10,20,30,40,50,60,70,80,90,100};
set<unsigned> all_k;
for (unsigned i=1; i<=5; i++) all_k.insert(i);
for (unsigned i=10; i<=100; i+=10) all_k.insert(i);

set<unsigned>::iterator it_k;

vector<unsigned> shared_anno;

//for each ranked co-expression group do:
for(it_rcgs = rcgs.begin(); it_rcgs != rcgs.end(); it_rcgs++) {

  //load GeneIDs
  vector<unsigned> rcg_i = rcgs[it_rcgs->first];
  shared_anno.push_back(rcg_i[0]); //save reference gene

  //count gene pairs with shared annotation for each k per RCG
  unsigned counter=0;
  unsigned i;

  for( unsigned j=1; j<rcg_i.size(); j++){
    for(i=0; i<j; i++){

      //select GeneIDs
      unsigned geneID1=rcg_i[i];
      unsigned geneID2=rcg_i[j];

      //both GeneIDs have an annotation
      //if(!is.na(geneID1) && !is.na(geneID2)){

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




cout<<"#"<<argv[0]<<endl;
cout<<"#File with the ranked co-expression groups: "<<rcg_file_name<<endl;
cout<<"#File with the annotation: "<<anno_file_name<<endl;
cout<<"#k=number of most correlated genes to each gene's RCG"<<endl;

//save results
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


/*for(unsigned i=0; i<shared_anno.size(); i++){
    cout<<shared_anno[i]<<" ";


}*/



    return 0;
}
