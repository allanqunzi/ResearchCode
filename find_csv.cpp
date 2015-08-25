# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <utility>
# include <map>
# include <algorithm> 
# include <cctype>
# include <cassert>
#include <iterator>
# include <ctime>

using namespace std;

const bool strict_flag = true; // true: turn on CG pairing;   false: turn off CG pairing
const string target("GTCT"), target1("TCTG");

//////////   function declarations   ///////////////////
string & uppercase(string & str);
void reverse_base(string & str, std::map<char, char> & basepair);
template<bool>
bool find(const std::pair<string, string> & p, std::map<string, std::vector<int> > & res, 
	std::map<char, char> & basepair, ofstream & outfile, const bool & CG_flag);
template<bool>
bool find_tar(const std::pair<string, string> & p, std::map<string, std::vector<int> > & res, 
    std::map<char, char> & basepair, ofstream & outfile, const bool & CG_flag);
////////// end of function of declarations   ///////////

int main(int argc, char * argv[])
{
  if(argc <= 1 || argc > 2){
    cout<<"Usage:"<<endl;
    cout<<argv[0]<<" filename"<<endl;
    cout<<"filename is a *.txt, and there shouldn't be spaces in *"<<endl;
    return 1;
  }

    string file(argv[1]);

    std::map<char, char> basepair;
    basepair['A'] = 'T';
    basepair['T'] = 'A';
    basepair['C'] = 'G';
    basepair['G'] = 'C';

    string line;
    string seqname;
    string temp;
    string temp_complement;

    ifstream infile(file.c_str(), std::ios::in|std::ios::ate);
    if(infile){
          try{
           std::ifstream::streampos filesize = infile.tellg();
           line.reserve(filesize);
           temp.reserve(filesize);
           temp_complement.reserve(filesize);

           infile.seekg(0);}catch(...){ assert("The file is too big!! " && 0); };

    }else { cout<<"bad file!!"<<endl; assert(0); }

    cout<<endl;
    cout<<file<<"  is being analyzed, ";
    if(strict_flag)cout<<"and CG pairing is turned on."<<endl;
    else cout<<"and CG pairing is turned off."<<endl;
    cout<<endl;

    string file1_name = "found_" + file.substr(0, file.rfind('.')) + ".dat";
    string file2_name = "result_" + file.substr(0, file.rfind('.')) + ".csv";
    string file3_name = "complementary_found_" + file.substr(0, file.rfind('.')) + ".dat";


    ofstream outfile1(file1_name.c_str(), ios ::out);
    ofstream outfile2(file2_name.c_str(), ios ::out);
    ofstream outfile3(file3_name.c_str(), ios ::out);
    //ofstream outfile4("complementary_NormalReverse_order.dat", ios ::out);

    std::pair<string, string>  data;

    std::map<string, std::vector<int> > res;
    std::map<string, std::vector<int> > res_tar;
    std::map<string, std::vector<int> >::iterator it;
    bool flag1, flag2;
    bool cflag1, cflag2;

    std::pair<string, string>  cdata;
    std::map<string, std::vector<int> > cres;
    std::map<string, std::vector<int> > cres_tar;

    time_t start = time(NULL); 
    while(getline(infile, line)){
    	if(line[0] == '>'){ 
    		//cout<<line<<endl;
    		if(!seqname.empty()){ 

          /////////////////////////////////////////////////////////
          //       analyze complementary sequence first          //
          /////////////////////////////////////////////////////////

          temp_complement.resize( temp.size());
          std::reverse_copy(temp.begin(), temp.end(), temp_complement.begin());

          for(auto & value : temp_complement){ 
            if(value == 'A' || value == 'T' || value == 'C' || value == 'G')
            value = basepair[value];
          }

          cdata = std::make_pair(seqname, std::move(temp_complement));
          cflag1 = find<true>(cdata, cres, basepair, outfile3, strict_flag);
          cflag2 = find_tar<true>(cdata, cres_tar, basepair, outfile3, strict_flag);
          if(cflag1 || cflag2)
          outfile3<<endl;
          //outfile3<<endl;
          temp_complement.clear();
          /////////////////////////////////////////////////////////
          //       end of complementary sequence                 //
          /////////////////////////////////////////////////////////

                  //data.push_back(std::make_pair(seqname, temp)); 
                  data = std::make_pair(seqname, std::move(temp));
                  flag1 = find<false>(data, res, basepair, outfile1, strict_flag);
                  flag2 = find_tar<false>(data, res_tar, basepair, outfile1, strict_flag);
                  if(flag1 || flag2)
                  outfile1<<endl;
                  //outfile1<<endl;
                  /*
                  if( flag1 || flag2  ){ 

    		      outfile1<<data.first<<endl;  
    		      //outfile1<<data.second<<endl; 
    	             }
                   */

                  temp.clear(); 
                }
    		seqname = line; 
    	}
    	else temp.append(uppercase(line));
    }




    temp_complement.resize( temp.size());
    //std::reverse_copy(std::begin(temp), std::end(temp), std::begin(temp_complement));
    std::reverse_copy(temp.begin(), temp.end(), temp_complement.begin());

    for(auto & value : temp_complement){ 
        if(value == 'A' || value == 'T' || value == 'C' || value == 'G')
            value = basepair[value];
          }

    cdata = std::make_pair(seqname, std::move(temp_complement));
    cflag1 = find<true>(cdata, cres, basepair, outfile3, strict_flag);
    cflag2 = find_tar<true>(cdata, cres_tar, basepair, outfile3, strict_flag);
    //outfile3<<endl;
    //outfile3<<endl;

    //data.push_back(std::make_pair(seqname, temp));
    data = std::make_pair(seqname, std::move(temp));
    flag1 = find<false>(data, res, basepair, outfile1, strict_flag);
    flag2 = find_tar<false>(data, res_tar, basepair, outfile1, strict_flag);
    /*
    if( flag1 || flag2  ){ 

       outfile1<<data.first<<endl;  
       outfile1<<data.second<<endl; 
    }
    */

    int counter1 = 0, counter2 = 0; 
////////// first merge the reverse complementary  ///////////////
    for(it = cres_tar.begin(); it != cres_tar.end(); ++it)
    {
      if(cres[it->first].size() == 0)cres[it->first] = it->second;
      else cres[it->first].insert( cres[it->first].end(), it->second.begin(), it->second.end());
      
    }
    for(it = cres.begin(); it != cres.end(); ++it)
    {
      std::sort((it->second).begin(), (it->second).end());
    }
/////////////////////////////////////////////////////////////////
////////////////// then merge the orginal  //////////////////////
    for(it = res_tar.begin(); it != res_tar.end(); ++it)
    {
      if(res[it->first].size() == 0)res[it->first] = it->second;
      else res[it->first].insert( res[it->first].end(), it->second.begin(), it->second.end());
      
    }
    for(it = res.begin(); it != res.end(); ++it)
    {
      std::sort((it->second).begin(), (it->second).end());

      outfile2<<it->first;
      for(int i = 0; i < it->second.size(); ++i)outfile2<<","<<(it->second)[i];
      outfile2<<endl;

      outfile2<<"RevComplement_"<<it->first;
      for(int i = 0; i < cres[it->first].size(); ++i)outfile2<<","<<cres[it->first][i];
      outfile2<<endl;
      cres.erase(it->first);
    }

    if(!cres.empty())outfile2<<" Something is wrong, please report to Aiqun Huang"<<endl;

    cout<<"The following files have been generated:"<<endl;
    cout<<file1_name<<endl;
    cout<<file2_name<<endl;
    cout<<file3_name<<endl;
    cout<<endl;

    time_t end = time(NULL);
    double seconds = difftime(end, start);
    cout<<"time elapsed: "<< int(seconds/60.00) <<" minutes and "<<(long int)(seconds)%60<< " seconds"<<endl;
    cout<<endl;

    return 0;
}


string & uppercase(string & str){
  for (int i = 0; i < str.size(); ++i)
  {
    str[i] = toupper(str[i]);
  }
  return str;
}
void reverse_base(string & str, std::map<char, char> & basepair)
{
  char a;
  for (int i = 0; i < str.size(); ++i)
  {
    a = str[i];
    if(a == 'A' || a == 'T' || a == 'C' || a == 'G')
    str[i] = basepair[a];
  }

}   
template<bool complementary_flag>
bool find(const std::pair<string, string> & p, std::map<string, std::vector<int> > & res, 
  std::map<char, char> & basepair, ofstream & outfile, const bool & CG_flag)
{
  int size = p.second.size(), bound = size - 6;
  string key = p.first;
  string pre, post, pre_pair0, pre_pair1, post_pair0, post_pair1;

  std::vector<int> temp;

    std::size_t found = p.second.find(target); 

    int sep_min = 4;   

    while(found != std::string::npos){
      if(found >= (2+sep_min+3+3) && found+3+3 < size  )
      {
        if(CG_flag)pre = p.second.substr(found-2, 3);
        else pre = p.second.substr(found-2, 2);
        post = p.second.substr(found+4, 3);
        reverse_base(pre, basepair);
        reverse_base(post,basepair);
        std::reverse(pre.begin(), pre.end());
        std::reverse(post.begin(), post.end());

            for(int j = sep_min; j <= 7; ++j){
                if(found-2-j-3-3 >= 0){
                       if(CG_flag)pre_pair0 = p.second.substr(found-(2+j+3), 3);
                       else pre_pair0 = p.second.substr(found-(2+j+2), 2);
                       post_pair0 = p.second.substr(found-(2+j+3+3), 3);

                       pre_pair1 =  pre_pair0;             
                       if( found-2-j-3-3-1 >= 0)post_pair1 = p.second.substr(found-(2+j+3+3+1), 3);
                       else post_pair1 = "";

                       if( (pre == pre_pair0 && post == post_pair0) || 
                           (pre == pre_pair1 && post == post_pair1 ) )
                         {
                            if(complementary_flag)temp.push_back(size-found);
                            else temp.push_back(found);
                            j += 5; 
                            outfile<<" normal order  "<<p.first<<", found at "<<found<<
                            ",  "<<(found-20)<<"-"<<(found+19)<<":  "<<p.second.substr(found-20, 40)<<endl;
                            

                         }
                }
            }

      }
      found = p.second.find(target, found+1);
    }
    res.insert( std::make_pair(key, temp) );
    if(!temp.empty()){ return true; }
    else return false;
}

template<bool complementary_flag>
bool find_tar(const std::pair<string, string> & p, std::map<string, std::vector<int> > & res, 
    std::map<char, char> & basepair, ofstream & outfile, const bool & CG_flag)
{
    int size = p.second.size(), bound = size - 18;
    string key = p.first;
    string pre, post, pre_pair0, pre_pair1, post_pair0, post_pair1;

    std::vector<int> temp;

    std::size_t found = p.second.find(target1); 

    int sep_min = 4;
    while(found != std::string::npos){
        if(found >= 3 && (found+3+2+sep_min+3+3) < size  )
        {
            pre = p.second.substr(found-3, 3);
            if(CG_flag)post = p.second.substr(found+3, 3);
            else post = p.second.substr(found+3+1, 2);
            reverse_base(pre, basepair);
            reverse_base(post,basepair);
            std::reverse(pre.begin(), pre.end());
            std::reverse(post.begin(), post.end());

            for(int j = sep_min; j <= 7; ++j){

                if(found+3+2+j+3+3 < size){
                    pre_pair0 = p.second.substr(found+3+2+j+3+1, 3);
                    if(CG_flag)post_pair0 = p.second.substr(found+3+2+j+1, 3);
                    else post_pair0 = p.second.substr(found+3+2+j+1, 2);

                    post_pair1 = post_pair0;

                    if(found+3+2+j+3+3+1 < size){ pre_pair1 = p.second.substr(found+3+2+j+3+1+1, 3); }
                    else pre_pair1 = "";

                    if(  (pre == pre_pair0 && post == post_pair0) ||
                         (pre == pre_pair1 && post == post_pair1) )
                    {
                        if(complementary_flag)temp.push_back(size-found);
                        else temp.push_back(found);
                        j += 5;
                        
                        outfile<<" reverse order  "<<p.first<<", found at "<<found<<
                        ",  "<<(found-20)<<"-"<<(found+19)<<":  "<<p.second.substr(found-20, 40)<<endl;

                    }

                }
            }
        }

        found = p.second.find(target1, found+1);
    }
    res.insert( std::make_pair(key, temp) );
    if(!temp.empty()){    return true; }
    else return false;
}

