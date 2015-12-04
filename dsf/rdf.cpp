#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <functional>
#include "Vector.h"

using namespace std;
constexpr unsigned dim = 3;
template<typename T>
using Vec = std::array<T, dim>;

const string output_folder = "./rdf_data/";
constexpr unsigned middle = 10, head= 4, tail = 4;
constexpr auto pep  = middle + head + tail;
constexpr unsigned configs = 30000, chain_num = 100, npart = pep * chain_num;
constexpr unsigned save_freq = 10, interval = 1;
constexpr double time_step = 0.005, dt = time_step * save_freq * interval;
//constexpr vector<double> L = {30.0, 30.0, 30.0}; // box size
const double L[3] = {30.0, 30.0, 30.0}; // box size
constexpr int rdf_average = 10; // the rdf is averaged over rdf_average configs
constexpr int locals_freq = rdf_average*10;
constexpr int locals_num = configs/locals_freq;
constexpr double dR = 1.0;
constexpr int shells_num = int(sqrt(30.0 * 30.0 * 3)) + 1;


void rdf(std::vector<std::vector<double> > & dist_matrix, double dr, std::vector<double> & rdf) {

	assert(npart == dist_matrix.size() && npart == dist_matrix[0].size());
	std::fill(rdf.begin(), rdf.end(), 0.0);
	auto row_num = dist_matrix.size();
	auto col_num = dist_matrix[0].size();
	for (int i = 0; i < row_num; ++i)
	{
		for(int j = 0; j < col_num; ++j){
			if(i == j)continue;
			assert(int(dist_matrix[i][j]/dr) < shells_num);
			rdf[int(dist_matrix[i][j]/dr)] += 1.0;
		}
	}

	double normalizer = 1.0/(npart * npart * 4.18879 * dr * dr * dr); // 4*Pi/3
	double tmp;
        rdf[0] *= normalizer;
	for (int i = 1; i < rdf.size(); ++i)
	{
		tmp = 1.0/(3*i*i + 3*i +1);
		tmp *= normalizer;
		rdf[i] *= normalizer;
	}
}

void construct_dist_matrix(std::vector<std::vector<double> > & dist_matrix,
	std::array< Vec<double>, npart > & r, int dimension)
{
	double tmp;

	assert(dist_matrix.size() == r.size());

	for (int i = 0; i < dist_matrix.size(); ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			dist_matrix[i][j] = 0.0;
			for(int k = 0; k < dimension; ++k){
				tmp = r[i][k] - r[j][k];
				if(tmp > L[k])tmp -= L[k];
				if(tmp <-L[k])tmp += L[k];
				dist_matrix[i][j] += tmp * tmp;
			}
			dist_matrix[i][j] = sqrt(dist_matrix[i][j]);
			dist_matrix[j][i] = dist_matrix[i][j];

		}

		dist_matrix[i][i] = 0.0;
	}
}

void save_rdf(string & file_name, std::vector<double> & rdf, double dr){

	ofstream outfile(file_name.c_str());
	for (int i = 0; i < rdf.size(); ++i)
	{
		outfile << setw(10) << i*dr << "   "<< setprecision(15) << rdf[i] <<'\n';
	}
}

int main(int argc, char const *argv[])
{
	ifstream infile("coordinates.txt", ios::in);
	std::array< Vec<double>, npart > r;
	std::vector<std::vector<double> > DistMatrix(npart, std::vector<double>(npart, 0.0));
	std::vector<double> single_rdf(shells_num, 0.0);
	std::vector<double> local_rdf (shells_num, 0.0);

	int id, type;
    string line, dump;
    std::array<double, dim>dump_r;

    // dump the 1st config
    for(int j=1;j<=(9+npart);j++){
        getline(infile,dump);
    }

    for(int i=0; i<configs; ++i){

        for(int j=1;j<=9;j++){
            getline(infile,line);
            //cout<<line<<endl;
        }

        // read each config into r[j][l]
        for(int j=1; j<=npart; ++j){
            infile>>id>>type;
            for (int l = 0; l < dim; ++l){
                infile>>r[j][l];
            }
            /*
            for (int l = 0; l < dim; ++l){
                infile>>dump_r[l];
            }
            */
         if(i%1000 == 0 && j == npart)std::cout<<i<<"  "<<id<<"  "<<r[j][0]<<",  "<<r[j][1]<<",  "<<r[j][2]<<std::endl;
        }

        if(i%locals_freq == 0){
        	std::fill(local_rdf.begin(), local_rdf.end(), 0.0);
        }

        if((i%locals_freq) < rdf_average){
        	construct_dist_matrix(DistMatrix, r, dim);
        	rdf(DistMatrix, dR, single_rdf);
        	for (int k = 0; k < single_rdf.size(); ++k)
        	{
        		local_rdf[k] += single_rdf[k];
        	}

        	if((i%locals_freq) == (rdf_average-1)){
        		//cout<<"----------------"<<endl;
	        	for (int k = 0; k < local_rdf.size(); ++k)
	        	{
	        		local_rdf[k] /= rdf_average;
	        	}
        		string fname = output_folder + std::to_string(i - rdf_average +1) + ".dat";
        		save_rdf(fname, local_rdf, dR);
        	}
        }

        infile.ignore(10, '\n');
    }
	return 0;
}
