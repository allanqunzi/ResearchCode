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
const string xyzfilename = "coordinates.txt";
constexpr unsigned middle = 10, head= 4, tail = 4;
constexpr auto pep  = middle + head + tail;
constexpr unsigned configs = 600000, chain_num = 100, npart = pep * chain_num;
constexpr unsigned save_freq = 10, interval = 1;
constexpr double time_step = 0.005, dt = time_step * save_freq * interval;
//constexpr vector<double> L = {30.0, 30.0, 30.0}; // box size
constexpr double L[3] = {30.0, 30.0, 30.0}; // box size
constexpr int rdf_average = 10; // the rdf is averaged over rdf_average configs
constexpr int locals_freq = rdf_average*10;
constexpr int locals_num = configs/locals_freq;
constexpr double dR = 1.0;
const int shells_num = int(sqrt(30.0 * 30.0 * 3)) + 1;
constexpr double rho = double(npart)/(L[0]*L[1]*L[2]);

void rdf(std::vector<std::vector<double> > & dist_matrix,
    std::array< Vec<double>, npart*26 > & images,
    std::array< Vec<double>, npart > real,
    double dr, std::vector<double> & rdf) {

	assert(npart == dist_matrix.size() && npart == dist_matrix[0].size());
	std::fill(rdf.begin(), rdf.end(), 0.0);
	auto row_num = dist_matrix.size();
	auto col_num = dist_matrix[0].size();
	for (int i = 0; i < row_num; ++i)
	{
		for(int j = i+1; j < col_num; ++j){
			//if(i == j)continue;
			assert(int(dist_matrix[i][j]/dr) < shells_num);
            rdf[int(dist_matrix[i][j]/dr)] += 2.0;
		}
	}

    double tmp_r = 0.0;

    for (int i = 0; i < real.size(); ++i)
    {
        for (int j = 0; j < images.size(); ++j)
        {
            tmp_r = 0.0;
            for (int k = 0; k < dim; ++k)
            {
                tmp_r += (real[i][k] - images[j][k])*(real[i][k] - images[j][k]);
            }
            tmp_r = sqrt(tmp_r);
            if(int(tmp_r/dr) < shells_num)rdf[int(tmp_r/dr)] += 1.0;
        }
    }

	double normalizer = 1.0/( (npart) * ((4.0 * 3.1415926/3.0)*(dR*dR*dR) * rho) );
	double tmp;

        rdf[0] *= normalizer;
        double sum = rdf[0];
	for (int i = 1; i < rdf.size(); ++i)
	{
		tmp = 1.0/(3*i*i + 3*i +1);
		tmp *= normalizer;
		rdf[i] *= tmp;
        //sum += rdf[i];
	}
        //cout<<"########## sum = "<<sum<<endl;
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
				//if(tmp > L[k])tmp -= L[k];
				//if(tmp <-L[k])tmp += L[k];
				//tmp = tmp - L[k]*round(tmp/L[k]);
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
	ifstream infile(xyzfilename.c_str(), ios::in);
	std::array< Vec<double>, npart > r;
	std::vector<std::vector<double> > DistMatrix(npart, std::vector<double>(npart, 0.0));
	std::vector<double> single_rdf(shells_num, 0.0);
	std::vector<double> local_rdf (shells_num, 0.0);
    std::array< Vec<double>, npart*26 > image_r;

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

            for (int l = 0; l < dim; ++l){
                infile>>dump_r[l];
            }

         if(i%1000 == 0 && j == npart)std::cout<<i<<"  "<<id<<"  "<<r[j][0]<<",  "<<r[j][1]<<",  "<<r[j][2]<<std::endl;
        }

        if((i%locals_freq) < rdf_average){
            int image_box = 0, cur_start;
            double r_shift[3];
            for (int ia = -1; ia <= 1 ; ++ia)
            {
                for (int ib = -1; ib <= 1; ++ib)
                {
                    for (int ic = -1; ic <= 1 ; ++ic)
                    {
                        if(ia == 0 && ib == 0 && ic == 0)continue;
                        cur_start = image_box * npart;
                        ++image_box;
                        r_shift[0] = ia * L[0];
                        r_shift[1] = ib * L[1];
                        r_shift[2] = ic * L[2];
                        for (int id = 0; id < npart; ++id)
                        {
                            for (int ie = 0; ie < dim; ++ie)
                            {
                                image_r[cur_start + id][ie] = r[id][ie] + r_shift[ie];
                            }

                        }

                    }
                }
            }
        }

        if(i%locals_freq == 0){
        	std::fill(local_rdf.begin(), local_rdf.end(), 0.0);
        }

        if((i%locals_freq) < rdf_average){
        	construct_dist_matrix(DistMatrix, r, dim);
        	rdf(DistMatrix, image_r, r, dR, single_rdf);
        	for (int k = 0; k < single_rdf.size(); ++k)
        	{
        		local_rdf[k] += single_rdf[k];
        	}

        	if((i%locals_freq) == (rdf_average-1)){
        		//cout<<"----------------"<<endl;
        		double sum = 0.0;
	        	for (int k = 0; k < local_rdf.size(); ++k)
	        	{
	        		local_rdf[k] /= rdf_average;
                    sum += local_rdf[k];
	        	}

                        //for (int k = 1; k < local_rdf.size(); ++k){
                        //     local_rdf[k] += local_rdf[k-1];

                        //}
                        //cout<<"-------- sum = "<<sum<<",  local_rdf.back()"<<local_rdf.back()<<endl;

        		string fname = output_folder + std::to_string(i - rdf_average +1) + ".dat";
        		save_rdf(fname, local_rdf, dR);
        	}
        }



        infile.ignore(10, '\n');
    }



	return 0;
}
