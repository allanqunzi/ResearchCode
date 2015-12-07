#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <complex>
#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <functional>
#include "Vector.h"
#include "fftw.cpp"

using namespace std;

constexpr unsigned dim = 3;
template<typename T>
using Vec = std::array<T, dim>;

const string output_folder  = "./dsf_time_check/";
const string output_folder2 = "./dsf_fixed_t/";
const string output_folder3 = "./dsf_fixed_k/";
constexpr unsigned middle = 10, head= 4, tail = 4;
constexpr auto pep  = middle + head + tail;
constexpr unsigned configs = 30000, chain_num = 100, npart = pep * chain_num;
constexpr unsigned save_freq = 100, interval = 10;
constexpr unsigned num_intervals = 2048; // this number should be smaller than configs/interval
constexpr double time_step = 0.005, dt = time_step * save_freq * interval, df = 1.00/(dt * num_intervals);
const double sqrt_npart = sqrt(double(npart));
constexpr Vec<double> L = {30.0, 30.0, 30.0}; // box size
constexpr double ave_par_dist = 1.0;
const Vec<double> dk = {2.0*twoPI/L[0], 2.0*twoPI/L[1], 2.0*twoPI/L[2]};
const Vec<unsigned> min_nk = {0, 0, 0};
constexpr Vec<unsigned> max_nk = {  20, //static_cast<unsigned>( (twoPI/ave_par_dist)/dk[0] ),
                                    20, //static_cast<unsigned>( (twoPI/ave_par_dist)/dk[1] ),
                                    20 //static_cast<unsigned>( (twoPI/ave_par_dist)/dk[2] )
                                  };
constexpr double max_k = max_nk[0] * twoPI/L[1];
constexpr int max_k_magnitude = int(sqrt( max_k * max_k * 3 ));
std::unique_ptr< std::array<std::complex<double>, configs> > k[max_nk[0]][max_nk[1]][max_nk[2]];
std::array<std::complex<double>, num_intervals> skt;
std::array<std::complex<double>, num_intervals> skf;

/*
constexpr std::complex<double> operator "" _i(long double d)
{
    return std::complex<double>{0.0, static_cast<double>(d)};
}
*/
template<unsigned input_size>
void complex_to_real(std::array<std::complex<double>, input_size> & in, double * out)
{
   for (int i = 0; i < input_size; ++i)
    {
        out[2*i] = in[i].real();
        out[2*i+1] = in[i].imag();
    }
}

template<unsigned out_size>
void real_to_complex(double * in, std::array<std::complex<double>, out_size> & out)
{
   for (int i = 0; i < out_size; ++i)
    {
        out[i] = std::complex<double>(in[2*i], in[2*i+1]);
    }
}

template<typename input_value_type, typename output_value_type, unsigned input_size, unsigned output_size>
void correlate(int ivl, int num_ivls, int cfgs, double deltat, //ofstream & of,
    std::array<input_value_type, input_size> & input_data,
    std::array<output_value_type, output_size> & output_data,
    std::function<output_value_type(input_value_type &, input_value_type &)> call_back)
//----------------------------------------------------------------------------
// input_size = cfgs
// output_size = num_ivls
//----------------------------------------------------------------------------
{
    double rl = 0.0, im = 0.0, ab = 0.0, md = 0.0, ab_temp = 0.0;

    for (int i = 0; i < num_ivls; ++i)
    {
        int dd = i * ivl;
        output_data[i] = output_value_type{};

        for (int j = dd; j < cfgs; ++j)
        {
            output_data[i] += call_back(input_data[j], input_data[j-dd]);
        }

        output_data[i] /= output_value_type(cfgs - dd);
/*

        if(i==0){
            rl = output_data[i].real();
            im = output_data[i].imag();
            ab = std::abs(output_data[i]);
            md = ab*ab;
         }

         ab_temp = std::abs(output_data[i])/ab;

         if(i != 0){
            of << setw(12) << (dt*i) <<
            "    "<<setprecision(15)<<output_data[i].real()/rl<<
            //"    "<<setprecision(15)<<output_data[i].imag()/im<<
            "    "<<setprecision(15)<<ab_temp<<
            "    "<<setprecision(15)<<ab_temp*ab_temp<<endl;
        }
*/
    }

}

int main(int argc, char const *argv[])
{
	ifstream infile("coordinates.txt", ios::in);

	std::array< Vec<double>, npart > r;
    int id, type;
    string line, dump;

    // allocate memory for k[][][], each element in k points to
    // an std::array<std::complex<double>, configs>.
    // dim related loop
    for (auto i = min_nk[0]; i < max_nk[0]; ++i)
    {
    	for (auto j = min_nk[1]; j < max_nk[1]; ++j)
    	{
    		for (auto l = min_nk[2]; l < max_nk[2]; ++l)
    		{
    			k[i][j][l].reset(new std::array<std::complex<double>, configs>());
    		}
    	}
    }

    // dump the 1st config
    for(int j=1;j<=(9+npart);j++){
        getline(infile,dump);
    }

    // read the configs, calculate exp{i(kx*x + ky*y + kz*z)}
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
         if(i%1000 == 0 && j == npart)std::cout<<i<<"  "<<id<<"  "<<r[j][0]<<",  "<<r[j][1]<<",  "<<r[j][2]<<std::endl;
        }

        infile.ignore(10, '\n');

        // loop through each k = (kx, ky, kz), dim related loop
        for (auto l = min_nk[0]; l < max_nk[0]; ++l){

	    	for (auto m = min_nk[1]; m < max_nk[1]; ++m)
	    	{
	    		for (auto n = min_nk[2]; n < max_nk[2]; ++n)
	    		{
	    			Vec<double> k_cur = {l*dk[0], m*dk[1], n*dk[2] }; // current k = (kx, ky, kz)
	    			double real = 0.0, imag = 0.0, temp = 0.0;
	    			for (int ii = 0; ii < npart; ++ii)
	    			{
	    				temp = k_cur * r[ii]; // kx*x + ky*y + kz*z
                                       //if(l == 0 && m == 0 && n==0)cout<<temp<<endl;
                                       //cout<<k_cur[0]<<"  "<<k_cur[1]<<"   "<<k_cur[2]<<endl;
                                       //cout<<r[ii][0]<<"  "<<r[ii][1]<<"   "<<r[ii][2]<<endl;
                                       //cout<<temp<<endl;
                                       //cout<<"----------------------"<<endl;
                                       //if(l == 0 && m == 0 && n==0)cout<<temp<<endl;
	    				real += cos(temp);
	    				imag += sin(temp);
	    			}
	    			real /= sqrt_npart;
	    			imag /= sqrt_npart;
	    			auto & arr = *(k[l][m][n]); // arr = std::array<std::complex<double>, configs>
	    			arr[i] = std::complex<double>(real, imag); // save this complex number

                               // if(l == 0 && m == 0 && n==0)cout<<i<<"   "<<arr[i]<<endl;
	    		}
	    	}
	    }


    } // end of loop for each config


    // calculate the correlation for each k = (kx, ky, kz)
    // dim related loop
    double fftw_real[2*num_intervals];
    std::vector<std::vector<std::complex<double> > > skt_k_avged(max_k_magnitude+1,
        std::vector<std::complex<double> >(num_intervals, std::complex<double>{}));
    std::vector<std::vector<std::complex<double> > > skf_k_avged(max_k_magnitude+1,
        std::vector<std::complex<double> >(num_intervals, std::complex<double>{}));
    std::vector<int> counter_k_avged(max_k_magnitude+1, 0);

    for (auto l = min_nk[0]; l < max_nk[0]; ++l){

    	for (auto m = min_nk[1]; m < max_nk[1]; ++m)
    	{
    		for (auto n = min_nk[2]; n < max_nk[2]; ++n)
    		{
    			std::string file =  output_folder+"kx"+std::to_string(l)+"_ky"+std::to_string(m)+"_kz"+std::to_string(n)+"_"
					    			+ std::to_string(dk[0]*l) + "_"
					    			+ std::to_string(dk[1]*m) + "_"
					    			+ std::to_string(dk[2]*n) + ".dat";
				//ofstream outfile(file.c_str(), ios::out);
                using complex = std::complex<double>;
                auto & arr = *(k[l][m][n]);

                correlate<complex, complex, configs, num_intervals>(
                    interval, num_intervals, configs, dt, //outfile,
                    arr, skt,
                    [](complex & c1, complex & c2) ->complex { return c1 * std::conj(c2); }
                    );

                complex_to_real<num_intervals>(skt, fftw_real);

                four1(fftw_real, num_intervals);

                real_to_complex<num_intervals>(fftw_real, skf);

                /*
                double df = 1.00/(dt * num_intervals);
                double skfsq;
                for (int ii = 0; ii < num_intervals; ++ii)
                {
                    skfsq = std::abs(skf[ii]);

                    outfile<<setw(15)<<setprecision(12)<<(df*ii)<<"   "<<setprecision(15)<<skfsq*skfsq<<std::endl;

                }
                */

				//outfile.close();

                Vec<double> k_tmp = { l*dk[0], m*dk[1], n*dk[2] };
                int k_magnitude = int(sqrt(k_tmp * k_tmp));
                assert(skt_k_avged[k_magnitude].size() == skt.size());
                ++counter_k_avged[k_magnitude];

                for (int i = 0; i < skt_k_avged[k_magnitude].size(); ++i)
                {
                    skt_k_avged[k_magnitude][i] += skt[i];
                    skf_k_avged[k_magnitude][i] += skf[i];
                }

    		}
    	}
    }

    double tmp_t, tmp_f;
    for (int i = 0; i <= max_k_magnitude; ++i)
    {
        std::string file =  output_folder3+"sk_t_f_fixed_k_"+std::to_string(i)+".dat";
        ofstream outfile(file.c_str(), ios::out);

        skt_k_avged[i][0] /= double(counter_k_avged[i]);
        skf_k_avged[i][0] /= double(counter_k_avged[i]);

        tmp_t = std::abs(skt_k_avged[i][0]);
        tmp_f = std::abs(skf_k_avged[i][0]);
        tmp_t *= tmp_t;
        tmp_f *= tmp_f;

        for (int j = 1; j < skt_k_avged[0].size(); ++j)
        {
            skt_k_avged[i][j] /= double(counter_k_avged[i]);
            skf_k_avged[i][j] /= double(counter_k_avged[i]);
            outfile<<
            setw(15)<<setprecision(12)<< dt*j << "  " << (std::pow(std::abs(skt_k_avged[i][j]), 2)/tmp_t) <<"  "<<
            setw(15)<<setprecision(12)<< df*j << "  " << (std::pow(std::abs(skf_k_avged[i][j]), 2)/tmp_f) <<"  "<<
            setw(15)<<setprecision(12)<< std::abs(skf_k_avged[i][j])/sqrt(tmp_f) << std::endl;
        }
        outfile.close();
    }

    for(int i = 0; i < skt_k_avged[0].size(); ++i)
    {
      if(i%20 == 0){
        std::string file =  output_folder2+"sk_t_f_fixed_t_or_f_"+std::to_string(i*dt)+"_"+std::to_string(i*df)+".dat";
        ofstream outfile(file.c_str(), ios::out);

        for (int j = 0; j <= max_k_magnitude; ++j)
        {
            outfile<<
            setw(15)<<setprecision(12)<<j<<"  "<< std::abs(skt_k_avged[j][i])<<"   "<<
            setw(15)<<setprecision(12)<<std::abs(skf_k_avged[j][i])<<std::endl;
        }
        outfile.close();
     }
    }

    return 0;
}


