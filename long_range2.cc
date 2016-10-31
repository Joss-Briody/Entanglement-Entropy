
//undiagonalmatrix.cc

// Numerical integration: Matrix integrals, diagonalization to follow.

#include <cmath>
#include <iostream>
#include <gsl/gsl_integration.h>
#include "cavlib/constants.hh"
#include <fstream>
#include <iomanip>
#include <vector>
#include <fftw3.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "gsl/gsl_rng.h"


using namespace Eigen;
using namespace std; 

//typedef double FUNC( double , void* );        // define the type of the function needed by GSL:


// This function returns a random number uniformally distributed between 0 and 1, using the
// GSL random number generator.

vector<double> random_uniform( unsigned long int seed = 0, unsigned int size = 1 )//20 )
{
	vector<double> random_v;
   static gsl_rng* rng = 0;

   if ( rng == 0 ) rng = gsl_rng_alloc( gsl_rng_default ); 

	if ( seed != 0 ) gsl_rng_set( rng, seed );

	for( unsigned int i = 0; i < size; i++) 
		random_v.push_back( gsl_rng_uniform( rng ) );

   gsl_rng_free( rng );								// Avoid memory leakage.

   return random_v;
}

void FFT( vector<double>* random_v ) {

  //size_t n = (*random_v).size() ;
	size_t n = 2;
  std::vector<double> inp( (2 * n) -1, 0 );
  std::vector<double> out( (2 * n) -1, 0 );

  // Set the real parts to something quasi-random:
	for(int i = -n + 1; i < n; i++)
	{	
		int k  = i + n - 1;
		inp[k] = i;  
		//inp[i] = (*random_v).size();
	}
  // FFTW wants the addresses of the input and output arrays, but has
  // its own own complex type (actually a typedef). Use a cast:
  fftw_complex* finp   = (fftw_complex*) &inp[0];
  fftw_complex* fout   = (fftw_complex*) &out[0];
  
  fftw_plan plan_forward = 
    fftw_plan_dft_1d (n, finp, fout, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute ( plan_forward );

  std::cout.precision( 3 );

  for ( int i = 0; i < n; i++ )
      std::cout << "(" << inp[2*i] << "," << inp[2*i+1] << ")  (" << 
          out[2*i] << ", " << out[2*i+1] << ")\n";

}

double populate(int idx, int sub_idx, int L, vector<double>* Output_t) {

	double entry = 0;

	if(sub_idx == 4)
		entry = (*Output_t)[2*(L+idx)];        

	else if(sub_idx == 3) 
		entry = -(*Output_t)[2*(L-idx)+1];	    // Output_t looks like: [...,g(-2), f(-1), g(-1), f(0), g(0), f(1), g(1), f(2), g(2), ... ]

	else if(sub_idx == 2) 
		entry = (*Output_t)[2*(L+idx)+1];

	else if(sub_idx == 1)
		entry = -(*Output_t)[2*(L+idx)];

	return entry;
}

struct quench_params { 
   int l;
   double h;
   double h_0;
   double t;
};

//All integrals to be done are in fact real so no need to include imaginary part of integrand

double g_init( double x ) {
	return 0.5*sin(x); //NB factor of i!!!
}

double g( double x ) {
	return 0.5*sin(x); //NB factor of i!!!
}

double m_init( double x) {
	return 300;
}

double m( double x ) {
	return 0.9999;
}

double omega_0_init( double x ) {
	//double alpha = 3; //0.05;
	//double y = ( pow(abs(x),alpha) );// + pow(m_init(x), alpha), 0.5) );
	return -cos(x);
}

double omega_0( double x ) {
	//double alpha = 3;//0.05;
	//double y = ( pow(abs(x),alpha) );//+ pow(m(x), alpha), 0.5) );
	return -cos(x);
	//return -cos(x);
}

double omega_init( double x ) {
	return pow( pow( omega_0_init(x) + m_init(x), 2 ) + pow( g_init(x), 2 ), 0.5);
}

double omega( double x ) {
		return pow( pow( omega_0(x) + m(x), 2 ) + pow( g(x), 2 ), 0.5);
}

double SinPhi( double x ) {
	return (1/(omega(x)*omega_init(x))) * ( g(x)*(omega_0_init(x) + m_init(x)) - g_init(x)*(omega_0(x) + m(x)) );
}

double CosPhi( double x ) {
	return (1/(omega(x)*omega_init(x))) * ( ( omega_0_init(x) + m(x) ) * ( omega_0(x) + m_init(x) ) + g(x)*g_init(x) );
}

double f_func(double y, double l, double t, const int N) {

	double x = 2*(C::pi)*y/N;
	return sin(l*x)*( SinPhi(x)*sin(2*omega(x)*t) );
}

double g_func(double y, double l, double t, const int N) {
	
	double x = 2*(C::pi)*y/N;
	double cos_theta = ( omega_0(x) + m(x) ) / omega(x);
	double sin_theta = -g(x)/omega(x);

	return  ( CosPhi(x)*(cos_theta*cos(l*x) - sin_theta*sin(l*x) ) - SinPhi(x)*cos(2*omega(x)*t)*( cos_theta*sin(l*x) +cos(l*x)*sin_theta ) );
}

double integrate( double (*func)(double, double, double, const int), double u, const int N, int l, double t, double h_rand ) {
	double sum = 0;
 	//struct quench_params Quench = {l, h_rand, 300, t};

	
	for( int x = -(N-1)/2; x < (N+1)/2; x ++  ) {
		double f = (*func)(x,l,t,N)/N;
		sum += f;
	}
	return sum;
}

//This function takes the data for the system at a particular time and creates the correlation matrix at that time, to be diagonalized.
void Correlation_matrix( vector< vector<double> >& correlation_mat_t, int t_band, int L, vector< vector<double> >* out )
{	
		vector<double> correlation_mat_t_a;
  	   vector<double> correlation_mat_t_b;
		correlation_mat_t_a.resize(0);
		correlation_mat_t_b.resize(0);
		correlation_mat_t.resize(0);	

		int sub_idx = 1;
		int j = 0;
		int row = 0;
		vector< vector<double> > Output = *(out);

		for(int i = 0; i < 4*L*L; i++) 											//sub_idx values correspond to different entries in 2x2 PI matrices
		{
			if( sub_idx == 1 || sub_idx == 2 )			
				correlation_mat_t_a.push_back( populate(j, sub_idx, L, &Output[t_band] ) );	

			if( sub_idx == 3 || sub_idx == 4 )			
				correlation_mat_t_b.push_back( populate(j, sub_idx, L, &Output[t_band] ) );				
																				
			sub_idx ++;																	//fill in one PI_l, then increment l

			if( sub_idx == 5) {
				sub_idx = 1;
				j--;
			}
				if(j == row-L) {
					correlation_mat_t.push_back(correlation_mat_t_a);     //Build up correlation matrix (at a given time) two rows at a time
					correlation_mat_t.push_back(correlation_mat_t_b);
					correlation_mat_t_a.resize(0);							   //Reset row vectors
					correlation_mat_t_b.resize(0);
					row++;
					j = row;		
				}	
		}	
}

double Entropy(vector<double>* eigenvalues_t, int t_band) {
	vector<double> eigenvals = *(eigenvalues_t);
	double entropy = 0;
	double epsilon = 1e-6;   //changing this gives some nans

	for(unsigned int i = 0; i < eigenvals.size(); i++ )
   {
		if(abs(eigenvals[i]-1) > epsilon && abs(eigenvals[i]+1) > epsilon && eigenvals[i] > -1)  {		//(1-x)*ln(1-x) not defined for x = 1. However, this is a rounding error 
			entropy += (-0.5*(1+eigenvals[i])*log(0.5*(1+eigenvals[i])) );											//and the limit x->1 is clearly 0 (physically this can be mapped to an uncoupled fermionic mode)
		}	
	}
	return entropy;
}

int main (void)
{
	vector<double> random_v = random_uniform(); // try a quench where the B-field is quenched to a random number.
															  // set number of random realizations taken inside random_uniform() 

	for(unsigned int i = 0; i < random_v.size(); i++ )
		cout << random_v[i] << endl;
		
	const int L     = 40;        //dimension of lattice - computing time limited by number of integrals, proportional to L*L. 
	const double T  = 50.0;		  //time to evolve system over.
	const double dt = 0.05;		  //size of time-step.
	const double N  = 1000;
   double err = 0.01;
   double u   = C::pi;
	
	FFT( &random_v );

	ofstream ENTROPY;
	ENTROPY.open("entropy_longrange_test2.txt");		

	vector< vector<double> > Entropy_v_rand;

	for(unsigned int rand_it = 0; rand_it < random_v.size(); rand_it++) 
	{
		double h_rand = 0.999; //random_v[rand_it];

		vector< vector < double > > Output;
	
		for( double t = 0.0; t < T; t += dt )	
		{
			vector<double> Output_t;

			for( int l = -L; l <= L; l++ ) 
			{
   			double f_U = integrate( &f_func, u, N, l, t, h_rand );
   			double g_U = integrate( &g_func, u, N, l, t, h_rand );
				Output_t.push_back(f_U);
				Output_t.push_back(g_U);
		   } 

			Output.push_back(Output_t);
   	}

// Have all the data stored in Output, one vector for each time-step. Now need to correctly index this data:

		vector<double> Entropy_v(0); //initialise to zero for each random number
	
		vector< vector<double> > correlation_mat_t;

		for( int t_band = 0; t_band < static_cast<int>(T/dt); t_band ++ ) 
		{
		
			Correlation_matrix( correlation_mat_t, t_band, L, &Output ); // produces the correlation matrix
															
			MatrixXd mat = MatrixXd::Constant(2*L, 2*L, 0);
	
			for (int i = 0; i < 2*L; i++) {
  				mat.row(i) = VectorXd::Map(&correlation_mat_t[i][0],correlation_mat_t[i].size());
			}

			EigenSolver<MatrixXd> es(mat, false);

			vector<double> eigenvalues_t;  //_t
	
			for(unsigned int i = 0; i < 2*L; i++ )
				eigenvalues_t.push_back(es.eigenvalues().imag()[i]);

			Entropy_v.push_back( Entropy(&eigenvalues_t, t_band) );
		}
	
		Entropy_v_rand.push_back( Entropy_v );
	}

	vector<double> Entropy_average;
	for(unsigned int t_band = 0; t_band < static_cast<int>(T/dt); t_band ++ )
	{
		double entropy_av_t = 0;	
		for( unsigned int rand_it = 0; rand_it < random_v.size(); rand_it++ )
		{
			entropy_av_t += Entropy_v_rand[rand_it][t_band]/random_v.size();
		}
	
		Entropy_average.push_back( entropy_av_t );
	}

//Print out the results.
	int count = 0;
	for( unsigned int t_band = 0; t_band < static_cast<int>(T/dt); t_band ++ ) 
	{
		if(t_band) {
			if(abs(Entropy_average[t_band]-Entropy_average[t_band-1]) < 0.1)     // do not plot numerical errors 
				 ENTROPY << t_band*dt << "\t" << Entropy_average[t_band] << endl;
		
			else if(abs(Entropy_average[t_band]-Entropy_average[t_band+1]) < 0.1) {
				ENTROPY << t_band*dt << "\t" << Entropy_average[t_band] << endl;
			}
			else
				count++;
		}
	//ENTROPY << t_band*dt << "\t" << Entropy_average[t_band] << endl;
	}	
	cout<<"count = "<<count<<endl;
	ENTROPY.close();

}

//print out correlaion matrix...
 /*int count = 0;
	for( unsigned int m = 0; m < correlation_mat.size(); m++) {
		for( unsigned int r = 0; r < correlation_mat[m].size(); r++) {
			count++;
			cout << m << "\t" << r << "\t" << correlation_mat[m][r] << endl;
		}
	}
	cout << "num of entries:" << "\t" << count << endl; 
 */
