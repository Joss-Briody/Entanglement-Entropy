//undiagonalmatrix.cc

// Numerical integration: Matrix integrals, diagonalization to follow.

#include <cmath>
#include <iostream>
#include <gsl/gsl_integration.h>
#include "cavlib/constants.hh"
#include <fstream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std; 


typedef double FUNC( double , void* );        // define the type of the function needed by GSL:

double populate(int idx, int sub_idx, int L, vector<double>* Output_t) {

	double entry = 0;

	if(sub_idx == 4)
		entry = (*Output_t)[2*(L+idx)];         ///can be made easily more efficient - actually probably not here but in tabulating integrals

	else if(sub_idx == 3) 
		entry = -(*Output_t)[2*(L-idx)+1];	    // Output_t looks like: [...,g(-2), f(-1), g(-1), f(0), g(0), f(1), g(1), f(2), g(2), ... ]

	else if(sub_idx == 2) 
		entry = (*Output_t)[2*(L+idx)+1];

	else if(sub_idx == 1)
		entry = -(*Output_t)[2*(L+idx)];

	return entry;
}

//All integrals to be done are in fact real so no need to include imaginary part of integrand

struct quench_params { 
   int l;
   double h;
   double h_0;
   double t;
};

double Epsilon( double x, double h ) {
	return pow((h-cos(x))*(h-cos(x)) + sin(x)*sin(x), 0.5);
}

double f_func(double x, void* p) {
   struct quench_params * params = (struct quench_params *)p;
    
   double l         = (params->l);
   double h         = (params->h);
   double h_0       = (params->h_0);
   double t         = (params->t);

   double E   = Epsilon(x,h); 
   double E_0 = Epsilon(x,h_0);

   return (1/(2*C::pi))*sin(l*x)*sin(x)*(h_0-h)*sin(2*E*t)/(E*E_0);    // odd in l

}

double g_func(double x, void* p) {
   struct quench_params * params = (struct quench_params *)p;

   double l         = (params->l);
   double h         = (params->h);
   double h_0       = (params->h_0);
   double t         = (params->t);
	
   double E   = Epsilon(x,h); 
   double E_0 = Epsilon(x,h_0);

   double C_x = ( 1 - cos(x)*(h+h_0) + h*h_0 )/(E*E_0); 
   double S_x = (sin(x)*(h_0-h)/(E*E_0));

   return (1/(2*C::pi))*(1/E)*( (cos(l*x)*( (cos(x)-h)*C_x - sin(x)*S_x*cos(2*E*t) ) - sin(l*x)*( (cos(x)-h)*S_x*cos(2*E*t) + sin(x)*C_x ) ) ); //one term odd in l, one term even in l - could halve #integrals performed.

}

double integrate( FUNC func, void* params, double x0, 
                  double x1, double& error, int l, double t )
{
   const int max_intervals = 1000; ////maybe change?
   gsl_integration_workspace* workspace = 
   gsl_integration_workspace_alloc ( max_intervals );
   double result;
   gsl_function F;
   struct quench_params Quench = {l, 0.99999, 100, t};
   F.function = func;
   F.params = &Quench;
 
   // Use a GSL routine for integrals:

   gsl_integration_qags( &F, x0, x1, 1e-4, 1e-7, max_intervals,
                        workspace, &result, &error); 
   gsl_integration_workspace_free (workspace);

   return result;
}

double Entropy(vector<double>* eigenvalues_t) {
	vector<double> eigenvals = *(eigenvalues_t);
	double entropy = 0;
	//double entropy_t = 0.0;
	double epsilon = 1e-4;   //changing this gives some nans
	for(unsigned int i = 0; i < eigenvals.size(); i++ ) {
		if(abs(eigenvals[i]-1) > epsilon && abs(eigenvals[i]+1) > epsilon)  {						//(1-x)*ln(1-x) not defined for x = 1. However, this is a rounding error and the limit x->1 is clearly 0 ( physically this is
			entropy += (-0.5*(1+eigenvals[i])*log(0.5*(1+eigenvals[i])) );	
		}
		else
			entropy += (0.0);
	}
		
	return entropy;
}

int main (void)
{
	const int L = 20;
	const double T = 40.0;
	const double dt = 0.05;

   double err = 0.1;
   double u   = C::pi;
   
	vector< vector < double > > Output;
	
	for( double t = 0.0; t < T; t += dt )	
	{
		vector<double> Output_t;
		for( int l = -L; l <= L; l++ ) 
		{
   		double f_U = integrate( f_func, &u, -u, u, err, l, t );
   		double g_U = integrate( g_func, &u, -u, u, err, l, t );
			Output_t.push_back(f_U);
			Output_t.push_back(g_U);
	   } 

		Output.push_back(Output_t);
   }

// Have all the data stored in Output, one vector for each time-step. Now need to correctly index this data:

 //vector< vector<double> > correlation_mat;


	int sub_idx = 1;
	int j = 0;
	int row = 0;
	int t_band = static_cast<int>(25.00/dt);
	cout<< "hello  " << Output.size();
	double entropy_t;
	//for( double t_band = 0; t_band < 15; t_band ++ ) 
	//{
	   vector<double> correlation_mat_t_a;
		vector<double> correlation_mat_t_b;
		vector< vector<double> > correlation_mat_t;

			correlation_mat_t_a.resize(0);
			correlation_mat_t_b.resize(0);
			correlation_mat_t.resize(0);	

		for(int i = 0; i < 4*L*L; i++) 
		{
			if( sub_idx == 1 || sub_idx == 2 )			
				correlation_mat_t_a.push_back( populate(j, sub_idx, L, &Output[t_band] ) );	

			if( sub_idx == 3 || sub_idx == 4 )			
				correlation_mat_t_b.push_back( populate(j, sub_idx, L, &Output[t_band] ) );				
																				
			sub_idx ++;

			if( sub_idx == 5) {
				sub_idx = 1;
				j--;
			}
				if(j == row-L) {
					//cout << row << "\t" << correlation_mat_t_a.size() << "\t" << correlation_mat_t_b.size()<< "\t" << sub_idx << endl;
								
					correlation_mat_t.push_back(correlation_mat_t_a);     //Build up correlation matrix (at a given time) two rows at a time
					correlation_mat_t.push_back(correlation_mat_t_b);
					correlation_mat_t_a.resize(0);							   //Reset row vectors
					correlation_mat_t_b.resize(0);
					//cout << row << "\t" << correlation_mat_t_a.size()<<endl;
					row++;
					j = row;		
				}	
		}	
																	
//print out some stuff...
 /*int count = 0;
	for( unsigned int m = 0; m < correlation_mat.size(); m++) {
		for( unsigned int r = 0; r < correlation_mat[m].size(); r++) {
			count++;
			cout << m << "\t" << r << "\t" << correlation_mat[m][r] << endl;
		}
	}
	cout << "num of entries:" << "\t" << count << endl; 
 */

//int SIZE = correlation_mat_t.size();



//MatrixXd mat(SIZE, SIZE);
MatrixXd mat = MatrixXd::Constant(2*L, 2*L, 0);


for (int i = 0; i < 2*L; i++) {
  mat.row(i) = VectorXd::Map(&correlation_mat_t[i][0],correlation_mat_t[i].size());
}
//cout<< mat << endl;
EigenSolver<MatrixXd> es(mat, false);

vector<double> eigenvalues_t;  //_t
for(unsigned int i = 0; i < 2*L; i++ )
	eigenvalues_t.push_back(es.eigenvalues().imag()[i]);

	
entropy_t = Entropy(&eigenvalues_t);

cout<<"entropy is: " << entropy_t<<endl;
//}

}

/////
//undiagonalmatrix.cc

// Numerical integration: Matrix integrals, diagonalization to follow.

#include <cmath>
#include <iostream>
#include <gsl/gsl_integration.h>
#include "cavlib/constants.hh"
#include <fstream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std; 


typedef double FUNC( double , void* );        // define the type of the function needed by GSL:

double populate(int idx, int sub_idx, int L, vector<double>* Output_t) {

	double entry = 0;

	if(sub_idx == 4)
		entry = (*Output_t)[2*(L+idx)];         ///can be made easily more efficient - actually probably not here but in tabulating integrals

	else if(sub_idx == 3) 
		entry = -(*Output_t)[2*(L-idx)+1];	    // Output_t looks like: [...,g(-2), f(-1), g(-1), f(0), g(0), f(1), g(1), f(2), g(2), ... ]

	else if(sub_idx == 2) 
		entry = (*Output_t)[2*(L+idx)+1];

	else if(sub_idx == 1)
		entry = -(*Output_t)[2*(L+idx)];

	return entry;
}

//All integrals to be done are in fact real so no need to include imaginary part of integrand

struct quench_params { 
   int l;
   double h;
   double h_0;
   double t;
};

double Epsilon( double x, double h ) {
	return pow((h-cos(x))*(h-cos(x)) + sin(x)*sin(x), 0.5);
}

double f_func(double x, void* p) {
   struct quench_params * params = (struct quench_params *)p;
    
   double l         = (params->l);
   double h         = (params->h);
   double h_0       = (params->h_0);
   double t         = (params->t);

   double E   = Epsilon(x,h); 
   double E_0 = Epsilon(x,h_0);

   return (1/(2*C::pi))*sin(l*x)*sin(x)*(h_0-h)*sin(2*E*t)/(E*E_0);    // odd in l

}

double g_func(double x, void* p) {
   struct quench_params * params = (struct quench_params *)p;

   double l         = (params->l);
   double h         = (params->h);
   double h_0       = (params->h_0);
   double t         = (params->t);
	
   double E   = Epsilon(x,h); 
   double E_0 = Epsilon(x,h_0);

   double C_x = ( 1 - cos(x)*(h+h_0) + h*h_0 )/(E*E_0); 
   double S_x = (sin(x)*(h_0-h)/(E*E_0));

   return (1/(2*C::pi))*(1/E)*( (cos(l*x)*( (cos(x)-h)*C_x - sin(x)*S_x*cos(2*E*t) ) - sin(l*x)*( (cos(x)-h)*S_x*cos(2*E*t) + sin(x)*C_x ) ) ); //one term odd in l, one term even in l - could halve #integrals performed.

}

double integrate( FUNC func, void* params, double x0, 
                  double x1, double& error, int l, double t )
{
   const int max_intervals = 1000; ////maybe change?
   gsl_integration_workspace* workspace = 
   gsl_integration_workspace_alloc ( max_intervals );
   double result;
   gsl_function F;
   struct quench_params Quench = {l, 0.99999, 100, t};
   F.function = func;
   F.params = &Quench;
 
   // Use a GSL routine for integrals:

   gsl_integration_qags( &F, x0, x1, 1e-4, 1e-7, max_intervals,
                        workspace, &result, &error); 
   gsl_integration_workspace_free (workspace);

   return result;
}

double Entropy(vector<double>* eigenvalues_t) {
	vector<double> eigenvals = *(eigenvalues_t);
	double entropy = 0;
	//double entropy_t = 0.0;
	double epsilon = 1e-4;   //changing this gives some nans
	for(unsigned int i = 0; i < eigenvals.size(); i++ ) {
		if(abs(eigenvals[i]-1) > epsilon && abs(eigenvals[i]+1) > epsilon)  {						//(1-x)*ln(1-x) not defined for x = 1. However, this is a rounding error and the limit x->1 is clearly 0 ( physically this is
			entropy += (-0.5*(1+eigenvals[i])*log(0.5*(1+eigenvals[i])) );	
		}
		else
			entropy += (0.0);
	}
		
	return entropy;
}

int main (void)
{
	const int L = 20;
	const double T = 20.0;
	const double dt = 0.05;

   double err = 0.1;
   double u   = C::pi;
   
	vector< vector < double > > Output;
	
	for( double t = 0.0; t < T; t += dt )	
	{
		vector<double> Output_t;
		for( int l = -L; l <= L; l++ ) 
		{
   		double f_U = integrate( f_func, &u, -u, u, err, l, t );
   		double g_U = integrate( g_func, &u, -u, u, err, l, t );
			Output_t.push_back(f_U);
			Output_t.push_back(g_U);
	   } 

		Output.push_back(Output_t);
   }

// Have all the data stored in Output, one vector for each time-step. Now need to correctly index this data:

 //vector< vector<double> > correlation_mat;


	int sub_idx = 1;
	int j = 0;
	int row = 0;
	int t_band = static_cast<int>(25.00/dt);
	cout<< "hello  " << Output.size();
	double entropy_t;
	//for( double t_band = 0; t_band < 15; t_band ++ ) 
	//{
	   vector<double> correlation_mat_t_a;
		vector<double> correlation_mat_t_b;
		vector< vector<double> > correlation_mat_t;

			correlation_mat_t_a.resize(0);
			correlation_mat_t_b.resize(0);
			correlation_mat_t.resize(0);	

		for(int i = 0; i < 4*L*L; i++) 
		{
			if( sub_idx == 1 || sub_idx == 2 )			
				correlation_mat_t_a.push_back( populate(j, sub_idx, L, &Output[t_band] ) );	

			if( sub_idx == 3 || sub_idx == 4 )			
				correlation_mat_t_b.push_back( populate(j, sub_idx, L, &Output[t_band] ) );				
																				
			sub_idx ++;

			if( sub_idx == 5) {
				sub_idx = 1;
				j--;
			}
				if(j == row-L) {
					//cout << row << "\t" << correlation_mat_t_a.size() << "\t" << correlation_mat_t_b.size()<< "\t" << sub_idx << endl;
								
					correlation_mat_t.push_back(correlation_mat_t_a);     //Build up correlation matrix (at a given time) two rows at a time
					correlation_mat_t.push_back(correlation_mat_t_b);
					correlation_mat_t_a.resize(0);							   //Reset row vectors
					correlation_mat_t_b.resize(0);
					//cout << row << "\t" << correlation_mat_t_a.size()<<endl;
					row++;
					j = row;		
				}	
		}	
//}
														
//print out some stuff...
 /*int count = 0;
	for( unsigned int m = 0; m < correlation_mat.size(); m++) {
		for( unsigned int r = 0; r < correlation_mat[m].size(); r++) {
			count++;
			cout << m << "\t" << r << "\t" << correlation_mat[m][r] << endl;
		}
	}
	cout << "num of entries:" << "\t" << count << endl; 
 */

//int SIZE = correlation_mat_t.size();


//}
//MatrixXd mat(SIZE, SIZE);
MatrixXd mat = MatrixXd::Constant(2*L, 2*L, 0);

for (int i = 0; i < 2*L; i++) {
  mat.row(i) = VectorXd::Map(&correlation_mat_t[i][0],correlation_mat_t[i].size());
}
//cout<< mat << endl;
EigenSolver<MatrixXd> es(mat, false);

vector<double> eigenvalues_t;  //_t
for(unsigned int i = 0; i < 2*L; i++ )
	eigenvalues_t.push_back(es.eigenvalues().imag()[i]);

	
entropy_t = Entropy(&eigenvalues_t);

cout<<"entropy is: " << entropy_t<<endl;
//}

}


