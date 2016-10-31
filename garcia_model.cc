#include <cmath>
#include <iostream>
#include <gsl/gsl_integration.h>
#include "cavlib/constants.hh"   //don't need thesse libraries, see if works without them and a different makefile.
#include <fstream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "gsl/gsl_rng.h"

//Entropy for quasiperiodic fibonacci Ising spin chain, using the general method given in "Entanglement entropy with localized and extended interface defects" appendix.

//IMPORTANT: Eigen stores eigenvectors as nth column in eigenvector matrix, with eigen-value nth value in eigenvalue vector.
//For both matrices of eigenvectors psi and phi, the corresponding vectors of eigenvalues are stored in the same order so can be accessed via the same element.
//eigenvectors are automatically normalized.

using namespace Eigen;
using namespace std;

vector<double> random_uniform( unsigned long int seed, unsigned int size )
{
	vector<double> random_v;
   static gsl_rng* rng = 0;
	double H = 20;

   if ( rng == 0 ) rng = gsl_rng_alloc( gsl_rng_default ); 

	if ( seed != 0 ) gsl_rng_set( rng, seed );

	for( unsigned int i = 0; i < size; i++) 
		random_v.push_back( H*gsl_rng_uniform( rng ) );

   gsl_rng_free( rng );								// Avoid memory leakage.

   return random_v;
}

double F_func( int r, const double alpha, const double lambda, const double N ) 
{
	double f = pow( 1 + pow(lambda,-2)*pow( (sin(C::pi*r/N))/(C::pi/N), 2*alpha), -0.5);

	return f;
}


vector<double> J_initiate(int Length, const double J, double r, const long double omega) {
	//long double omega = 0.5*(sqrt(5L) + 1);
	vector<double> J_i(Length);

	for(int i = 1; i <= Length; i++) {
		J_i[i-1] = J*pow( r, 1 + static_cast<int>(floor(i/omega)) - static_cast<int>(floor((i+1)/omega)) );
		//cout<<1 + static_cast<int>(floor(i/omega)) - static_cast<int>(floor((i+1)/omega) )<<endl;
	}

	return J_i;
}

double A_pop( int i, int j, double h, vector<double>* J, const int N_SIGN) {
	int L = (*J).size();	
	if(i==j)
		return -h;

	//else if(i-j == 1 )//&& i != L-1 )&& j !=L-1)
	//	return -0.5*(*J)[i];

	//else if(i == L-1 && j == 0) 
	//	return 0.5*(*J)[L-1]*N_SIGN;

	//else
	//	return (-0.5*(*J)[i]/( pow( (abs(i-j)), 0.8 ) ) ); //0;

	else
		return ( F_func(i-j, 0.8, 2, 987 ) );
}
double B_pop(int i, int j, double h, vector<double>* J, const int N_SIGN) {
	/*int L = (*J).size();
	if(i-j == 1 )//&& i!= L-1) && j!= L-1)
		return (*J)[i]/2;
	
	//else if(i == L-1 && j == 0)
	//	return 0.5*(*J)[L-1]*N_SIGN;
	else if( i==j)
		return 0;	

	else
		return (0.5*(*J)[i]/( pow(abs(i-j), 0.8 ) ) ); //0; */
	return 0;
}

//G(0)_(i,j)
void G_0(MatrixXd* psi_mat_p, MatrixXd* phi_mat_p, MatrixXd& static_correlation, const int L) {
	
	double z;
	for( int i = 0; i < L; i++ ) {
		for(int j = 0; j < L; j++ ) {
			z = 0;
			for(int q = 0; q < L; q++) {
   			z += ( (*psi_mat_p).col(q)(i) )*( (*phi_mat_p).col(q)(j) );
			}	
			static_correlation(i,j) = -z;
		}
	}
}
void Contractions_init( MatrixXd& A_A, MatrixXd& B_B, MatrixXd& A_B, 
				   MatrixXd* phi_mat_p, MatrixXd* psi_mat_p, 
					VectorXd* eigenval_v_p, 
					const double t, const int L)  { 											//need to check that using correct eigenvalues/vectors
	
	double x,y,z;
	for( int i = 0; i < L; i++ ) {
		for(int j = 0; j < L; j++ ) {
			x = y = z = 0;
			for(int q = 0; q < L; q++) {
				if( (*eigenval_v_p)(q) < 0 && abs((*eigenval_v_p)(q)) < 0.0001 )  // have some eigenvalues ~ -1e-15 i.e = 0.
					(*eigenval_v_p)(q) = 0;

				x += sin(sqrt((*eigenval_v_p)(q))*t)*( (*phi_mat_p).col(q)(i) )*( (*psi_mat_p).col(q)(j) );  //eigenvalues all positive (as spectrum bounded)
				y += cos(sqrt((*eigenval_v_p)(q))*t)*( (*psi_mat_p).col(q)(i) )*( (*psi_mat_p).col(q)(j) );
				z += cos(sqrt((*eigenval_v_p)(q))*t)*( (*phi_mat_p).col(q)(i) )*( (*phi_mat_p).col(q)(j) );
	
			}	
			A_A(i,j) = z;
			B_B(i,j) = y;
			A_B(i,j) = x;
		}
	}
}

void Correlation_mat_create( MatrixXd* A_A, MatrixXd* B_B, MatrixXd* A_B, MatrixXd* static_correlation, MatrixXd& correlation_mat, const int L ) {

	MatrixXd prod_1, prod_2, prod_3, prod_4;
	
	prod_1 = (*A_B)*(*static_correlation)*(*A_A);  //in matrix notation, only 4 different matrices used to build correlation matrix in (eqA16 of paper)
	prod_2 = (*B_B)*(*static_correlation)*(*A_A);
	prod_3 = (*A_B)*(*static_correlation)*(*A_B);
	prod_4 = (*B_B)*(*static_correlation)*(*A_B);
	
	for(int l=1; l <= L; l++) {
		for(int m = 1; m <= L; m++) {
			correlation_mat(2*l-2,2*m-2) =  ( prod_1(l-1,m-1) - prod_1(m-1,l-1) )*pow(1, (l+m-2)); //relative signs different to equations in paper as have cancelled 
			correlation_mat(2*l-2,2*m-1) = -( prod_2(m-1,l-1) + prod_3(l-1,m-1) )*pow(1, (l+m-2)); //i's because matrix is ultimately real so don't need to introduce complex matrices.
			correlation_mat(2*l-1,2*m-2) = ( prod_3(m-1,l-1) + prod_2(l-1,m-1) )*pow(1, (l+m-2));
			correlation_mat(2*l-1,2*m-1) =  ( prod_4(m-1,l-1) - prod_4(l-1,m-1) )*pow(1, (l+m-2));
		}
	}
}

double Entropy(VectorXd* eigenvalues ) {

	double entropy = 0;
	//double epsilon = 1e-6;   //changing this gives some nans
	double x;
	for(unsigned int i = 0; i < (*eigenvalues).size(); i++ )
   {
		x = (-0.5*(1+(*eigenvalues)(i))*log(0.5*(1+(*eigenvalues)(i))) );
		//if(abs((*eigenvalues)(i)-1) > epsilon && abs((*eigenvalues)(i)+1) > epsilon && (*eigenvalues)(i) > -1)  {		//(1-x)*ln(1-x) not defined for x = 1. However, this is a rounding error 
		//	entropy += (-0.5*(1+(*eigenvalues)(i))*log(0.5*(1+(*eigenvalues)(i))) );											//and the limit x->1 is clearly 0 (physically this can be mapped to an uncoupled fermionic mode)
		if( x <10000)
			entropy += x;		
	//}	
	}
	return entropy;
}

int main() {

	const int L      = 987;//377;//377; //144; //377;  // L = F_14
	const int sub_L  = 34;//377;//233;	
	const double h   = 0.5; //1;//0.50;
	const double h_0 = 1000000000;//0;//1000000000000;
	const double r   = 1;//0.25;//0.75; 						 // r = 1 -> homogeneous chain.  
	//const double t   = 600.0;
	const long double omega = 0.5L*(sqrt(5L) + 1L);
	const double J   = pow(r,-(1-(1/omega)));
	const int N_SIGN = 0;//-1;//pow((-1),L); 				//exp(i*pi*N) factor in Jordan-Wigner.

	MatrixXd A(L,L), B(L,L), A_0(L,L), B_0(L,L);

	vector<double> J_i = J_initiate(L, J, r, omega);

	for(int i = 0; i < L; i++) {               //Hamiltonian after the quench
		for(int j = 0; j < L; j++){
			A(i,j) = A_pop(i, j, h, &J_i, N_SIGN);
			A(j,i) = A(i,j);						    //A symmetric (cf. normal modes)

			B(i,j) = B_pop(i, j, h, &J_i, N_SIGN);
			B(j,i) = -B(i,j); 					    //B is antisymmetric
			}
	}

	for(int i = 0; i < L; i++) {               //Hamiltonian before the quench.
		for(int j = 0; j < L; j++){
			A_0(i,j) = A_pop(i, j, h_0, &J_i, N_SIGN);
			A_0(j,i) = A_0(i,j); 					 //A_0 symmetric

			B_0(i,j) = B_pop(i, j, h_0, &J_i, N_SIGN);
			B_0(j,i) = -B(i,j);  					 //B is antisymmetric
			}
	}

	MatrixXd D(L,L), E(L,L), D_0(L,L), E_0(L,L);

	D = (A + B)*(A - B);
	E = (A - B)*(A + B);
	
	D_0 = (A_0 + B_0)*(A_0 - B_0);
	E_0 = (A_0 - B_0)*(A_0 + B_0);

	/*cout<<"commute= "<< endl;
	cout<<A*B-B*A<<endl;
	cout<<"A= "<< endl;
	cout<<A<<endl;
	cout<<"B= "<< endl;
	cout<<B<<endl;
	cout<<"A_0= "<< endl;
	cout<<A_0<<endl;
	cout<<"B_0= "<< endl;
	cout<<B_0<<endl;
	cout<<"D= "<< endl;
	cout<<D<<endl;
	cout<<"E= "<< endl;
	cout<<E<<endl;
	cout<<"E_0= "<< endl;
	cout<<E_0<<endl;
	cout<<"D_0= "<< endl;
	cout<<D_0<<endl;
*/
	/*cout<<"A0= "<< endl;
	cout<<A_0<<endl;
	cout<<"B_0= "<< endl;
	cout<<B_0<<endl;
	cout<<"D_0= "<< endl;
	cout<<D_0<<endl;
	cout<<"E_0= "<< endl;
	cout<<E_0<<endl;
*/

	SelfAdjointEigenSolver<MatrixXd> eD(D); 												//both D and E obviously self-adjoint as A^T=A, B^T=-B and both real.
	SelfAdjointEigenSolver<MatrixXd> eE(E);

	SelfAdjointEigenSolver<MatrixXd> eD_0(D_0);
	SelfAdjointEigenSolver<MatrixXd> eE_0(E_0);
	
	MatrixXd static_correlation = MatrixXd::Constant(L,L,0);

	MatrixXd psi_mat_0 = eD_0.eigenvectors();
	MatrixXd phi_mat_0 = eE_0.eigenvectors();

	G_0( &psi_mat_0, &phi_mat_0, static_correlation, L );								//static correlation matrix G_0 calculated with initial Hamiltonian.

	//calculate time-dependent contractions.
	MatrixXd psi_mat = eD.eigenvectors();
	MatrixXd phi_mat = eE.eigenvectors();
	VectorXd eigenval_v = eD.eigenvalues(); 												//this is the same vector as eE.eigenvalues().
	VectorXd eigenval_v2 = eE.eigenvalues();
	

	double T_end = 485165196;
	double dt 	 = 1.04;
	ofstream ENTROPY;
	ENTROPY.open("Garcia_test.txt"); //("Fibonacci_boudarytermsnominus.txt");
	for( double time = 0.05; time < T_end; time *= dt ) {
	
		//make the matrices used to calculate the correlation matrix.	
		MatrixXd A_A, B_B, A_B;
		A_A = A_B = B_B = MatrixXd::Constant(L,L,0);  												//B_A(L,L) is just A_B(l,L)^T

		Contractions_init(A_A, B_B, A_B, &phi_mat, &psi_mat, &eigenval_v, time, L);   				//calculate the three different time-dependent contractions.

		/*
		for(int i =0; i<L; i++) {
			for(int j=0; j<L; j++) {
				VectorXd v1 = psi_mat.col(i);
				VectorXd v_1 = D*v1;
				VectorXd v2 = phi_mat.col(i);
				VectorXd v_2 = E*v2;
				cout<< "test1=   " << v_1(j)/(eigenval_v(i)) - psi_mat.col(i)(j) <<endl;
				cout<< "test2=   " << v_2(j)/(eigenval_v(i)) - phi_mat.col(i)(j) <<endl;
			}
		}*/

		//make the correlation matrix. This is a SKEW-symmetric matrix
		MatrixXd correlation_mat = MatrixXd::Constant(2*L, 2*L, 0);
		Correlation_mat_create( &A_A, &B_B, &A_B, &static_correlation, correlation_mat, L );

		MatrixXd reduced_correlation_mat(2*sub_L,2*sub_L);

		for(int i = 2*233; i < 2*sub_L + 2*233; i++) {
				for(int j = 2*233; j < 2*sub_L + 2*233; j++) {
					reduced_correlation_mat(i-466,j-466) = correlation_mat(i,j);
				}
		}	

	/*	if( (L - sub_L) % 2 == 0 ) {
			for(int i = L-sub_L; i < L+sub_L; i++) {
				for(int j = L-sub_L; j < L+sub_L; j++) {
					reduced_correlation_mat(i-(L-sub_L),j-(L-sub_L)) = correlation_mat(i,j);
				}
			}	
		}
		else if( (L - sub_L) % 2 == 1) {
		//	cout<<"testing"<<endl;
			for(int i = L - sub_L + 1; i < L + sub_L + 1; i++) {
				for(int j = L - sub_L + 1; j < L + sub_L + 1; j++) {
					reduced_correlation_mat(i-(L - sub_L + 1),j-(L - sub_L + 1)) = correlation_mat(i,j);
				}
			}	
		} 
	*/
	/*MatrixXd reduced_correlation_mat(2*sub_L,2*sub_L);
	for(int i = 0; i < 2*sub_L; i++) {
		for(int j = 0; j < 2*sub_L; j++) {
			reduced_correlation_mat(i,j) = correlation_mat(i,j);
		}
	}	*/

	/*cout<<correlation_mat<<endl;
	cout<<"bla"<<endl;
	cout<<reduced_correlation_mat<<endl;
*/
	EigenSolver<MatrixXd> es(reduced_correlation_mat, false);
	
	VectorXd Correlation_eigen_vals = es.eigenvalues().imag();  //eigenvalues are all pure imaginery.
	//cout<<es.eigenvalues().real()<<endl;

	double entropy = Entropy(&Correlation_eigen_vals);
	
		ENTROPY << time << "\t" << entropy << endl;

	}
	ENTROPY.close();
}





