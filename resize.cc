#include<iostream>
#include<Eigen/Dense>

using namespace Eigen;
using namespace std;
int main() {

	MatrixXd m(3,3);
	 m(0,0) = 0;
m(0,1) = 1;
m(0,2) = 2;
m(1,0)=3;
m(1,1)=4;
m(1,2)=5;
m(2,0)=6;
m(2,1)=7;
m(2,2)=8;

cout<<m<<endl;

cout<<m.col(1)(2)<<endl;
cout<<m(2,1)<<endl;


}

