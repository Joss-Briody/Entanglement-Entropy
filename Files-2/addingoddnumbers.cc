//Add first N odd numbers

#include <iostream>
#include <cmath>
using namespace std;

int main()
{
	float number,i, total=0.0;
	int N;
	cout<< "enter the integer, n for which the first n odd numbers will be summed"<<endl;
	cin>>N;

	for(int i=0; i<N; i++)
	{
	number=2*i+1;
	
	total=total+number;
	}
//output sum of odd numbers
cout<<total<<endl;
return 0;
}
