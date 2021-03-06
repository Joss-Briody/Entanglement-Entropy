//FFT_test
// Trivial FFT example using FFTW
#include <vector>
#include <iostream>
#include <cmath>
#include <fftw3.h>
int main()
{
  const int n = 800;
  std::vector<double> inp( 2 * n, 0 );
  std::vector<double> out( 2 * n, 0 );

  // Set the real parts to something quasi-random:
	for(int i = 0; i <inp.size(); i++)
		inp[i] = i;

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
  return 0;
}
