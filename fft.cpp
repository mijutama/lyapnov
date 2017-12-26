
#include <vector>
#include <fftw3.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main()
{

	auto split = [](const string &s)
	{
		vector<string> ret;
		stringstream ss(s);
		string sss;
		while(getline(ss, sss, ' '))
			ret.push_back(sss);
		return ret;
	};

	ifstream fs("plot.data");
	ofstream of("fft.data");
	string s;
	vector<string> r;

	// N
	getline(fs,s);
	r = split(s);
	int N = stoi(r[1]);

	fftw_complex *x, *xx;
	x  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
	xx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

	fftw_plan planx;
	planx = fftw_plan_dft_1d(N, x, xx, FFTW_FORWARD, FFTW_ESTIMATE);

	for(int i = 0; i < N && getline(fs,s); i++)
	{
		r = split(s);
		x[i][0] = stoi(r[0]);
		x[i][1] = 0.0;
	}

	fftw_execute(planx);
	for(int i = 0; i < N; i++)
		of << xx[i][0]*xx[i][0] << " " << xx[i][1]*xx[i][1] << endl;

	if(planx)
		fftw_destroy_plan(planx);
	fftw_free(x);
       	fftw_free(xx);

	return 0;
}
