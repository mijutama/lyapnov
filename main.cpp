
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include "Eigen/Dense"
#include "rungekutta.h"

using namespace std;
using namespace Eigen;

// v to new v and get norm
void gramschmidt(Matrix3d &m, double *norm)
{

  Vector3d v[3];
  v[0] << m(0,0),m(1,0),m(2,0);
  v[1] << m(0,1),m(1,1),m(2,1);
  v[2] << m(0,2),m(1,2),m(2,2);
	Vector3d vv;
	for(int i = 0; i < 3; i++)
	{
		vv << 0.0,0.0,0.0;
		for(int j = 0; j < i; j++)
			vv += (v[i].dot(v[j])/(v[j].norm()*v[j].norm()))*v[j];
		v[i] -= vv;		
		norm[i] += log(v[i].norm());
	}
	for(int i = 0; i < 3; i++)
		v[i] = v[i]/v[i].norm();

  m << v[0], v[1], v[2];
}

int main(int argv, char* argc[])
{
	const int n = 60000;
	const double dt = 0.001;

	////////////////lorenz//////////////////
	Rungekutta<double> lorenz;
	double p=39.8, r=45.92, b=4.0;
	lorenz.v.push_back(1.0);
	lorenz.v.push_back(1.0);
	lorenz.v.push_back(1.0);
	lorenz.f.push_back([&p](vector<double> xs){return -p*xs[0]+p*xs[1];});
	lorenz.f.push_back([&r](vector<double> xs){return -xs[0]*xs[2]+r*xs[0]-xs[1];});
	lorenz.f.push_back([&b](vector<double> xs){return xs[0]*xs[1]-b*xs[2];});
	lorenz.dt = dt;
	///////////////////w////////////////////
	Eigen::Matrix3d Jacobi;
	Matrix3d w;
	w<<1,0,0, 0,1,0, 0,0,1;
	Rungekutta<Matrix3d> rw;
  rw.dt = dt;
  rw.v.push_back(w);
  rw.f.push_back([&Jacobi](vector<Matrix3d> xs){return Jacobi*xs[0];});
	////////////////////////////////////////

	ofstream fp ("plot.data");
	double lyapnov[3]={0.0,0.0,0.0};
	double alyapnov[3]={0.0,0.0,0.0};

	fp << "# " << n << endl;

	for(int i = 0; i < n; i++)
	{
		rungekutta<double>(lorenz);
		Jacobi << -p,p,0.0, -lorenz.v[2]+r,-1.0,-lorenz.v[0], lorenz.v[1],lorenz.v[0],-b;	

    rw.v[0] = w;

		//// rungekutta method
    rungekutta<Matrix3d>(rw);
    w = rw.v[0];

		gramschmidt(w, lyapnov);

		double t = (double)i*dt;

		// x y z lyapnov1 lyapnov2 lyapnov3  
		if( i > 0)
		{
			fp << fixed << setprecision(numeric_limits<double>::max_digits10) << t << " " << lorenz.v[0] << " " << lorenz.v[1] << " " <<  lorenz.v[2] << " " << lyapnov[0]/t << " " << lyapnov[1]/t << " " << lyapnov[2]/t << endl;
		}

		if(i >= n-10000)
			for(int j = 0; j < 3; j++)
				alyapnov[j] += lyapnov[j]/t;


	}

	for(int j = 0; j < 3; j++)
		cout << alyapnov[j]/10000.0 << endl;

}
