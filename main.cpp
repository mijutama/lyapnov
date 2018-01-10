
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "rungekutta.h"

using namespace std;
using namespace Eigen;

// v to new v and get norm
void gramschmidt(Vector3d *v, double *norm)
{

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

}

int main(int argv, char* argc[])
{
	const int n = 60000;
	const double dt = 0.001;

	////////////////lorenz//////////////////
	Rungekutta<double> lorenz;
	double p=16.0, r=45.92, b=4.0;
	lorenz.v.push_back(1.0);
	lorenz.v.push_back(1.0);
	lorenz.v.push_back(1.0);
	lorenz.f.push_back([&p](vector<double> xs){return -p*xs[0]+p*xs[1];});
	lorenz.f.push_back([&r](vector<double> xs){return -xs[0]*xs[2]+r*xs[0]-xs[1];});
	lorenz.f.push_back([&b](vector<double> xs){return xs[0]*xs[1]-b*xs[2];});
	lorenz.dt = dt;
	///////////////////w////////////////////
	Eigen::Matrix3d Jacobi;
	Vector3d w[3];
	w[0]<<1,0,0;
	w[1]<<0,1,0;
	w[2]<<0,0,1;
	Rungekutta<Vector3d> rw[3];
	for(int i = 0; i < 3; i++)
	{
		rw[i].dt = dt;
		rw[i].v.push_back(w[i]);
		rw[i].f.push_back([&Jacobi](vector<Vector3d> xs){return Jacobi*xs[0];});
	}
	////////////////////////////////////////

	ofstream fp ("plot.data");
	double lyapnov[3]={0.0,0.0,0.0};
	double alyapnov[3]={0.0,0.0,0.0};

	fp << "# " << n << endl;

	for(double i = 0.0; i < n; i+=1.0)
	{
		rungekutta<double>(lorenz);
		Jacobi << -p,p,0.0, -lorenz.v[2]+r,-1.0,-lorenz.v[0], lorenz.v[1],lorenz.v[0],-b;	

		for(int j = 0; j < 3; j++)
			rw[j].v[0] = w[j];

		//// rungekutta method
		for(int j = 0; j < 3; j++)
		{
			rungekutta<Vector3d>(rw[j]);
			w[j] = rw[j].v[0];
		}

		gramschmidt(w, lyapnov);

		double t = (double)i*dt;

		// x y z lyapnov1 lyapnov2 lyapnov3 
		if(i != 0.0)
			fp << t << " " << lorenz.v[0] << " " << lorenz.v[1] << " " <<  lorenz.v[2] << " " << lyapnov[0]/t << " " << lyapnov[1]/t << " " << lyapnov[2]/t << endl;

		if(i >= 50000)
			for(int j = 0; j < 3; j++)
				alyapnov[j] += lyapnov[j]/t;


	}

	for(int j = 0; j < 3; j++)
		cout << alyapnov[j]/10000.0 << endl;

}
