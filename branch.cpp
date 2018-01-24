
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "Eigen/Dense"
#include "rungekutta.h"

#define digits numeric_limits<double>::digits10

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

vector<string> split(string s)
{
  vector<string> r;
  stringstream ss(s);
  string buf;
  while(getline(ss, buf, ' '))
    r.push_back(buf);
  return r;
}

string double_to_string(double d)
{
  ostringstream oss;
  oss << fixed <<  setprecision(digits) << d;
  return oss.str();
}

int branch(vector<vector<double>> &d, vector<double> &l, const int n, const double dt, double p, double r, double b)
{

	////////////////lorenz//////////////////
	Rungekutta<double> lorenz;
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

	double lyapnov[3]={0.0,0.0,0.0};
	double alyapnov[3]={0.0,0.0,0.0};

  double old_p[3];
  double old_v[3];
  bool first = true;
  ifstream ifs("data/lorenz_"+double_to_string(p)+"_"+double_to_string(r)+"_"+double_to_string(b)+".data");
  vector<double> ly;
  auto check = [](vector<double> &t, double tt){
    for(int i = 0; i < t.size(); i++)
      if(t[i] == tt)
        return false;
    return true;
  };
	double v[3];
  if(ifs.is_open())
  {
    string s;
    int c;
    for(c = 0;getline(ifs,s);c++)
    {
      vector<string> ss = split(s);
			for(int j = 0; j < 3; j++)
			{
				v[j] = stod(ss[j]);
				if(!first && old_v[j] > 0.0 && v[j] - old_p[j] <= 0.0)
					d[j].push_back(v[j]);
				old_v[j] = v[j]-old_p[j];
				old_p[j] = v[j];
				first = false;
			}

      for(int j = 0; j < 3; j++)
        alyapnov[j] += stod(ss[3+j]);
    }
      
    for(int i = 0; i < 3; i++)
      l.push_back(alyapnov[i]/(c-1));
    ifs.close();
    return 0;
  }
  else
  {

    ofstream ofs("data/lorenz_"+double_to_string(p)+"_"+double_to_string(r)+"_"+double_to_string(b)+".data");
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

      if(i >= n-10000)
      {
        for(int j = 0; j < 3; j++)
          alyapnov[j] += lyapnov[j]/t;

        for(int j = 0; j < 3; j++)
        {
          if(!first && old_v[j] > 0.0 && lorenz.v[j] - old_p[j] <= 0.0)
            d[j].push_back(lorenz.v[j]);
          old_v[j] = lorenz.v[j]-old_p[j];
          old_p[j] = lorenz.v[j];
          first = false;
        }
        // write
        ofs << fixed << setprecision(digits) << lorenz.v[0] << " " << lorenz.v[1] << " " << lorenz.v[2] << " " << lyapnov[0]/t << " " << lyapnov[1]/t << " " << lyapnov[2]/t << endl;
      }

    }

    for(int j = 0; j < 3; j++)
      l.push_back(alyapnov[j]/10000.0);
    ofs.close();
  }

  return 0;
}

int main()
{
	const int n = 60000;
	const double dt = 0.001;
  vector<vector<double>> d;
  vector<double> lyapnov;
  ofstream f("data/branch_p-x.data");
  d.resize(3);
  double p = 16.0;
  double r = 45.92;
  double b = 4.0;
  double &x = p;
	ofstream ff("lyapnov.data");
	double dx = 0.001;
	bool branchf = true;
  for(x = 39.8; x < 39.9;x +=dx)
  {
    branch(d,lyapnov,n,dt,p,r,b);
    // p max lyapnov0 lyapnov1 lyapnov2 
    for(int j = 0; j < d[0].size(); j++)
		{
      f << fixed << setprecision(digits) << x << " " << d[0][j] << endl;
			/* find chaos point
			if(j > 0)
			{
				if(d[0][j-1]+1.0 < d[0][j])
				{
					if(branchf)
					{
						cout << fixed << setprecision(digits)<< "branch:" << x-dx << "<p<" << x << endl;
						branchf = false;
					}
				}
			}
			*/
		}
		ff << fixed << setprecision(digits) <<  x << " " << r << " " << b << " " << lyapnov[0] << " " << lyapnov[1] << " " << lyapnov[2] << endl;
    
		if(lyapnov[0] < 0.01 && lyapnov[0] > -0.01 && lyapnov[1] < 0.01 && lyapnov[1] > 0.01)
			cout << fixed << p << " " << r << " " << b << endl;
    d[0].clear();
    d[1].clear();
    d[2].clear();
    lyapnov.clear();
  }

	f.close();
	ff.close();
}
