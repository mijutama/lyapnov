
#ifndef _RUNGEKUTTA_H_
#define _RUNGEKUTTA_H_

#include <initializer_list>
#include <vector>

using namespace std;

template<class T>
class Rungekutta
{
private:
public:
	double dt;
	vector<function<T(vector<T>)>> f;
	vector<T> v;
};

template<class T>
void rungekutta(Rungekutta<T> &t)
{
	const long n =  t.v.size();
	vector<T> var;
	vector<T> k[5];
	var.resize(n);
	for(int i = 0; i < 5; i++)
		k[i].resize(n);

  double a[4] {0.0, 0.5, 0.5, 1.0};
  for(int i = 1 ; i < 5; i++)
  {
    for(int j = 0; j < n; j++)
      var[j] = t.v[j]+a[i-1]*k[i-1][j];
    for(int j = 0; j < n; j++)
      k[i][j] = t.dt*t.f[j](var);
  }

	for(int i = 0; i < n; i++)
		t.v[i] += (k[1][i]+2.0*k[2][i]+2.0*k[3][i]+k[4][i])/6.0;

}

#endif

