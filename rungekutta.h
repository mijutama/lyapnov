
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
	vector<T> k[4];
	var.resize(n);
	for(int i = 0; i < 4; i++)
		k[i].resize(n);

	// 1
	for(int i = 0; i < n; i++)
		var[i] = t.v[i];
	for(int i = 0; i < n; i++)
		k[0][i] = t.dt*t.f[i](var);

	// 2
	for(int i = 0; i < n; i++)
		var[i] = t.v[i]+k[0][i]/2.0;
	for(int i = 0; i < n; i++)
		k[1][i] = t.dt*t.f[i](var);

	// 3
	for(int i = 0; i < n; i++)
		var[i] = t.v[i]+k[1][i]/2.0;
	for(int i = 0; i < n; i++)
		k[2][i] = t.dt*t.f[i](var);

	// 4
	for(int i = 0; i < n; i++)
		var[i] = t.v[i]+k[2][i];
	for(int i = 0; i < n; i++)
		k[3][i] = t.dt*t.f[i](var);
	
	for(int i = 0; i < n; i++)
		t.v[i] += (k[0][i]+2.0*k[1][i]+2.0*k[2][i]+k[3][i])/6.0;

}

#endif

