
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

int main()
{
	
	auto split = [](string s){
		
		vector<string> ret;
		string buf;
		stringstream  ss(s);
		while(getline(ss, buf, ' '))
			ret.push_back(buf);
		return ret;
	};
	string s;
	vector<string> ss;
	vector<string> xs, ys, zs;
	ifstream ifs("plot.data");
	ofstream xfs("xreat.data");
	ofstream yfs("yreat.data");
	ofstream zfs("zreat.data");

  getline(ifs,s);
	while(getline(ifs, s))
	{
		ss = split(s);
		xs.push_back(ss[1]);
		ys.push_back(ss[2]);
		zs.push_back(ss[3]);
	}

	ifs.close();

	// x
	int t1 = 50;
	int t2 = 50;
	for(int i = t1+t2; i+t1+t2 < xs.size(); i++)
		xfs << xs[i-t1-t2] << " " << xs[i-t2] << " " << xs[i] << endl;
	// y
	//t1 = 1;
	//t2 = 1;
	for(int i = t1+t2; i+t1+t2 < ys.size(); i++)
		yfs << ys[i-t1-t2] << " " << ys[i-t2] << " " << ys[i] << endl;
	// z
	//t1 = 1;
	//t2 = 1;
	for(int i = t1+t2; i+t1+t2 < zs.size(); i++)
		zfs << zs[i-t1-t2] << " " << zs[i-t2] << " " << zs[i] << endl;
	xfs.close();
	yfs.close();
	zfs.close();
	return 0;
}
