#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void function(double* f, double* y, double eta)
{
	f[0]=y[1];
	f[1]=(eta-pow(y[0],2))*y[0];	
}

void rk3(double*,double,const double,int);

int main()
{
	double xmax = 100;
	double dx;
	cout << "Enter dx= "<< endl;
	cin >> dx;
	//double dx = 0.1;
	int n = xmax/dx;
	double y[2];
	const double eta = 1/4;
	double initv = 1e-5;
	y[0]=initv;
	y[1]=sqrt(eta)*initv;
	double f[2];
	
	rk3(y,dx,eta,n);

}

void rk3(double* y,double dx,const double eta,int n )
{
	
	
	double k1[2],k2[2],k3[2];
	double ytemp[2];	
	
	ofstream output("data.txt");
	output << 0 << "\t" << y[0] << "\t" << y[1] << endl;
	
	for (int i = 1; i <=n; i++)
	{
	function(k1,y,eta);
		
	ytemp[0]=y[0]+(1./2.)*dx*k1[0];
	ytemp[1]=y[1]+(1./2.)*dx*k1[1];
	function(k2,ytemp,eta);
	
	
	ytemp[0]=y[0]+2.*k2[0]*dx+(-1.)*dx*k1[0];
	ytemp[1]=y[1]+2.*k2[1]*dx+(-1.)*dx*k1[1];
	function(k3,y,dx);
	
	y[0]=y[0]+(dx/6.)*(k1[0]+4*k2[0]+k3[0]);
	y[1]=y[1]+(dx/6.)*(k1[1]+4*k2[1]+k3[1]);
	output << i*dx << "\t" << y[0] << "\t" << y[1] << endl;
	}

	output.close();
}
