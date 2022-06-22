#include "asl.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;
void test(){
	cout<< "ASL: Alimonti Statistics Library  Copyright (C) 2022  Alimonti Davide."<<endl;
    cout<< "This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'. ";
    cout<< "This is free software, and you are welcome to redistribute it ";
    cout<< "under certain conditions."<<endl<<endl;
}

double avg(int d,double a[]){
	double y=0;
	for(int i=0;i<d; i++){
		y+= a[i];
	}
	y = y/float(d);
	return y;
}

double avg2(int d,double a[]){
	double y=0;
	for(int i=0;i<d; i++){
		y+= pow(a[i],2.);
	}
	y = y/float(d);
	return y;
}

double sq_mod(int d,double a[]){
	double r=0.;
	for(int i = 0; i < d; i++){
		r += a[i]*a[i];
		}
	return r;
	
}

double pearson_rho(int d,double a[],double b[]){
	double rho=0.,c[d]={0.};
	double sigmax=0.,sigmay=0.;
	
	//arrays for covariance
	for(int i=0;i<d;i++){
		c[i] = a[i]*b[i];
	}
	
	//compute sigmax,sigmay;
	sigmax = avg2(d,a)-pow(avg(d,a),2.);
	sigmay = avg2(d,b)-pow(avg(d,b),2.);
	sigmax = sqrt(sigmax);
	sigmay = sqrt(sigmay); 
	rho = (avg(d,c) - (avg(d,a)*avg(d,b)))/(sigmax*sigmay);

	return rho;
}
