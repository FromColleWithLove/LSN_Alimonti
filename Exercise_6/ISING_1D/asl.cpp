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

double arraysum(int d,double a[]){
	double y=0;
	for(int i=0;i<d; i++){
		y+= a[i];
	}
	return y;	
}

double arraysum2(int d,double a[]){
	double y=0;
	for(int i=0;i<d; i++){
		y+= pow(a[i],2.);
	}
	return y;
}
double avg(int d,double a[]){
	double y=0;
	y = arraysum(d,a);
	y = y/float(d);
	return y;
}

double avg2(int d,double a[]){
	double y=0;
	y =arraysum2(d,a) ;
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
double autocorr_disc(double a[], int tmax,int t){
	//We take for granted that the array starts from 0 and goes to tmax
	double autocor , lag_prod_sum=0.,reduced_sum=0.,lag_sum=0.,sum=0., sum2=0.;
	double norm;
	
	norm = 1./double(tmax-t);
	sum = arraysum(tmax,a);
	sum2 = arraysum2(tmax,a);
	reduced_sum=arraysum(tmax-t,a);
	for(int i=0; i< tmax - t; i++){
		lag_prod_sum += a[i]*a[i+t];
		lag_sum += a[i+t];		
	}
	
	autocor = (norm*lag_prod_sum) - (pow(norm,2.)*reduced_sum*lag_sum);
	autocor /= (sum2/double(tmax)) - pow(sum/double(tmax),2.);
	return autocor;
}
