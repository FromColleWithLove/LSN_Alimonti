/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// excersise nr 2 by Alimonti Davide
 
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "asl.h"

using namespace std;
double g(double x), f(double x);

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   test();
   int M=stoi(argv[1]),N=stoi(argv[2]),L=int(M/N);
   double x,y;
   double sumy[N]= { 0 },sumy2[N]= { 0 },av,sigma2;
   //double const a=M_PI/2.;
   int c=-1;
   //optimization: reducing M to a multiple of N
   if ( (M%N) != 0){
	   M = L*N;
   }
   cout<< "M=" << M << endl;
   cout<<"N= " << N <<endl;
   cout <<"L= "  << L <<endl <<"**********************" <<endl;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   ofstream res2("res2.dat");
   ofstream imps("imps.dat");
   ofstream out_num("numbers.dat");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   for(int i=0; i<=M; i++){
     if (i%L == 0){
		 
		 if (i > 0){
			 //perform on-the-fly calculations and write on file
			 sumy2[c] = pow(sumy[c],2.);
			 av = avg(c+1,sumy);
			 if (c==0){
				 sigma2=0.;
			 }
			 else{
				 sigma2= (1./float(c))*(avg(c+1,sumy2)- pow(av,2.));
			}
			 res2 <<av <<","<<pow(sigma2,.5) <<endl;
	 }
		//increase counter
		 c += 1;
     }
	 if(i < M){
      x = rnd.Rej_Sampling(g,g(0),0.,1.);
      y= f(x)/g(x);
      //this performs brute force average integration (over a nice pdf)gg
      sumy[c] += y/float(L);
      //we now want to importance-sample with the taylor series
      
      
	}
   }
   rnd.SaveSeed();
   return 0;
   
}
double f(double x){
	double const a=M_PI/2.;
	return a * cos(a * x);
}
double g(double x){
	double const norm= 0.639500384868153;
	return ( 1 - (pow(M_PI*x,2.)/8.) + (pow(M_PI*x,4.)/384.) )/norm;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
