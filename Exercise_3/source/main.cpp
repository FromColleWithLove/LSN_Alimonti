/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "asl.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   int M,N,L;
   double w,T;
   double sigma,mu,si,price,uncert=0.;
   cout<<"Insert # of prices generated"<<endl;
   cin >> M;
   cout<<"Insert # of blocks"<<endl;
   cin >> N;
   L= M/N;
   double c[N]={0.} , p[N]={0.};
   cout<<"***********************"<<endl;
   cout<<"# prices generated : "<<M<<endl;
   cout<<"# blocks           : "<<N<<endl;
   cout<<"length of blocks   : "<<L<<endl;
   cout<<"***********************"<<endl;
   cout<<"Insert T,mu,sigma,S(O)" <<endl;
   cin >>T>>mu>>sigma>>si;
   
   ifstream Primes("primes32001.in");
   ofstream call("rescall.dat");
   ofstream put("resput.dat");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
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

   for(int i=0; i<M; i++){
      //generate prices at time T
      w = rnd.Gauss( 0. , T);
      price= si*exp((mu- 0.5*sigma*sigma)*T + (sigma*w*sqrt(T)));
      c[i/L] +=exp(-mu*T) * max(0., price - si) / double(L) ;
      p[i/L] +=exp(-mu*T) * max(0., si - price) / double(L) ;
      cout<<i<<" "<<c[i/L]<<endl;
   }
   for(int i=0; i<N ; i++){
	   if(i>0){
		   uncert = (avg2(i+1,c) - pow(avg(i+1,c),2.))/double(i);
		   uncert = sqrt(uncert);
	   }
	   call<<avg(i+1,c)<<" "<<uncert<<endl;
	   uncert=0.;
	   if(i>0){
		   uncert = (avg2(i+1,p) - pow(avg(i+1,p),2.))/double(i);
		   uncert = sqrt(uncert);
	   }
	    put<<avg(i+1,p)<<" "<<uncert<<endl;
   }
   rnd.SaveSeed();
   return 0;
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
