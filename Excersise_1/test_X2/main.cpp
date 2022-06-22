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
#include <cmath>
#include <fstream>
#include <string>
#include "random.h"
#include "asl.h"

/*
 * Program that generates statistics from the RANNYU RNG in random.cpp
 * 
 * 
 * M :: # of intervals
 * N :: # of number/test
 * L :: # of tests
 *  
 * chi_squared test
 * 
 */
 

using namespace std;
 
int main (int argc, char *argv[]){
   Random rnd;
   int seed[4];
   int p1, p2;
   //don't touch before
   //GPL header for the ASL
   test();
   cout << " Chi squared test with :" <<endl;
   int const M= 100, N = 10000 , L = 100;
   int l;
   double r,chisq[L]={ 0. };
   //numbers that fall in the Mth interval during the Lth test!
   int n[M] = { 0 };
      
   cout<< "# of sub-intervals M=" << M << endl;
   cout<<"# of numbers generated per test = " << N <<endl;
   cout <<"# of tests ="  << L <<endl;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   ofstream resout("resout2.dat");
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
   /* we make a binning within several tries
    * i-th try , j-th number of the try
    * 
    * for each try we need the value of chisq
    */
   double z = float(N)/float(M);
   for(int i=0; i<L; i++){
	
		for(int j=0; j<M ; j++){
		n[j] = 0;
	}
	   for(int j=0 ; j<N; j++){
		   r = rnd.Rannyu();
		   l = floor( r*M);
		   n[l] ++;
		   cout <<"n=	"<<i*N + j <<"	"<< l << endl;
	   }
		   
	   //now the numbers have been generated for the ith try
	   // let's calculate the ith chisq
	   for(int t=0 ; t<M; t++){
		   chisq[i] += (pow((n[t] - z),2.))/z;
		   
	   }
	   resout << chisq[i]<<endl;
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
