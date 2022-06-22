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
#include "random.h"
#include "asl.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   int M=stoi(argv[1]),N=stoi(argv[2]), L=(M/N);
   double a[N-2] = {0.}, b[N-2]={0.},c[N-2]={0.},x;
   cout<<"# numbers generated: "<<M<<endl;
   cout<<"# blocks			 : "<<N<<endl;
   cout<<"# per block        : "<<L<<endl;
   ifstream Primes("primes32001.in");
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
	   x = rnd.Rannyu();
	   if(i/L < N-1){
	   a[i/L] += x/double(L);
	   
   }
	   
	   if(i/L > 0){
		   b[(i/L)-1] += x/double(L);
	   }
	   if(i/L > 1){
		   c[(i/L)-2] += x/double(L);
	   }
	   
   }
   for(int i=0;i<N-2;i++){
	   cout<<a[i]<<" "<<b[i]<<" "<<c[i]<<endl;
   }
	cout << pearson_rho(N-2,a,b)<<endl;
	cout << pearson_rho(N-2,a,c)<<endl;
	
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
