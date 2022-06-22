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
 * takes argument from CL: M,N
 * M :: quantity of numbers generated
 * N :: length of block
 * 
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
   int M=stoi(argv[1]),N=stoi(argv[2]),L=int(M/N);
   //optimization: reducing M to a multiple of N
   if ( (M%N) != 0){
	   M = L*N;
   }
   double av[N]= { 0. },av2[N]= { 0. },x,sigma;
   double rq[N] = { 0. },rq2[N] = { 0. },sigmaq;
   
   cout<< "M=" << M << endl;
   cout<<"N= " << N <<endl;
   cout <<"L= "  << L <<endl;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   ofstream resout("resout.dat");
   ofstream resout2("resout2.dat");
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
   /* generatore*
    * We want to subdivide data in block and link some average values
    * to each block, namely the average and the average square 
    * 
    */
    int c=-1;
   for(int i=0; i<=M; i++){
     if (i%L == 0){
		 
		 if (i > 0){
			 //perform on-the-fly calculations and write on file
			 av[c] /= float(L);
			 rq[c] /= float(L);
			 av2[c] = pow(av[c],2.);
			 rq2[c] = pow(rq[c],2.);
			 if (c==0){
				sigma=0.;
				sigmaq=0.;
			 }
			else{
				sigma = pow((avg(c+1,av2) - pow(avg(c+1,av),2.))/float(c+1),.5);
				sigmaq = pow((avg(c+1,rq2) - pow(avg(c+1,rq),2.))/float(c+1),.5);
			}	
			resout << avg(c+1,av) <<","<<sigma <<endl;
			resout2 << avg(c+1,rq) <<","<<sigmaq <<endl;
	 }
		//increase counter
		
		 c += 1;
     }
	 if(i < M){
      x = rnd.Rannyu();
      av[c] += x;
      rq[c] += pow((x-.5),2.);
	}
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
