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
   //optimization: reducing M to a multiple of N
   int M = 1000000;
   double x, s1 = 0., s2 = 0., s10= 0., s100= 0.;
   
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   ofstream lorentz1("lorentz1.dat");
   ofstream lorentz2("lorentz2.dat");
   ofstream lorentz10("lorentz10.dat");
   ofstream lorentz100("lorentz100.dat");
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
     //throw standard dice
     
     x = rnd.Cauchy(1.);
     
     if ( i < 10000){
		 
		 lorentz1 << x << endl;
	 }
     
     s2 += x;
     s10 += x;
     s100 += x; 
     
     
     
     if (i > 0){
		 
		 if ((i+1)%2 == 0 and i < 20000){
			 lorentz2 << s2/2. << endl;
			 s2 = 0.;
		 }
		 if ((i+1)%10 == 0 and i < 100000){
			 lorentz10 << s10/10. << endl;
			 s10=0. ;
		 }
		 if ((i+1)%100 == 0){
			 lorentz100 << s100/100. << endl;
			 s100 = 0. ;
		 }
		 
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
