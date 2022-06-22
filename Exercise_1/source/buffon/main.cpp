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
   
   
   int const M=10000, L=100 , N= M/L;
   int C=0;
   double theta,y ,a[L] = { 0. },a2[L] = { 0. },P,var;
   double const r=1., d=2.01;
   
   
   test();
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   ofstream resout("buffon.dat");
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

   for(int i=0; i<=M; i++){
	   
	   if(i<M){
	   //sampling angles
	   //like Giovanni Muciaccia, i've already done that
	   theta = rnd.RandomAngle_4();
	   cout<< theta <<"	";
	   //now we define the distance of the center from the closest line
	   y =  (d/2.)*rnd.Rannyu();
	   //now we check if it intercepts the line
	   if (( y + r*sin(theta) > d ) or ( y - r*sin(theta) < 0.)){
		   C += 1;
	   }
	   cout<< y << "	" <<y + r*sin(theta) <<"	"<<  y - r*sin(theta) <<endl;

   }
	   if((i+1)%N == 0){
		   //blocking!
		   P = (4*r*N)/(C*d);
		   a[i/N]= P ;
		   a2[i/N] = P * P;
		   C=0;
		   if (i/N != 0){
			   var = (avg((i/N)+1,a2) - pow(avg((i/N)+1,a),2.))/float((i/N));
		   }
		   else if(i/N == 0){
			   var=0.;
			   
		   }
		   resout << avg((i/N)+1,a) << "," << pow(var,.5)<<endl;
		   resout2<< P <<"	"<<P*P <<endl;
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
