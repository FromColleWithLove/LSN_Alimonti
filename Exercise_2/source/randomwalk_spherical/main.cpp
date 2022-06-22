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

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   test();
   //don't touch before

   int M=10000 ,N=100, T=100,L=int(M/N);
   double pos[3];
   double ar2[T+1][N]={0.}, r2tot[T]={0.},r22tot=0,sigma2=0.;
   double theta,phi;
   cout<<"Discrete random walkerer"<<endl;
   cout<< "# of RWs=" << M << endl;
   cout<<"# of blocks= " << N <<endl;
   cout<<"# of steps =" << T <<endl;
   cout <<"length of blocks = "  << L <<endl <<"**********************" <<endl;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   ofstream r2disc("r2disc.dat");

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

   //cycling over (100) steps, but the first coordinate is surely (0,0,0)!
   
     //we need indication for each coordinate
     //we iterate over the walks and extract r2
     for(int j=0; j<M ;j++){
		for(int i=0; i<3;i++){
			pos[i]=0.;
		}
		for(int i=1; i<=T; i++){
			//simulate the ith step
			theta = rnd.Rannyu()* M_PI;
			phi = rnd.Rannyu()*2*M_PI;
			pos[0] += sin(theta)*cos(phi);	//x
			pos[1] += sin(theta)*sin(phi);	//y
			pos[2] += cos(theta);	//z
			//sum the sqmod of the ith step of the correct block
			ar2[i][j/L] += sq_mod(3,pos)/double(L);
		} 
		 
		
	 }
	 //we should by now have a table ar2 which contains sums of the sqmod
	 //for each  block, averaged per block. We now average over all the blocks
	 //we also determine the statistical uncertainty
   cout<<"lol";
   for(int i=0;i<=T;i++){
	   r22tot=0.;
	   sigma2=0.;
	   for(int j=0;j<N;j++){
	   r2tot[i] += ar2[i][j]/double(N);
	   //average of the squares(of the squares)
	   r22tot+= pow(ar2[i][j],2.)/double(N);

   }
   
   if(i<2){
	   sigma2=0;
	   cout<<"LOL";
   }
   else{
	   sigma2 = (r22tot - pow(r2tot[i],2.))/float(N-1);
   
   sigma2 = sigma2/(2. * sqrt(r2tot[i]));
   }
   r2disc <<i<<" "<<pow(r2tot[i],.5)<<" "<<sqrt(sigma2)<<endl;
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
