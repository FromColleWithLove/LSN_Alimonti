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
   
   
   Random rnd;
   
   //simulation settings
   int const M=100000,L=200,N=M/L;
   int const nbins=200;
   double const bin_space = 8./(double)nbins;
   int const nstep = 120,ntemp= 0;
   double const mustep = 0.1, sigmastep=0.1,betastep=4.125;
   //common variables
   
   int bin[200]={0};
   
   
   int iblk,iann,itmp;
   
   double x,mu,sigma,muold,sigmaold;
   double V,psi,psi2,anterm,H_loc,H_tmp;
   double beta;
   double aH, glob_aH, glob_aH2, stima_aH,blk_aH,err_aH;
   double blk_norm;
   double accepted,attempted,annaccepted,annattempted;
   
   //annealing megaloop vars
   
   double Hmin,Hnew,Hold_ann,mumin,sigmamin,ErrHmin;
   
   //functions
   void Measure();   
   void Accumulate();
   void Reset(int);
   void Move();
   void Compare();
   void Averages(int);
   void MeasureEnergy();
   void AnnealMove();




int main (int argc, char *argv[]){
 
 
   mu=stof(argv[1]);
   sigma=stof(argv[2]);

   int seed[4];
   int p1, p2;


   ofstream Test("test.out");
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
   rnd.SaveSeed();
   
   Measure();
   MeasureEnergy();//starting en measure
   //setting Compare()
   mumin=mu;
   sigmamin=sigma;
   Hmin=Hnew;   
   for(iann=0;iann<(nstep*ntemp);iann++){//annealing megaloop
   beta = 5 + (iann/nstep)*betastep;    //setting beta
   AnnealMove();
   Compare(); //store params if H is minimum
   cout<<"Annealing move "<<iann+1<<" out of "<<nstep*ntemp<<". Actual beta is = "<<beta<<" Fictemp is = "<< 1/beta <<" it takes some time"<<endl;
   }
   
   cout<<"Enmin variazionale = "<<Hmin<<"   mu= "<<mumin<<"  sigma= "<<sigmamin<<endl;
   
   return 0;
}


void Measure(){
	//measures "everything"
	double s2= pow(sigma,2.),s4=pow(sigma,4.);
	double piu2= pow((x+mu),2.),meno2=pow((x-mu),2.);
	double d2psi;
	psi= exp(-piu2/2/s2) + exp(-meno2/2/s2);
	V = pow(x,4.) - 2.5*(pow(x,2));
	d2psi= (exp(-piu2/2/s2)*(piu2-s2)/s4)+(exp(-meno2/2/s2)*(meno2-s2)/s4);
	anterm= -d2psi/2/psi;
	H_tmp = anterm + V;
	return;
}
void Reset(int iblk){
	if(iblk==1){
		glob_aH2=0;
		glob_aH=0;
	}
	blk_aH=0;
	blk_norm =0;
	accepted = 0;
	attempted = 0;
}


void Accumulate(){
	blk_aH += H_tmp;
	blk_norm +=1;	
}
void Averages(int iblk){
	//modified for annealing
	
	stima_aH= blk_aH/blk_norm;
	glob_aH += stima_aH;
	glob_aH2 += pow(stima_aH,2.);
	ofstream H_out;
	ofstream sampsi;
	
	err_aH= Error(glob_aH,glob_aH2,iblk);
	
	if(iblk==N){
	sampsi.open("sampsi.out",ios::app);
	for(int i=0;i<nbins; i++){
		sampsi<< (-4+(bin_space/2.) + (double(i)*bin_space)) <<" "<< bin[i] <<endl;
	}
	
	}
	H_out.open("output_H.out",ios::app);
	Hnew = glob_aH/(double)iblk;
	H_out<< iblk << " "<<mu <<" "<<sigma<<" "<< Hnew <<" "<<err_aH<<endl;  	
	H_out.close();
	//}
	
	// /*
	cout << "Block number " << iblk << endl;
    	cout << "Acceptance rate " << accepted/attempted << endl << endl;
    	cout <<"*************************************************"<<endl;
    	//*/ //nice to have this uncommented if possible
}

void Compare(){//store minimum parameters and writes to text
	
	ofstream HT_out; 
	HT_out.open("output_HT.out",ios::app);
	if(Hnew < Hmin) sigmamin = sigma; //got to put it here because of mysterious behaviour
	if(Hnew < Hmin){//set new minimum
	Hmin = Hnew;
	mumin= mu;
	ErrHmin= err_aH;
	}
	HT_out << iann+1 <<" "<<beta<<" "<<Hmin<<" "<<ErrHmin<<endl;
	HT_out.close();
	return;
}

void Move(){
	double xold,psiold,hold, ptrans;
	xold = x;
	psiold = psi;
	hold = H_tmp;
	ofstream xs;
	xs.open("xs.out",ios::app);
	x += -2 + 4*rnd.Rannyu();
	Measure();
	ptrans = pow(psi,2.)/pow(psiold,2.);
	if (rnd.Rannyu() < ptrans){
	
	accepted +=1.;
	
	}

	else{
	x=xold;
	psi = psiold;
	H_tmp = hold;
	
	}
	attempted +=1.;
	//binning
	
	if(abs(x) < 4.) bin[int((x+4.)/bin_space)] +=1; //throw outside -4.4
	xs << x <<endl;
	xs.close();

	return;
}

void MeasureEnergy(){
//return a value for <H>;
x=0;
for(iblk=1;iblk <= N; iblk++){//given parameters
   Reset(iblk);
   
   
   for(int j=0;j<L;j++){//inside a block
   	Move();
   	Accumulate();
   	}
	
   Averages(iblk);   
   }
  return;

}

void AnnealMove(){
	
	double ptrans;
	//save old
	Hold_ann=Hnew;
	muold = mu;
	sigmaold = sigma;
	//set new params
	mu += -mustep + 2*mustep*rnd.Rannyu();
  	sigma += -sigmastep + 2*sigmastep*rnd.Rannyu();
  	mu = abs(mu);
  	sigma = abs(sigma);
	//measure energy
	MeasureEnergy();
	//decide if move
	ptrans = exp(-(Hnew-Hold_ann)*beta);
	if (rnd.Rannyu() > ptrans){//reject move,reinstate old params
		mu=muold;
		sigma=sigmaold;
		Hnew=Hold_ann;
	}//otherwise, keep the new params saved
	
	return;
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
