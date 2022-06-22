/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//ex. nr. 5, sampling orbital psi with metropolis
//we use a0 units

#include "main.h"

using namespace std;
int main (int argc, char *argv[]){
   double tests;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   ofstream pos("out_pos.dat");
   ofstream testout("out_test.dat");
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
	

   cout<<"Insert the name of the input file, of the form input.name"<<endl;
  cin >> filename;
  strcpy(nametmp,".in");
  strcat(filename,nametmp);
  cout <<"Using " << filename <<" as input file"<<endl;
  
  for(int l=0;l<20;l++){
  
  Input(); //Inizialization
  Resetblk();
  tests= 2.2 + l*0.05;
  sigmax= tests;
  sigmay= tests;
  sigmaz= tests;
  
 //fool cycles to equilibrate
  for(int i=0;i<nstepeq;i++){
	  Move();
  }

	  
  for(int i=0; i<nblk; i++){
	   for(int j=0;j<nstep;j++){
	  Move();
	  Measure();
  }
	  Output(i+1);
	  testout<< tests <<" "<<double(acc)/double(att)<<endl;
	  Resetblk();
  }
  
}

   rnd.SaveSeed();
   return 0;
}

double psi1q(double x,double y,double z){
	//pdf of psi1 orbital
	double p;
	r = sqrt(x*x+y*y+z*z);
	theta = acos(z/r);
	phi= acos(x/(r*sin(theta)));
	p= exp(-2. * r)/M_PI;
	return p;
}	
double psi2q(double x,double y,double z){
	//pdf of psi2 orbital
	double p;
	r = sqrt(x*x+y*y+z*z);
	theta = acos(z/r);
	phi= acos(x/(r*sin(theta)));
	p=(pow(r,2.)*exp(-r)*pow(cos(theta),2.))/(32.*M_PI);
	return p;
}	

void Input(){
	ifstream ReadInput;
	ReadInput.open(filename);
	cout<<"Metropolis sampling of orbital functions"<<endl;
	cout<<"Reading from "<<filename <<endl;
	ReadInput >> norb;
	if (norb==0){
		cout<<"Orbital type: Psi 1,0,0"<<endl;
	}
	if (norb==1){
		cout<<"Orbital type: Psi 2,1,0"<<endl;
	}
	ReadInput >> GU;
	if (GU==0){
		cout<<"Uniform transition probability"<<endl;
		ReadInput >> ax;
		ReadInput >> ay;
		ReadInput >> az;
		cout <<"Steps are: "<<ax<<","<<ay<<","<<az<<endl;
	}
	if (GU==1){
		cout<<"Multivariate normal transition probability"<<endl;
		ReadInput >> sigmax;
		ReadInput >> sigmay;
		ReadInput >> sigmaz;
		cout <<"Variances are: "<<sigmax<<","<<sigmay<<","<<sigmaz<<endl;
	}
	ReadInput>> xi;
	ReadInput>> yi;
	ReadInput>> zi;
	x=xi;
	y=yi;
	z=zi;
	ReadInput>> nblk;
	ReadInput>> nstep;
	ReadInput>> nstepeq;
	if(norb==0){
		psio=psi1q(x,y,z);	
	}
	if(norb==1){
		psio=psi2q(x,y,z);
	}
	acc=0;
	att=0;
	
}
void Move(void){
	double probtr;
	xo=x;
	yo=y;
	zo=z;
	if(GU==0){
	//uniform move	
	x= x + (ax*rnd.Rannyu()) - (ax/2.);
	y= y + (ay*rnd.Rannyu()) - (ay/2.);
	z= z + (az*rnd.Rannyu()) - (az/2.);
}
if(GU==1){
	//gaussian move	
	x= x + (rnd.Gauss(0.,0.));
	y= y + (rnd.Gauss(0.,sigmay));
	z= z + (rnd.Gauss(0.,sigmaz));
}
	if(norb==0){
		psi= psi1q(x,y,z);
	}
	if(norb==1){
		psi= psi2q(x,y,z);
	}
	probtr=psi/psio;
	if(rnd.Rannyu() < probtr){
		//accept the move
		xo=x;
		yo=y;
		zo=z;
		psio=psi;
		acc +=1;
	}
	else{
		//reject the move
		x=xo;
		y=yo;
		z=zo;
		psi=psio;
	}
	att +=1;
}
void Resetblk(void){
	acc=0;
	att=0;
	accrate=0;
	cont_r=0.;
	
}

void Measure(){
	
	r= sqrt(x*x+y*y+z*z);	
	
	cont_r += r;
	
	
	rglob += r;
	
	
	
}
void Output(int ordinal){
	
   ofstream out_accrate;
   ofstream out_r;
   
   accrate = double(acc)/double(att);
   
   
   accratebulk += accrate;
   
   accrate2bulk += pow(accrate, 2.);
   
   stimaccrate = accratebulk/double(ordinal);
   stimr=cont_r/double(nstep);
   
   r2glob += stimr*stimr;
   globav_r = rglob/double(nstep*ordinal);
   
   if(ordinal>1){
	   varaccrate=((accrate2bulk/double(ordinal))-pow(accratebulk/double(ordinal),2))/double(ordinal-1);
	   varr = ((r2glob/double(ordinal))-pow(globav_r,2.))/double(ordinal-1);
	  }
	else{
		varaccrate=0.;
		varr=0.;
	}
   
   out_accrate.open("out_accrate.dat",ios::app);
   out_r.open("out_r.dat",ios::app);
	
	out_accrate <<ordinal << " "<< accrate <<" "<< stimaccrate <<" "<<sqrt(varaccrate)<<endl;
	out_r << ordinal <<" "<<globav_r<<" "<<stimr<<" "<<sqrt(varr)<<endl;
	
	out_accrate.close();
	out_r.close();
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
