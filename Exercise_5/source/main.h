

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include "asl.h"
#include "random.h"

//settings

Random rnd;
char filename[80],nametmp[20];
double ax=0.,ay=0.,az=0.;
double sigmax=0.,sigmay=0.,sigmaz=0.;
double xi,yi,zi;
int norb, GU,nblk,nstep,nstepeq;

//variables and observables

double r,theta,phi,x,y,z,xo,yo,zo;
double psio,psi,r2;
double accratebulk=0.,accrate=0.,accrate2bulk=0.,stimaccrate=0.,varaccrate=0.;
double rblk=0.,rglob=0.,r2glob=0.,stimr=0.,varr=0.,globav_r=0.,cont_r=0.;
int acc,att;



//functions
double psi1q(double r,double theta,double phi);
double psi2q(double r,double theta,double phi);
void Input(void);
void Move(void);
void Resetblk(void);
void Measure();
void Output(int ordinal);
