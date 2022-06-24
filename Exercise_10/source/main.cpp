/* Genetic Algorithm solver for the Traveling Salesman Problem
 * Ex nr. 9 of LSN course
 * By Davide Alimonti
 */


//libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include "random.h"
#include "asl.h"
#include "mpi.h"
//common variables


//initialization
int const n_cities=50,n_cromosoma=4500,n_generations=50000,n_migr=30;
int geni[n_cromosoma][n_cities], best_gene[n_cities], wbest_gene[n_cities];
int migrbuf1[n_cities],migrbuf2[n_cities];
int genitmp[n_cromosoma][n_cities];
int const starting_mutation=n_cromosoma*n_cities*5;
int const single_mutation = n_cromosoma*5;
char filename[40];
char nametmp[10];

//parallelization


//define the 2 types of norm
double l1[n_cromosoma]={0.};
double l2[n_cromosoma]={0.};
double lord[n_cromosoma];
//matrix with cities'coordinates
double cityxy[n_cities][2];
double beta = 0.1,bias=0;
double sqdist,bestl1,worstl1,tg_length,glob_bestl1,world_bestl1=0.;
bool popok;
int ch_index,nr_best,nr_best_world,igen,offspring_ptr;
int parent1,parent2;
int pop_best;
double ptried,pacc;
double lmin,lmax;
//functions,methods


void ImportMap();
void CheckPopulation();
void GeneratePopulation();
void MeasureLengths();
void CitiesSqDist(int,int);
void PrintChromosome(int);
void MutateChromosome(int);
void CrossOver();
void GenerateChromosome(int);
void OutputBest();
void PurgeUnfit();
void CosmicRain();
void OutputBestOverall(int);
void ParentSelector();
void ParentSelectorLin();
void OffspringReplace();
void JustCopy();
void MeasureBestHalf();
double Appeal(double);
using namespace std;
   
Random rnd;

int main (int argc, char *argv[]){
	  
	  int size,rank;

   //parallelizer inizializator
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   double gathl1[size];
   
   //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   
   //rng initializator
   int seed[4];
   int p1, p2;
   ifstream Primes("primes32001.in");
   if (Primes.is_open()){
      for(int i=0;i<=rank;i++) Primes >> p1 >> p2 ; //each thread has its primes
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
   	double tstart= MPI_Wtime();
   	parent1=floor(rnd.Rannyu()*n_cromosoma);
   	do parent2=floor(rnd.Rannyu()*n_cromosoma); while (parent1==parent2);
   	cout<<parent1<< " "<<parent2<<endl;
 
 	cout<<"PARALLELIZED VERSION"<<endl;
	//just setting up!
	string nameout="bestdist."+to_string(rank)+".out";
 	ofstream bdist(nameout);
 	ofstream wbestl1("world_bestdist.out");
 	sprintf(filename,"%s",argv[1]);
	strcpy(nametmp,".map.in");
	strcat(filename,nametmp);
	ImportMap(); 
	GeneratePopulation();
	cout<<"FIRST GENERATION"<<endl;
	//for(int i=0; i<n_cromosoma;i++) PrintChromosome(i);
	cout<<"Best  : "<<bestl1<<endl;
	glob_bestl1=bestl1;
	cout<<"Worst  : "<<worstl1<<endl;
	//now we start evolving
	
	/* 
	for(int igen=2;igen<=n_generations;igen++){
		PurgeUnfit(); //destroy unfit chromosome
		CosmicRain(); //mutate everything
		MeasureLengths();
		cout <<"GENERATION NR "<<igen<<endl;
		cout<<"Best overall : "<<glob_bestl1<<endl;
		cout<<"Best : "<<bestl1 <<endl;
		cout<<"Worst : "<<worstl1 <<endl;
		tg_length = (5*bestl1+worstl1)/6.;
		cout<<"Threshold next : "<<tg_length<<endl;
		cout<<"--------------------------------------------------"<<endl;
		bdist<<igen<<" "<<glob_bestl1<<endl;
	}*/ //Implementation with no crossing-over, just for showing
	
	for(int igen=2;igen<=n_generations;igen++){
		
		CheckPopulation();
		//start generating offspring
		offspring_ptr=0;
		do{
			ParentSelector();//choose 2 parents
			if (rnd.Rannyu()<0.8) CrossOver();//cross the 2 parents over
			else JustCopy();
			offspring_ptr +=2; //go on with repopulation!
		}while(offspring_ptr <= n_cromosoma-2);
		OffspringReplace();
		MeasureLengths(); //useful to measure
		CosmicRain(); //ozone layer who?
		CheckPopulation();
		MeasureLengths();
		MeasureBestHalf();
		//for(int i=0; i<n_cromosoma;i++) PrintChromosome(i);
		cout <<"GENERATION NR "<<igen<<" OF POP. NR "<<rank<<endl;
		cout<<"Best overall : "<<glob_bestl1<<endl;
		cout<<"Best : "<<bestl1 <<endl;
		cout<<"Worst : "<<worstl1 <<endl;
		cout<<"Acceptance rate so far: "<<pacc/ptried<<endl;
		cout<<"--------------------------------------------------"<<endl;
		bdist<<igen<<" "<<glob_bestl1<<endl;
		
		if(igen%n_migr == 0 and size>1){//è tempo di migrare
			int cont1,cont2;
			if(not rank){//let p0 select two random continents
			 cont1=floor(rnd.Rannyu()*size);
			do  cont2=floor(rnd.Rannyu()*size); while(cont2==cont1);
			}
			//cont 1, cont 2 to exchange their best
			MPI_Bcast(&cont1,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(&cont2,1,MPI_INT,0,MPI_COMM_WORLD);
			
			if(rank==cont1){
				for(int i=0;i<n_cities;i++){
					migrbuf1[i]=best_gene[i];
				}
			}
			if(rank==cont2){
				for(int i=0;i<n_cities;i++){
					migrbuf2[i]=best_gene[i];
				}

			}
			MPI_Bcast(migrbuf1,n_cities,MPI_INT,cont1,MPI_COMM_WORLD);
			MPI_Bcast(migrbuf2,n_cities,MPI_INT,cont2,MPI_COMM_WORLD);
			if(rank==cont1){//on cont1, put the best of cont2
				for(int i=0;i<n_cities;i++){
					geni[nr_best][i]=migrbuf2[i];
				}
			}
			if(rank==cont2){//on cont2, put the best of cont1
				for(int i=0;i<n_cities;i++){
					geni[nr_best][i]=migrbuf1[i];
				}
			}
			
			//migration hopefully completed
		}
		MeasureLengths();//update
		 
		MPI_Gather(&glob_bestl1,1,MPI_DOUBLE,gathl1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		//we now have all gathered up and nice! let's pick the best
		if (not rank){//core 0 will take care of this

			//we find the minimum, and find the corresponding index
			//then we'll let each continent know
			world_bestl1=arraymin(size,gathl1);
			for(int i=0;i<size;i++) if(gathl1[i]==world_bestl1) pop_best=i;
			}
			MPI_Bcast(&world_bestl1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&pop_best,1,MPI_INT,0,MPI_COMM_WORLD);
		//comunicate the best chromosoma to the world
		if (rank==pop_best)for(int i=0;i<n_cities;i++) wbest_gene[i]=best_gene[i];
		MPI_Bcast(wbest_gene,n_cities,MPI_INT,pop_best,MPI_COMM_WORLD);
		//hopefully, now each continent know who's the best gene and how
		// long it has it
		//in core 0, write the best of world length
		if (not rank) wbestl1<<igen<<' '<<world_bestl1<<endl;
		
	}
	
 OutputBestOverall(rank);
 bdist.close();
 rnd.SaveSeed();
 wbestl1.close();
 cout<<"Time taken : "<<MPI_Wtime() - tstart<<endl;
 cout<<pacc/ptried<<endl;
 MPI_Finalize();  
 return 0;
 }

//FUNCTIONS




void ImportMap(){//takes a map from .in file and puts it into cityxy 

//CAREFUL... check # of cities!
	ifstream MapIn(filename);
	for (int i=0;i<n_cities;i++){
		MapIn >> cityxy[i][0] >> cityxy[i][1];
		cout << cityxy[i][0]<<" "<<cityxy[i][1]<<endl;	
	}
	cout <<"SETTINGS: "<<endl;
	cout <<"Nr. of cities : "<<n_cities<<"  Nr. of chromosoma : "<<n_cromosoma <<endl;
	cout <<"Reading map from: " << filename <<endl; 
	MapIn.close();
	return;


}
void CheckPopulation(){ //returns true if population satisfies bonds
	//idea: start from first city of each cromosome
	//if we find the same again, problem, otherwise np!
	popok = true;
	int tmp,buf;	
	for(int i=0;i<n_cromosoma;i++){ //over cromosoma
		for(int j=0; j<n_cities-1;j++){//over geni
			tmp=geni[i][j]; //store city
			for(int k=j+1;k<n_cities;k++){//check next geni
				if(tmp == geni[i][k] or geni[i][k] >= n_cities or geni[i][k] < 0) popok=false; //problem
			}
		}
	}

	if (not popok){
	 cout<< " !! PROBLEM! Population non satisfying bonds !!"<<endl; 
	 cin >> buf;
	}
	return;
}

void GeneratePopulation(){
	int sel;
	for(int i=0;i<n_cromosoma;i++){ //cycle over cromosoma
		//we will generate ordered sequencies of cities,
		// then we'll mix them with suitable mutation operators
		//this ensures we can have a working population easily
		for(int j=0;j<n_cities;j++){
			geni[i][j]=j;
		}
	}
	for(int i=0;i<starting_mutation;i++){//initial mutation
		sel = floor(rnd.Rannyu()*n_cromosoma);
		MutateChromosome(sel);
	
	}
	CheckPopulation();
	MeasureLengths();
	return;
	
}

void GenerateChromosome(int index){//generate a single chromosome at index "index"
	for(int j=0;j<n_cities;j++) geni[index][j]=j;
	for(int i=0;i<single_mutation;i++) MutateChromosome(index);
	
	return;
	
}

void MeasureLengths(){ //measures l1, l2 lengths of paths and returns best (using l1 norm)
	int oldn,newn,htown;
	double newl1;
	for(int i=0;i<n_cromosoma;i++){
		l1[i]=0.;
		l2[i]=0.;
		oldn=geni[i][0];//stores first city
		htown=oldn;
		for (int j=1 ; j<n_cities ; j++){
			newn=geni[i][j];//reads next city
			CitiesSqDist(newn,oldn);
			l2[i]+= sqdist;
			l1[i]+= sqrt(sqdist);
			oldn=newn; //switches old<->new
		}
		//we need to mesure the last trip to return to hometown
		CitiesSqDist(htown,oldn);
			l2[i]+= sqdist;
			l1[i]+= sqrt(sqdist);
			newl1 = l1[i]; 
			if(not i){
			 bestl1 = newl1;
			 worstl1 = newl1;
			 nr_best = 0;
			 }
			if(newl1 < bestl1){
			bestl1 = newl1;
			nr_best= i;
			}
			if(newl1 > worstl1) worstl1 = newl1;
			if(newl1 <= glob_bestl1){
			 for(int k=0;k<n_cities;k++) best_gene[k]=geni[i][k];
			 glob_bestl1 = newl1;
			 } 
	}
	return;
}

void CitiesSqDist(int c1,int c2){//measures square distance between cities c1 and c2
	sqdist=0;
	sqdist += pow((cityxy[c1][0]-cityxy[c2][0]),2.);
	sqdist += pow((cityxy[c1][1]-cityxy[c2][1]),2.);
	return; 
	
}

void PrintChromosome(int index){//print a given cromosome
	cout<<"[ ";
	for(int i=0;i<n_cities;i++){
	cout<<geni[index][i]<<" ";
	}
	cout<<"] l1: "<<l1[index]<<" l2: "<<l2[index]<<endl;
	return;

}

void MutateChromosome(int index){//mutates a chromosome by switching 2 cities
	
	int loc1,loc2;
	int tmp;
	//select 2 loci excepting the 0-th 
	loc1 = 1+floor(rnd.Rannyu()*(n_cities-1));
	loc2 = 1+floor(rnd.Rannyu()*(n_cities-1));
	//switch the corresponding geni
	tmp=geni[index][loc1];
	geni[index][loc1]=geni[index][loc2];
	geni[index][loc2]=tmp;
	return;
}
void CrossOver(){//performs crossing over of 2 chromosoma
	int tmp1[n_cities],tmp2[n_cities];
	int xloc,A1,A2;
	xloc = 1+floor(rnd.Rannyu()*n_cities);	
	//first of all copy chr.. in buffer
	for (int i=0;i<n_cities;i++){
		tmp1[i]=geni[parent1][i];
		tmp2[i]=geni[parent2][i];
		//cout<<tmp1[i]<<" "<<tmp2[i]<<endl;
		if(i >= xloc){
			geni[parent1][i] = 0;
			geni[parent2][i] = 0;
		}
	}
	//write over parent1
	for(int i=0;i<n_cities;i++){
	A2 = tmp2[i];
	int k=1;
	do{
		A1=geni[parent1][k];
		if(not A1) geni[parent1][k]=A2;
		k++;
	}while(A1 and (A1 !=A2)); //cycle until gene itself is found or empty locus
	
	}
	//write over parent2
	for(int i=0;i<n_cities;i++){
	A1 = tmp1[i];
	int k=1;
	do{
		A2=geni[parent2][k];
		if(not A2) geni[parent2][k]=A1;
		k++;
	}while(A2 and (A1 !=A2)); //cycle until gene itself is found or empty locus
	
	}
	//now the geni are recombined... let's put them into genitmp and reset the geni
	for(int i=0;i<n_cities;i++){
		genitmp[offspring_ptr][i] = geni[parent1][i];
		genitmp[offspring_ptr+1][i] = geni[parent2][i];
		geni[parent1][i] = tmp1[i];
		geni[parent2][i] = tmp2[i];
	}
	return;
}

void OutputBest(){//takes shortest route and prints to file (uses l1 norm)
	int cityid;
	ofstream xybest;
	xybest.open("xybest.out");
	for(int i=0; i<n_cities; i++){
		cityid=geni[nr_best][i];
		xybest << cityxy[cityid][0] <<' '<<cityxy[cityid][1]<<endl;
	}
	xybest.close();
	

	return;
}
void OutputBestOverall(int thread_id){//takes shortest route so far and prints to file (uses l1 norm)
	int cityid;
	string nametmp;
	//strcat(st,"out");
	nametmp="xybest." + to_string(thread_id) + ".out";
	ofstream xybest;
	xybest.open(nametmp);
	for(int i=0; i<n_cities; i++){
		cityid=best_gene[i];
		xybest << cityxy[cityid][0] <<' '<<cityxy[cityid][1]<<endl;
	}
	//add trip back home
	xybest << cityxy[0][0] <<' '<<cityxy[0][1]<<endl;
	xybest.close();
	

	return;
}

void PurgeUnfit(){ //replaces all chromosoma under certain fitness with new chromosoma

	for(int i=0;i<n_cromosoma;i++){
		if(l1[i]>tg_length) GenerateChromosome(i);
	}
	MeasureLengths();
}
void CosmicRain(){ //mutates a little bit of eveything
	double mutfun;
	//loop over chromosoma and mutate them with some prob
	for (int i=0;i<n_cromosoma;i++){
		mutfun = (l1[i]-bestl1)/(worstl1-bestl1);
		for (int j=0;j<n_cities;j++){
			//choose whether to mutate
			
			if (rnd.Rannyu() <= 0.10) MutateChromosome(i);
		}
	}
	
	return;	
}
void ParentSelector(){//metropolis choice of new parents
	int newp,check;
	double p;
	
	check=0;
	newp = floor(rnd.Rannyu()*n_cromosoma);
	p = exp(-beta*(l1[newp]+bias-l1[parent1]));
	if (rnd.Rannyu() <= p) parent1 = newp;
	do{
	newp = floor(rnd.Rannyu()*n_cromosoma);
	p = exp(-beta*(l1[newp]-l1[parent2]));
	if (rnd.Rannyu() <= p) parent2 = newp;
	check++;
	ptried+=1;
	if (check > 1000*n_cromosoma) {
	cout<<"WARNING:PARENTS TRIED SURPASSED NR CHROMOSOMA*1000 TRY REDUCE BETA"<<endl;
	cout<<"SWITCHING TO RANDOM SELECTION"<<endl;
	pacc-=1;
	parent2=newp;//problem
	}
	}while(parent2==parent1);
	pacc+=1;
	//if(parent2==parent1) parent2= (parent1+2)%n_cromosoma;
	return;
}

void ParentSelectorLin(){//metropolis choice of new parents
	int newp,check;
	double p;
	
	check=0;
	newp = floor(rnd.Rannyu()*n_cromosoma);
	p = Appeal(l1[newp])/Appeal(l1[parent1]);
	if (rnd.Rannyu() <= p) parent1 = newp;
	do{
	newp = floor(rnd.Rannyu()*n_cromosoma);
	p = Appeal(l1[newp])/Appeal(l2[parent2]);
	if (rnd.Rannyu() <= p) parent2 = newp;
	check++;
	ptried+=1;
	if (check > 1000*n_cromosoma) {
	cout<<"WARNING:PARENTS TRIED SURPASSED NR CHROMOSOMA*1000 TRY REDUCE BETA"<<endl;
	cout<<"SWITCHING TO RANDOM SELECTION"<<endl;
	pacc-=1;
	parent2=newp;//problem
	}
	}while(parent2==parent1);
	pacc+=1;
	//if(parent2==parent1) parent2= (parent1+2)%n_cromosoma;
	return;
}

double Appeal(double l){
	double ap;
	ap=l-lmin;
	ap = ap/(lmax-lmin);
	ap = 1.2-ap;
	return ap;
}	



void OffspringReplace(){ //substitutes last generation with new one
	for(int i=0;i<n_cromosoma;i++){
		for(int j=1;j<n_cities;j++){
			geni[i][j]=genitmp[i][j];
			genitmp[i][j]=0;
		}
	}
	return;
}
void JustCopy(){//copies geni parent1,parent2 in the offspring
	for(int i=0;i<n_cities;i++){
		genitmp[offspring_ptr][i]=geni[parent1][i];
		genitmp[offspring_ptr+1][i]=geni[parent2][i];
	}

	return;
}

void MeasureBestHalf(){//takes measure from the best half of the population
	ofstream outl1;
	outl1.open("besthalf_l1.out",ios::app);
	double tmpl[n_cromosoma];
	double tmin,buf;
	double a,a2,err;
	for(int i=0; i<n_cromosoma;i++) tmpl[i]=l1[i];
	for(int i=0; i<n_cromosoma;i++){
		tmin=arraymin(n_cromosoma,tmpl);
		lord[i]=tmin;
		int j=0;
		do{
			buf=tmpl[j];
			if(buf==tmin) tmpl[j]=500000;
			j++;

		}while(buf != tmin);
	}
	a=avg(n_cromosoma/2 , lord);
	a2=avg2(n_cromosoma/2,lord);
	err= a2 - pow(a,2.);
	err=sqrt(err);
	outl1 << a <<" "<< err <<endl;
	outl1.close();
	return;
}
