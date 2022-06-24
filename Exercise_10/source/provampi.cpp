#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include "random.h"
#include "asl.h"
#include "mpi.h"

using namespace std;

int main(int argc,char *argv[]){
	int size,rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	cout<<"sono il "<<rank<<" su "<<size<<endl;
	
	int my_values[3];
	for (int i=0; i<3 ; i++)
		if(rank == 0) my_values[i]=i+1;
		else my_values[i]=0;
		cout<<"Prima :"<< my_values[0]<<" " << my_values[1] <<" " <<my_values[2] <<" per il processo "<<rank<<endl; 
	MPI_Bcast(my_values,3,MPI_INT,0,MPI_COMM_WORLD);
	cout<<"Dopo :"<< my_values[0]<<" " << my_values[1] <<" " <<my_values[2] <<" per il processo "<<rank<<endl; 
	
	//close program
	MPI_Finalize();
	return 0;
}
