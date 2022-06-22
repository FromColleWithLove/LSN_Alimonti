/*Program that computes the autocorrelation function of a set of data from an
  input file at given times */
  #include<cmath>
  #include<iostream>
  #include<cstring>
  #include<fstream>
  using namespace std;
  void logspace(int xs[],int start, int stop, int npts);
  void linspace(int xs[],int start, int stop, int npts);
  
  int main(int argc, char *argv[]){
  
  string line;	
  char filename[50],filename2[50];
  double ac ;
  double norm;
  int pts= stoi(argv[2]),stop=stoi(argv[3]), size= stoi(argv[4]);
  int ts[pts]={0};
  int x,t,k;
  double y,z,c,dati[size];
  long double lagprodsum,reduced_sum,lag_sum=0.,sum=0., sum2=0.;
  sprintf(filename,"%s",argv[1]);
  strcpy(filename2,filename);
  strcat(filename,".dat");
  strcat(filename2,"_autocorr.dat");
  cout <<"Reading data from " <<filename <<endl;
  cout <<"Writing data on   " <<filename2 <<endl;
  
  ifstream data(filename);
  ofstream autocorr(filename2);
  
  k=0;
  
  linspace(ts,1,1500,pts);
  cout<<"uca"<<endl;

  for(int k=0;k<size;k++){
  data >> x >> y >> z >> c;
  sum2 += pow(y,2.);
  sum  += y;
  dati[k] = y;
}
cout<<ts[5];
cout<<"suca"<<endl;
 for(int i=0;i<pts;i++){
 	lag_sum = 0.;
 	lagprodsum = 0.;
 	reduced_sum=0.;
 	t= ts[i];
 	for(int j=0;j<(size-t);j++){
 		lagprodsum += dati[j]*dati[j+t];
 		reduced_sum += dati[j];
 		lag_sum += dati[j+t];
 	
 	}
 	ac = ((lagprodsum/double(size-t)) - (reduced_sum*lag_sum)/pow(double(size-t),2))/((sum2/double(size)) -pow((sum/double(size)),2.));
 	cout<<"Calculated " <<i+1 << " out of " << pts<<endl;
 	autocorr << t << " "<< ac <<endl; 		
 
 }
  //We generate an array of logarithmically spaced numbers
  //we need to feed the function an array full of zeroes
  
  
  
  
}

void logspace(int xs[],int start, int stop, int npts){
	double expi,expf,spacing;
	expi=log10((double)start);
	expf=log10((double)stop);
	spacing=(expf-expi)/double(npts-1);
	for(int i=0; i<npts; i++){
		xs[i]=int(pow(10,(expi + i*spacing)));
		cout<<xs[i]<<endl;
	}

	return;
}

void linspace(int xs[],int start, int stop, int npts){
	double spacing;
	spacing=(stop-start)/double(npts-1);
	for(int i=0; i<npts; i++){
		xs[i]= int(start + (i*spacing)) ;
		cout<<xs[i]<<endl;
	}

	return;
}


