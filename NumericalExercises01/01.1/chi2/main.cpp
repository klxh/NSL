#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
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

	// chi2 test

	int N = 100; // Number of chi2 experiments
	double chi2[N] = {0}; 
	int M = 100; // Number of sub-intervals of [0,1]
	int n = 10000; // Number of throws for each chi2 test
	double delta = 1./M; // Sub-interval width

	ofstream myfile;
	myfile.open("./data/chi2.dat");
	//myfile.open("./data/chi2_M10.dat");

	for(int i = 0; i < N; i++){
		int n_i[M] = {0}; // Number of points drawn for each sub-interval
		for(int j = 0; j < n; j++){
			double r = rnd.Rannyu();
			double b = 0;
			int k = -1; // sub-interval index
			// Identifying the sub-interval for the point r
			while(r >= b){
				k++;
                // Since by construction of the PRNG r is always >= 0
                // the while condition is at first always satisfied.
                // Therefore at the end k is always at least 0.
				b += delta;
			}
			n_i[k]++;
		}
        // Computation of the chi2 for the i-th iteration
		for(int l = 0; l < M; l++){
			chi2[i] += (double) (n_i[l] - n/M) * (n_i[l] - n/M) / (n/M);
		}
		myfile << i + 1 << " " << chi2[i] << endl;
	}

	myfile.close();
}
