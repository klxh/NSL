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

	// Exercise 01.1

	int M = 1000000; // Number of samplings
	int N = 100; // Number of sampling blocks
	int L = M/N; // Number of samplings for each block
	double ave[N]; // Averages of sampled values for each block
	double ave2[N]; // Square averages

    // For each block I compute the block average and the square block average
	for(int i = 0; i < N; i++){
		double sum = 0;
		for(int k = 0; k < L; k++){
			sum += rnd.Rannyu(); // Adding up the sampled numbers
		}
		ave[i] = sum / L; // Block average
		ave2[i] = sum * sum / (L * L); // Square block average (the block average squared)
	}

	double sum_prog[N] = {0};
	double sum2_prog[N] = {0};
	double err_prog[N] = {0};

	ofstream myfile1;
	myfile1.open("./data/avg.dat");

	// Evaluating <r> and <r^2> as new block values become available
	for(int i = 0; i < N; i++){
		for(int j = 0; j < i + 1; j++){
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] /= i + 1; // Progressive block average
		sum2_prog[i] /= i + 1; // Progressive  block square average

        // Progressive block statistical error
	    err_prog[i] = sqrt((sum2_prog[i] - sum_prog[i] * sum_prog[i]) / (i + 1));

		myfile1 << (i + 1) * L << " "
					 << sum_prog[i] - .5 << " "
					 << err_prog[i] << endl;
	}

	myfile1.close();

	// Part 2, <(r-.5)^2>

	double sigma[N];
	double sigma2[N];

	for(int i = 0; i < N; i++){
		double sum = 0;
		for(int j = 0; j < L; j++){
			double a = rnd.Rannyu();
			a = (a - .5) * (a - .5);
			sum += a;
		}
		sigma[i] = sum / L;
		sigma2[i] = sigma[i] * sigma[i];
	}

	double sigma_prog[N] = {0};
	double sigma2_prog[N] = {0};
	double sigmaerr_prog[N] = {0};

	ofstream myfile2;
	myfile2.open("./data/sig.dat");

	for(int i = 0; i < N; i++){
		for(int j = 0; j < i + 1; j++){
			sigma_prog[i] += sigma[j];
			sigma2_prog[i] += sigma2[j];
		}
		sigma_prog[i] /= i + 1;
		sigma2_prog[i] /= i + 1;
		sigmaerr_prog[i] = sqrt((sigma2_prog[i] - sigma_prog[i] * sigma_prog[i]) / (i + 1));
		myfile2 << (i + 1) * L << " "
					 << sigma_prog[i] - (double) 1/12 << " "
					 << sigmaerr_prog[i] << endl;
	}

	myfile2.close();

  rnd.SaveSeed();
  return 0;
}
