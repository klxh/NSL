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
#include <math.h>
#include "random.h"

using namespace std;

#define CALL_TH 14.975790778311286 // Theoretical call price
#define PUT_TH 5.4595325819072364 // Theoretical put price
 
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

	// Option pricing sampling asset price path

	double S_0 = 100; // Asset price a t = 0
	double T = 1; // Delivery time
	double K = 100; // Strike price
	double r = .1; // Risk-free interest rate
	double sigma = .25; // Volatility
	int t = 100; // # time steps
	double dt = T/t; // Time interval between steps

	int M = 1000000; // MC steps
	int N = 100; // # of blocks
	int L = M/N;
	double C[N] = {0}; // Call option price for each block
	double C2[N] = {0}; // C squared

    // Loop for computing block variables
	for(int i = 0; i < N; i++){
		for(int j = 0; j < L; j++){
			double S = S_0;
			for(int l = 0; l < t; l++){
                // Sampling path to final price
				S *= exp((r - .5 * sigma * sigma) * dt + sigma * sqrt(dt) * rnd.Gauss(0,1));
			}
			if(S > K) { C[i] += exp(-r * T) * (S - K) / L; }
		}
		C2[i] = C[i] * C[i];
	}

    // Computing block averages and statistical errors

	double sum_prog[N] = {0};
	double sum2_prog[N] = {0};
	double err_prog[N] = {0};

	ofstream myfile1;
	myfile1.open("./data/call_path.dat");

	for(int i = 0; i < N; i++){
		for(int j = 0; j < i + 1; j++){
			sum_prog[i] += C[j] / (i+1);
			sum2_prog[i] += C2[j] / (i+1);
		}
		err_prog[i] = sqrt((sum2_prog[i] - sum_prog[i] * sum_prog[i]) / (i + 1));
		myfile1 << (i + 1) << " "
					 << sum_prog[i] - CALL_TH << " "
					 << err_prog[i] << endl;
	}

	myfile1.close();

    // Clearing old variables
	for(int i = 0; i < N; i++){
		sum_prog[i] = 0;
		sum2_prog[i] = 0;
		err_prog[i] = 0;
	}

	double P[N] = {0}; // Put option prices for each block
	double P2[N]; // P squared

	// Put option pricing
    // Repeating the same for put options

	for(int i = 0; i < N; i++){
		for(int j = 0; j < L; j++){
			double S = S_0;
			for(int l = 0; l < t; l++){
				S *= exp((r - .5 * sigma * sigma) * dt + sigma * sqrt(dt) * rnd.Gauss(0,1));
			}
			if(S < K) { P[i] += exp(-r * T) * (K - S) / L; }
		}
		P2[i] = P[i] * P[i];
	}

	ofstream myfile2;
	myfile2.open("./data/put_path.dat");

	for(int i = 0; i < N; i++){
		for(int j = 0; j < i + 1; j++){
			sum_prog[i] += P[j] / (i+1);
			sum2_prog[i] += P2[j] / (i+1);
		}
		err_prog[i] = sqrt((sum2_prog[i] - sum_prog[i] * sum_prog[i]) / (i + 1));
		myfile2 << (i + 1) << " "
					 << sum_prog[i] - PUT_TH << " "
					 << err_prog[i] << endl;
	}

	myfile2.close();

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
