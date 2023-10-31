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


	// Buffon's experiment

	double d = 1.; // Distance between lines
	double L = .8; // Needle length
	double N = 5000000; // Number of needle throws
	int M = 100; // Number of MC blocks
	double pi[M] = {0}; // Estimates of pi for each MC step
	double pi2[M] = {0}; // Estimates squared

	// I can sample angles uniformly on [0, 2pi] by uniformly sampling rays on a disk
	for(int k = 0; k < M; k++){
		double N_hits = 0;
        // N/M samplings for each block
		for(int i = 0; i < N/M; i++)
		{
			double x_c = rnd.Rannyu() * .5 * d; // Center of the needle drawn uniformly in [0, d/2]
			double dx, dy;
			double gamma;
			do{
				dx = rnd.Rannyu();
				dy = rnd.Rannyu();
				gamma = sqrt(dx*dx + dy*dy);
			} while(gamma > 1); // Sampling rays on a quarter of a disk; rejecting values outside the disk
			double x_t = x_c - .5 * L * dx / gamma; // x-coordinate of the top of the needle
                                                    // Note that dx / gamma is the cosine of the angle
                                                    // the needle forms with the x direction
			if(x_t < 0) // Intersection condition
				N_hits++;
		}
		pi[k] = 2 * L * N / (M * d * N_hits); // Pi estimate for k-th MC block
		pi2[k] = pi[k] * pi[k];
	}

	// Error evaluation
	double sum_prog[M] = {0};
	double sum2_prog[M] = {0};
	double err_prog[M] = {0};

	ofstream myfile;
	myfile.open("./data/pi.dat");

	for(int i = 0; i < M; i++){
		for(int j = 0; j < i + 1; j++){
			sum_prog[i] += pi[j];
			sum2_prog[i] += pi2[j];
		}
		sum_prog[i] /= i + 1;
		sum2_prog[i] /= i + 1;
		err_prog[i] = sqrt((sum2_prog[i] - sum_prog[i] * sum_prog[i]) / (i + 1));
		myfile << (i + 1) << " "
					 << sum_prog[i] - M_PI << " "
					 << err_prog[i] << endl;
	}

	myfile.close();

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
