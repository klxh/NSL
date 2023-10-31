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

    // RW on a cubic lattice

    int T = 101; // Observation time
    double a = 1.0; // Lattice constant
    double M = 10000; // Number of RWs in the ensemble
    int N = 100; // Number of blocks
    int L = M/N; // Number of RWs in each block
    double MSD[N][T] = {0}; // Mean Square Displacement [block][time step]

    for(int i = 0; i < N; i++){
        double RWL[L][3] = {0}; // Positions in 3D space of each RW in the block
        for(int t = 1; t < T; t++){
            for(int j = 0; j < L; j++){
				int jump_dir = floor(6*rnd.Rannyu()); // Sampling the direction of the step
				double jump[3] = {0};

				switch(jump_dir){
					case 0: jump[0] += a; break;
					case 1: jump[0] -= a; break;
					case 2: jump[1] += a; break;
					case 3: jump[1] -= a; break;
					case 4: jump[2] += a; break;
					case 5: jump[2] -= a; break;
				}
                // Updating position of RW
				RWL[j][0] += jump[0];
				RWL[j][1] += jump[1];
				RWL[j][2] += jump[2];				

                // Contribution to ensemble average of square displacement of j-th RM in the block
				MSD[i][t] += ( pow(RWL[j][0], 2) +  pow(RWL[j][1], 2) +  pow(RWL[j][2], 2) ) / L;
                }
            }
        }
	ofstream file1;
	file1.open("./data/latticeRW.dat");

	// Computing global averages and errors at each time step
	double avg[T] = {0};
	double avg2[T] = {0};
	double err[T] = {0};
	for(int t = 0; t < T; t++){
		for(int i = 0; i < N; i++){
			avg[t] += sqrt(MSD[i][t]) / N;
			avg2[t] += MSD[i][t] / N;
		}
		err[t] = sqrt(avg2[t] - avg[t] * avg[t]) / sqrt(N-1);
		file1 << t << " " << avg[t] << " " << err[t] << endl;
	}
	file1.close();


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
