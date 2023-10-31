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

// Defining integrand function
double g(double x){
    return .5*M_PI*cos(.5*M_PI*x);
    }

// Defining sampling probability distribution for importance sampling
double q(double x){
    return 1.5 * (1 - x * x) ;
    }

/* Taylor expansion
double q(double x){
    return (1 - M_PI*M_PI/8 * x * x) / (1 - M_PI*M_PI / 24);
    }
*/

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

/*
   for(int i=0; i<20; i++){
      cout << rnd.Rannyu() << endl;
   }

   rnd.SaveSeed();
*/

    // MC Integration

    // 1. Sampling uniform distribution

    int M = 1000000; // Number of throws
    int N = 100; // Number of blocks
    int L = M/N; // Number of throws per block

    double I[N] = {0}; // Block estimate of integral
    double I2[N] = {0}; // I^2

    ofstream file1;
    file1.open("./data/uniform.dat");

    for(int i = 0; i < N; i++){
        for(int j = 0; j < L; j++){
            I[i] += g(rnd.Rannyu()) / L; // Computing block value of integral
            }
        I2[i] = I[i] * I[i];
        }

    // Progressive average values
    double avg_prog[N] = {0}; // Progressive average of I
    double avg2_prog[N] = {0}; // Progressive average of I^2
    double err_prog[N] = {0}; // Progressive std dev of the mean
    double th_value = 1.; // Theoretical value of the integral

    // Computing averages with data blocking method
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){
            avg_prog[i] += I[j] / (i+1);
            avg2_prog[i] += I2[j] / (i+1);
            }
        err_prog[i] = sqrt(avg2_prog[i] - avg_prog[i] * avg_prog[i] ) / sqrt(i+1);

        // Saving data to file
        file1 << i+1 << " " << avg_prog[i] - th_value << " " << err_prog[i] << std::endl;
        }

    file1.close();

    // Importance sampling

    // Clearing old variables
    for(int i = 0; i < N; i++){
        I[i] = 0;
        I2[i] = 0;
        avg_prog[i] = 0;
        avg2_prog[i] = 0;
        err_prog[i] = 0;
        }

    ofstream file2;
    file2.open("./data/importance.dat");

    int count = 0; // Counts number of samplings (just for information purpose)

    double q_max = q(0); // Maximum of function q(x)
    for(int i = 0; i < N; i++){
        for(int j = 0; j < L; j++){
            double x, r;
            // Sampling function q(x) with rejection technique
            do { x = rnd.Rannyu(); r = rnd.Rannyu(); count++; } while(r > q(x)/q_max);
            I[i] += g(x)/q(x)/L;
            }
        I2[i] = I[i]*I[i];
        }

    // Computing averages with data blocking method
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){
            avg_prog[i] += I[j] / (i+1);
            avg2_prog[i] += I2[j] / (i+1);
            }
        err_prog[i] = sqrt(avg2_prog[i] - avg_prog[i] * avg_prog[i]) / sqrt(i+1);

        // Saving data to file
        file2 << i+1 << " " << avg_prog[i] - th_value << " " << err_prog[i] << std::endl;
        }

    file2.close();

    cout << "# PRNG calls (rejected + accepted) = " << 2*count << endl;

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
