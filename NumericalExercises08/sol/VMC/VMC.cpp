#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include "VMC.h"
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

    // Starting configuration
    double x0=1;
    double m=0.81618876;
    double s=0.59989305;
	
    double w = 5.0 ; // maximum width of a transition x1 -> x2
    
    int M = 2000000; // Number of throws
    int N = 100; // Number of blocks
    int L = M/N; // Number of throws per block

    double AR = 0; // Metropolis acceptance rate

    double avg[N] = {0};
    double avg2[N] = {0};
    double err[N] = {0};

    FILE *fp1, *fp2, *fp3;
    fp1 = fopen("./data/energy.dat", "a");
    if (fp1 == NULL) {
        printf("Can't open energy.dat\n");
        exit(EXIT_FAILURE);
    }

    fp2 = fopen("./data/hst.dat", "w");
    if (fp2 == NULL) {
        printf("Can't open hst.dat\n");
        exit(EXIT_FAILURE);
    }

    fp3 = fopen("./data/final_config.xyz", "w");
    if (fp3 == NULL) {
        printf("Can't open final_config.xyz\n");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j < L; j++){
            // Sampling new position from uniform transition rate with width w
            double x = x0 + w * (rnd.Rannyu() - .5);

            // Acceptance rate
            double alpha = 1.0;
            if (alpha > p(x,m,s) / p(x0,m,s)) alpha = p(x,m,s) / p(x0,m,s);

            // Acceptance test
            if (rnd.Rannyu() <= alpha){
                // Accepted
                x0 = x;
                AR += 1./M; // Update acceptance rate
                }
            // Write positions to file
            fprintf(fp2, "%.6f\n", x0);

            // Update block average
            avg[i] += g(x0,m,s) / L;
            }
        avg2[i] = avg[i] * avg[i];
        }
    
    // Write final configuration
    fprintf(fp3, "%.6f\n", x0);

    double prog_avg[N] = {0};
    double prog_avg2[N] = {0};
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){
            prog_avg[i] += avg[j] / (i+1);
            prog_avg2[i] += avg2[j] / (i+1);
            }
        err[i] = sqrt(prog_avg2[i] - prog_avg[i] * prog_avg[i]) / sqrt(i+1);
        // Write average values to file
        fprintf(fp1, "%3.0d\t%.6f\t%.6f\n", i+1, prog_avg[i], err[i]);
        }

    cout << "Acceptance rate  = " << AR << endl;

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

    rnd.SaveSeed();


    return 0;
}

/***************************************************************************************************
 ***************************************************************************************************
 ***************************************************************************************************/

// Mathematical functions

double e(double x, double m, double s){
    return exp(-.5*pow((x-m)/s, 2));
}

// Trial function
double psi(double x, double m, double s){
    return e(x,m,s) + e(-x,m,s);
}

// Probability density (not normalized)
double p(double x, double m, double s){
    return psi(x,m,s)*psi(x,m,s);
}

// Potential energy
double V(double x){
    return pow(x, 4) - 2.5*pow(x, 2);
}

// Intermediate function for computing the kinetic energy
double f(double x, double m, double s){
    return -pow(hbar, 2) / (2*mass*pow(s,4)) * (pow(x-m, 2) - pow(s, 2));
}

// (Hpsi)/psi
double g(double x, double m, double s){
    return V(x) + (f(x,m,s)*e(x,m,s) + f(-x,m,s)*e(-x,m,s))/(e(x,m,s) + e(-x,m,s));
}
