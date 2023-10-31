#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include "SAN.h"

using namespace std;

int main (int argc, char *argv[]){

    START();

    // SIMULATED ANNEALING 

    for(int n=0; n<nA; n++){
        // Proposing new parameters
        m = m0 + dm * (rnd.Rannyu() - 0.5);
        s = s0 + ds * (rnd.Rannyu() - 0.5);
        
        // Computing new energy with proposed parameters
        Energy_new = ENERGY(n);
        dE = Energy_new - Energy_old;
        beta = 1./T[n];
        // Metropolis test for the proposed move
        alpha = 1.;
        if(exp(-beta*dE) <= alpha) alpha = exp(-beta*dE);
        // If the move is accepted, update parameters
        if (rnd.Rannyu() <= alpha) { m0 = m; s0 = s; Energy_old = Energy_new; }

        // Print the results
        WRITE(n);
    }

    END();

    return 0;
}

/***************************************************************************************************
 ***************************************************************************************************
 ***************************************************************************************************/

void START(){

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


    fp = fopen("./data/SAN.dat", "a");
    if (fp == NULL) {
        printf("Can't open ./data/SAN.dat\n");
        exit(EXIT_FAILURE);
    }

    // Energy initialization
    m = m0, s = s0;
    Energy_old = ENERGY(-1);
    fprintf(fp, "%4d\t%.8f\t%.8f\t%.8f\t%.8f\n", 0, m, s, E, eps);

    // Temperature initialization
    for(int i=0; i<nA; i++){
        T[i] = 1/(a+b*log(i+1));
    }

    return;
}

void WRITE(int n){
    printf("Annealing cycle # %3.3d\nTemperature = %.8f\n", n+1, T[n]);
    printf("VMC acceptance rate = %.4f\n****************************\n", VMC_AR);
    fprintf(fp, "%4d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", n+1, m0, s0, E, eps, T[n]);

    return;
}

double ENERGY(int n){
    
    if(n==nA-1 || n==-1){
        sprintf(name, "./data/hst%4.4d.dat", n+1);
        gp = fopen(name, "w");
    }

    // Clearing old variables
    VMC_AR = 0; E = 0; eps = 0;
    for(int i=0; i<N; i++){ avg_E[i]=0; }

    // MCMC computation of average energy with Variational Monte Carlo
    for(int i=0; i<N; i++){
        for(int j=0; j<L; j++){
            // Metropolis move
            x = x0 + w * (rnd.Rannyu()-.5);
            // Test for the proposed move
            alpha = 1.0;
            if(alpha > p(x,m,s)/p(x0,m,s)) alpha = p(x,m,s)/p(x0,m,s);
            if(rnd.Rannyu() <= alpha){
                x0 = x;
                VMC_AR += 1./M;
            }
            if(n==nA-1 || n==-1){
                fprintf(gp, "%.8f\n", x0);
                // Updating block average
            }
            avg_E[i] += g(x0,m,s)/L;
        }
        E += avg_E[i]/N;
        eps += pow(avg_E[i], 2)/N;
    } 
    eps = sqrt(eps - E*E)/sqrt(N);

    if(n==nA-1 || n==-1) {fclose(gp);}

    return E;
}

void END(){
    fclose(fp);
    rnd.SaveSeed();
    return;
}

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
