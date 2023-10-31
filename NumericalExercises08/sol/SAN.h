#ifndef SAN_H
#define SAN_H

// Random numbers
#include "random.h"
Random rnd;
int seed[4];

// Variables

// Constants
const double hbar = 1.;
const double mass = 1.;

double w=5.0; // Transition width for sampling of trial function

const int M = 200000; // Number of MC steps (for each evaluation of the energy)
const int N = 100; // Number of blocks
const int L = M/N; 

// Trial parameters (starting, updated, transition width)
double m0=2.0, m, dm=.05, s0=2.0, s, ds=.05;
double x0=1, x; // Starting position and updated position

// Temperature
const int nA = 10000; // Number of annealing cycles
double a=1.0; // Temperature formula parameters
double b=3.75;
double T[nA]; // Temperature for each cycle
double beta; // 1/T

// Energy
// VMC energy: final average, block averages, final error
double E, avg_E[N], eps;
// Energy entering a simulated annealing cycle, exiting a cycle, and energy difference
double Energy_old, Energy_new, dE;

double alpha; // MCMC acceptance

double VMC_AR; // Variational MC acceptance rate

// I/O
FILE *fp;
FILE *gp;
char name[8*sizeof(long int)];

// Functions
void START(void);
void END(void);
void WRITE(int n);
double ENERGY(int n);

// Mathematical functions
double e(double x, double m, double s);
double psi(double x, double m, double s);
double p(double x, double m, double s);
double V(double x);
double f(double x, double m, double s);
double g(double x, double m, double s);


#endif
