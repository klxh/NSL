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

#define a0 1.0 // Bohr radius

using namespace std;

// Definition of |psi_100|^2
double psi2(double x, double y, double z){
    return exp(-2 * sqrt(x*x + y*y + z*z) / a0) / (M_PI * pow(a0, 3));
    }

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

    // Psi 100 with gaussian transition rate

    // Starting configuration
    double x0, y0, z0;
    bool restart = 1;
    if(restart == 0){
        x0 = 200.0 * a0;
        y0 = 200.0 * a0;
        z0 = 200.0 * a0;
    }
    else{
        ifstream infile("./data/final_config.xyz");
        // Check if the file was opened successfully
        if (!infile.is_open()) {
            cout << "Error opening file 'final_config.xyz'" << endl;
            return 1;
        }
        infile >> x0 >> y0 >> z0;
        infile.close();        
    }

    double w = 0.85 * a0; // maximum width of a transition x1 -> x2
    
    int M = 1000000; // Number of throws
    int N = 100; // Number of blocks
    int L = M/N; // Number of throws per block

    double AR = 0; // Metropolis acceptance rate

    double avg[N] = {0};
    double avg2[N] = {0};
    double err[N] = {0};

    FILE *fp1, *fp2, *fp3;
    fp1 = fopen("./data/gaussian100.dat", "a");
    if (fp1 == NULL) {
        printf("Can't open gaussian100.dat\n");
        exit(EXIT_FAILURE);
    }

    fp2 = fopen("./data/gaussian100.xyz", "w");
    if (fp2 == NULL) {
        printf("Can't open gaussian100.xyz\n");
        exit(EXIT_FAILURE);
    }

    fp3 = fopen("./data/final_config.xyz", "w");
    if (fp3 == NULL) {
        printf("Can't open final_config.xyz\n");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j < L; j++){
            // Sampling new position from gaussian transition rate with std deviation w
            double x = x0 + w * rnd.Gauss(0,w) ;
            double y = y0 + w * rnd.Gauss(0,w) ;
            double z = z0 + w * rnd.Gauss(0,w) ;

            // Acceptance rate
            double alpha = 1.0;
            if (alpha > psi2(x,y,z) / psi2(x0,y0,z0)) alpha = psi2(x,y,z) / psi2(x0,y0,z0);

            // Acceptance test
            double r = rnd.Rannyu();
            if (r <= alpha){
                // Accepted
                x0 = x;
                y0 = y;
                z0 = z;
                AR += 1./M; // Update acceptance rate
                }
            // Write positions to file
            fprintf(fp2, "%.4f\t%.4f\t%.4f\n", x0, y0, z0);

            // Update block average
            avg[i] += sqrt(x0*x0 + y0*y0 + z0*z0) / L;
            }
        avg2[i] = avg[i] * avg[i];
        }
    
    // Write final configuration
    fprintf(fp3, "%.4f\t%.4f\t%.4f\n", x0, y0, z0);

    double prog_avg[N] = {0};
    double prog_avg2[N] = {0};
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){
            prog_avg[i] += avg[j] / (i+1);
            prog_avg2[i] += avg2[j] / (i+1);
            }
        err[i] = sqrt(prog_avg2[i] - prog_avg[i] * prog_avg[i]) / sqrt(i+1);
        // Write average values to file
        fprintf(fp1, "%3.0d\t%.4f\t%.4f\n", i+1, prog_avg[i], err[i]);
        }

    cout << "Acceptance rate  = " << AR << endl;

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

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
