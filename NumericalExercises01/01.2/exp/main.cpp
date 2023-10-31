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


	// How to call PRNG: 1. Uniform, rnd.Rannyu(); 2. Exponential, rnd.Exp(lambda); 3. Cauchy, rnd.Cauchy(mu, gamma).

	// Exponential distribution with lambda = 1
    double lmbd = 1.;
	int samples = 10000;

	ofstream myfile1, myfile2, myfile3, myfile4;
	myfile1.open("./data/exponentialN1.dat");

    // Writing samples from U[0,1] to file
	for(int i = 0; i < samples; i++){
		double r = rnd.Exp(lmbd);
		myfile1 << r << endl;
	}

	myfile1.close();
	
	// Now I will sample (X1+X2)/2, where X1 and X2 are iid ~ U[0,1]
	myfile2.open("./data/exponentialN2.dat");

	for(int i = 0; i < samples; i++){
        double r1 = rnd.Exp(lmbd);
        double r2 = rnd.Exp(lmbd);
        myfile2 << (r1 + r2)/2 << endl;
	}


    // Now I will do the same for (X1 + ... + X10)/10
	myfile3.open("./data/exponentialN10.dat");

	for(int i = 0; i < samples; i++){
        double r = 0;
		for(int j = 0; j < 10; j++) { r += rnd.Exp(lmbd)/10; }
        myfile3 << r << endl;
	}

	myfile3.close();

    // Now I will do the same for (X1 + ... + X100)/100
	myfile4.open("./data/exponentialN100.dat");

	for(int i = 0; i < samples; i++){
        double r = 0;
		for(int j = 0; j < 100; j++) { r += rnd.Exp(lmbd)/100; }
        myfile4 << r << endl;
	}

	myfile4.close();

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
