#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include "specimen.h"
#include "population.h"
#include "tsp.h"
#include "vertices.h"
#include "random.h"

using namespace std;

/***************************************************************************************************
*****************************           MAIN          **********************************************
***************************************************************************************************/

int main(){

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

    // Reading data on positions of the cities
    ifstream capitals;
    capitals.open("./data/pos_US_caps.dat");
    if (!capitals.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return 1;
    }

    vector<double> x, y;
    vector<int> tags;
    double x_val, y_val;
    int tag = 1;
    while(capitals >> x_val >> y_val){
        x.push_back(x_val);
        y.push_back(y_val);
        tags.push_back(tag);
        tag++;
    }

    capitals.close();

    // Initialization of data structure with information on cities
    vertices cities(tags, x, y);

    fp1 = fopen("./data/fitness.dat", "w");
    fp2 = fopen("./data/route_in.dat", "w");
    fp3 = fopen("./data/route_out.dat", "w");

    // Initialization of population with individuals randomly built
    P.palingenesis(Np, cities, rnd);
    P.rank(cities, norm); // Ranking of the popuolation based on fitness

    fprintf(fp1, "%4d\t%.6f\t%.6f\n", 0, P.fittest().fitness(cities, norm), P.top_half_avg(cities,norm));
    for(size_t i=0; i<cities.size(); i++){
        fprintf(
            fp2, "%.6f\t%.6f\t%3d\n", 
            cities.get_x()[P.fittest().nucleotide(i)-1],
            cities.get_y()[P.fittest().nucleotide(i)-1],
            cities.get_vertex(P.fittest().nucleotide(i)-1)
        );
    }

    // EVOLUTION
    specimen parent1, parent2;
    std::vector<specimen> children;

    for(int n=0; n<gen; n++){
        for(size_t i=0; i<P.size()/2; i++){

            parent1 = P.selection(rnd);
            parent1.inversion(p_i, rnd);
            parent1.translation(p_t, rnd);
            parent1.switching(p_s, rnd);
            parent1.check(cities);

            parent2 = P.selection(rnd);
            parent2.inversion(p_i, rnd);
            parent2.translation(p_t, rnd);
            parent2.switching(p_s, rnd);
            parent2.check(cities);

            children = breed(parent1, parent2, p_c);
            children[0].check(cities);
            children[1].check(cities);

            P.elimination(children[0], rnd);
            P.elimination(children[1], rnd);
        }
        P.rank(cities, norm);
        fprintf(fp1, "%4d\t%.6f\t%.6f\n", n+1, P.fittest().fitness(cities, norm), P.top_half_avg(cities,norm));
    }
    for(size_t i=0; i<cities.size(); i++){
        fprintf(
            fp3, "%.6f\t%.6f\t%3d\n", 
            cities.get_x()[P.fittest().nucleotide(i)-1],
            cities.get_y()[P.fittest().nucleotide(i)-1],
            cities.get_vertex(P.fittest().nucleotide(i)-1)
        );
    }

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

    rnd.SaveSeed();

    return 0;
}

/***************************************************************************************************
*****************************        FUNCTIONS        **********************************************
***************************************************************************************************/

// Genetic crossover
std::vector<specimen> breed(specimen parent1, specimen parent2, double p_c){

    std::vector<specimen> children;

    double r = rnd.Rannyu();
    std::vector<int> child1_dna(parent1.size()), child2_dna(parent2.size());
    if(r <= p_c) {
        int C1 = 1 + rnd.Integer(parent1.size() - 1); // Cross spot 1
        int C2; // Cross spot 2
        do{ C2 = 1 + rnd.Integer(parent1.size() - 1); }
        while(C2 == C1);
        if(C1 > C2) { int a = C2; C2 = C1; C1 = a; } // To make C1 the smallest of the two

        std::vector<int> unavailables1, unavailables2;
        for(int i=0; i<C1; i++){
            child1_dna[i] = parent1.nucleotide(i); // Copying the preserved portion of parent DNA
            unavailables1.push_back(parent1.nucleotide(i)); // These elements can not be taken from other parent
            child2_dna[i] = parent2.nucleotide(i);
            unavailables2.push_back(parent2.nucleotide(i));
        }
        for(size_t i=C2+1; i<parent1.size(); i++){
            child1_dna[i] = parent1.nucleotide(i);
            unavailables1.push_back(parent1.nucleotide(i));
            child2_dna[i] = parent2.nucleotide(i);
            unavailables2.push_back(parent2.nucleotide(i));
        }

        std::vector<int> buffer1, buffer2;
        // Elements from other parent DNA for crossover
        for(size_t i=C1; i<parent1.size(); i++){
            buffer1.push_back(parent2.nucleotide(i));
            buffer2.push_back(parent1.nucleotide(i));
        }
        for(int i=0; i<C1; i++){
            buffer1.push_back(parent2.nucleotide(i));
            buffer2.push_back(parent1.nucleotide(i));
        }

        // Deleting candidates for crossover already fixed in the new DNA
        buffer1.erase(std::remove_if(buffer1.begin(), buffer1.end(), [&](int v) {
            return std::find(unavailables1.begin(), unavailables1.end(), v) != unavailables1.end();
            }), buffer1.end());

        buffer2.erase(std::remove_if(buffer2.begin(), buffer2.end(), [&](int v) {
            return std::find(unavailables2.begin(), unavailables2.end(), v) != unavailables2.end();
            }), buffer2.end());

        // Crossover of the legal genes
        for(int i=C1; i<=C2; i++){
            child1_dna[i] = buffer1[i-C1];
            child2_dna[i] = buffer2[i-C1];
        }
        specimen child1(child1_dna);
        specimen child2(child2_dna);
        children.push_back(child1);
        children.push_back(child2);
    }
    else { children.push_back(parent1); children.push_back(parent2); }
    return children;
}
