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

using namespace std;

/***************************************************************************************************
*****************************           MAIN          **********************************************
***************************************************************************************************/

int main(){

    //srandom(time(NULL));
    srandom(seed);

    fp1 = fopen("./data/fitness.dat", "w");
    fp2 = fopen("./data/route_in.dat", "w");
    fp3 = fopen("./data/route_out.dat", "w");

    // Initialization of data structure with information on cities
    cities.init(Nc, geometry);

    // Initialization of population with individuals randomly built
    P.palingenesis(Np, cities);
    P.rank(cities, norm); // Ranking of the popuolation based on fitness

    fprintf(fp1, "%4d\t%.6f\t%.6f\n", 0, P.fittest().fitness(cities, norm), P.top_half_avg(cities,norm));
    for(int i=0; i<Nc; i++){
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

            parent1 = P.selection();
            parent1.inversion(p_i);
            parent1.translation(p_t);
            parent1.switching(p_s);
            parent1.check(cities);

            parent2 = P.selection();
            parent2.inversion(p_i);
            parent2.translation(p_t);
            parent2.switching(p_s);
            parent2.check(cities);

            children = breed(parent1, parent2, p_c);
            children[0].check(cities);
            children[1].check(cities);

            P.elimination(children[0]);
            P.elimination(children[1]);
        }
        P.rank(cities, norm);
        fprintf(fp1, "%4d\t%.6f\t%.6f\n", n+1, P.fittest().fitness(cities, norm), P.top_half_avg(cities,norm));
    }
    for(int i=0; i<Nc; i++){
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
    
    return 0;
}

/***************************************************************************************************
*****************************        FUNCTIONS        **********************************************
***************************************************************************************************/

// Genetic crossover
std::vector<specimen> breed(specimen parent1, specimen parent2, double p_c){

    std::vector<specimen> children;

    double r = (double) random()/RAND_MAX;
    std::vector<int> child1_dna(parent1.size()), child2_dna(parent2.size());
    if(r <= p_c) {
        int C1 = 1 + random() % (parent1.size() - 1); // Cross spot 1
        int C2; // Cross spot 2
        do{ C2 = 1 + random() % (parent1.size() - 1); }
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
