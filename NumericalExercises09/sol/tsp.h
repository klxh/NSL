#ifndef __PARAMS__
#define __PARAMS__

#include <vector>

#include <vector>
#include <iostream>
#include "vertices.h"
#include "population.h"

FILE *fp1;
FILE *fp2;
FILE *fp3;

const std::string geometry = "CIRCLE"; // Possibilities: {"CIRCLE", "SQUARE"}
const std::string norm = "L1"; // Possibilities: {"L1", "L2"}
const int Nc = 34; // Number of cities
const int Np = 1000; // Individuals in the population
const int gen = 200; // Number of generations

// Mutation probabilities
const double p_i = 0.033; // Inversion
const double p_t = 0.033; // Translation
const double p_s = 0.033; // Switching

const double p_c = 0.80; // Crossover probability

const int seed = 33;

vertices cities;
population P;

// Breed function
std::vector<specimen> breed(specimen parent1, specimen parent2, double p_c);

#endif
