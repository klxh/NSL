#ifndef __POPULATION__
#define __POPULATION__

#include <vector>
#include "specimen.h"
#include "vertices.h"
#include "random.h"

class population{
    private:
    std::vector<specimen> m_pop; // Elements of the population
    std::vector<double> m_prob_selection; // Selection probabilities
    std::vector<double> m_prob_elimination; // Elimination probabilities

    public:
    population(std::vector<specimen> spec);
    population();
    ~population();

    specimen spec(int elem); // Get element elem of the population

    void populate(specimen spec); // Adds element spec to the population
    void palingenesis(int Np, vertices cities, Random& rnd); // Clear and initialize new population
    void clear();
    size_t size();

    // Reorder the elements of the population on the basis of fitness level ([0]=fittest)
    void rank(vertices cities, const std::string& norm); 
    specimen fittest();
    double top_half_avg(vertices cities, const std::string& norm);

    // Selection operator
    specimen selection(Random& rnd);
    void elimination(specimen young, Random& rnd);


};

#endif
