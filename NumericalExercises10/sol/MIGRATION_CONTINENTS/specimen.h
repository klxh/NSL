#ifndef __SPECIMEN__
#define __SPECIMEN__

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "vertices.h"
#include "random.h"

class specimen{
    public:
    specimen(std::vector<int> spec);
    specimen();
    ~specimen();

    void print();
    void check(vertices cities);
    void init(vertices cities, Random& rnd);
    size_t size();
    int nucleotide(int spot);

    std::vector<int> dna(void) { return m_spec; }

    double fitness(vertices cities, const std::string& norm);
    double fitness(void) { return m_fitness; }

    // Mutations
    void inversion(double p_i, Random& rnd);
    void translation(double p_t, Random& rnd);
    void switching(double p_s, Random& rnd);

    private:
    std::vector<int> m_spec;
    double m_fitness;

};

#endif
