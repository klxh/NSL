#ifndef __SPECIMEN__
#define __SPECIMEN__

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "vertices.h"

class specimen{
    public:
    specimen(std::vector<int> spec);
    specimen();
    ~specimen();

    void print();
    void check(vertices cities);
    void init(vertices cities);
    size_t size();
    int nucleotide(int spot);

    double fitness(vertices cities, const std::string& norm);
    double fitness(void) { return m_fitness; }

    // Mutations
    void inversion(double p_i);
    void translation(double p_t);
    void switching(double p_s);

    private:
    std::vector<int> m_spec;
    double m_fitness;

};

#endif
