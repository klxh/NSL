#include <math.h>
#include "population.h"
#include "specimen.h"
#include "random.h"

population::population(std::vector<specimen> spec){
    m_pop = spec;
}

population::population() {}

population::~population() {}

void population::populate(specimen spec){
    m_pop.push_back(spec);
}

void population::palingenesis(int Np, vertices cities, Random& rnd){
    m_pop.clear(); // Global deluge
    m_prob_selection.clear();

    // A new beginning
    for(int i=0; i<Np; i++){
        specimen c;
        c.init(cities, rnd);
        c.check(cities);
        m_pop.push_back(c);
    }
    return;
}

void population::clear() { 
    m_pop.clear(); 
    m_prob_selection.clear();
    return;
}

size_t population::size(){
    return m_pop.size();
}

specimen population::spec(int elem){
    if(elem < 0 || (size_t)elem >= m_pop.size()) {
        std::cerr << 
        "population::spec ERROR MESSAGE.\nTrying to access element out of bounds. Lower bound = 0. Upper bound = " 
        << m_pop.size() - 1 << ". EXIT PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }
    return m_pop[elem];
}

// Ranking the population by fitness
void population::rank(vertices cities, const std::string& norm){

    // Definition of a comparator function
    auto comparator = [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return a.second < b.second; // Sort based on the second column
    };

    std::vector<std::pair<int, double>> sf; // First entry: specimen index; Second entry: fitness level
    m_prob_selection.clear();
    m_prob_elimination.clear();
    double L=0; // Total fitness of the population
    for(size_t i=0; i<m_pop.size(); i++){
        m_prob_selection.push_back(m_pop[i].fitness(cities,norm));
        L+=m_prob_selection[i];
        sf.push_back(std::make_pair(i, m_prob_selection[i]));
    }

    // Constructing probabilities of elimination
    for(int i=0; i<m_prob_selection.size(); i++){ m_prob_elimination.push_back(m_prob_selection[i]/L); }
    // Constructing probabilities of selection
    for(int i=0; i<m_prob_selection.size(); i++) 
    { m_prob_selection[i] = (1 - m_prob_selection[i]/L) / (m_prob_selection.size() -1); }

    std::sort(sf.begin(), sf.end(), comparator);

    // Ordering population
    std::vector<specimen> tmp;
    std::vector<double> tmp_prob;
    for(size_t i=0; i<m_pop.size(); i++){
        tmp.push_back(m_pop[sf[i].first]); 
        tmp_prob.push_back(m_prob_selection[sf[i].first]);
    }

    m_pop = tmp;
    m_prob_selection = tmp_prob;
/*
    for(int i=0; i<m_pop.size(); i++){
        std::cout << "Individual " << i+1 << "\t Fitness = " << m_pop[i].fitness() 
        << "\t Probability of selection = " << m_prob_selection[i] << std::endl;
    }
*/
    return;
}

specimen population::fittest() {
    return m_pop[0];
}

double population::top_half_avg(vertices cities, const std::string& norm) {
    double avg=0;
    double top_half = (size_t)(double)m_pop.size()/2;
    for(size_t i=0; i<top_half; i++) {
        avg += m_pop[i].fitness(cities, norm) / top_half;
    }
    return avg;
}

specimen population::selection(Random& rnd) {
/*    
    double p=m_prob_selection[0];
    int selected = 0;
    double r = (double)random()/RAND_MAX;
    while(r > p){
        selected++;
        p += m_prob_selection[selected];
    }
*/
/**/
    double r = rnd.Rannyu();
    int selected = (int) m_pop.size() * pow(r, 5);
/**/
    return m_pop[selected];
}

void population::elimination(specimen young, Random& rnd){
    double p=m_prob_elimination[0];
    int eliminated = 0;
    double r = rnd.Rannyu();
    while(r > p){
        eliminated++;
        p += m_prob_elimination[eliminated];
    } 
    m_pop[eliminated] = young;
    return;
}
