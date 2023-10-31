#include <stdlib.h>
#include <math.h>
#include "specimen.h"

using namespace std;

specimen::specimen(std::vector<int> spec) {
    m_spec = spec;
}

specimen::specimen() {}

specimen::~specimen() {}

void specimen::print(){
    for (int i : m_spec) printf("%2d ", i);
    printf("\n");
}

// Check legal solutions
void specimen::check(vertices cities){
    std::vector<int> sorted;
    sorted = m_spec;
    std::sort(sorted.begin(), sorted.end());
    // When sorted in increasing order the listing should match the reference list
    if(sorted != cities.get_vertices()){
        std::cerr << "Ill-defined element\nEXITING PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program
    }
    // The first element should always be 1 in order to avoid degenerate solutions
    if(m_spec[0] != 1){
        std::cerr << "First vertex is not 1\nEXITING PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }
    return;
}

// Initialization (random route)
void specimen::init(vertices cities){
    std::vector<int> remaining = cities.get_vertices();
    m_spec.push_back(1); // Starting vertex of the route always 1
    remaining.erase(remaining.begin()); // To complete the route we have one less vertex
    for(size_t i=1; i<cities.get_vertices().size(); i++){
        // Next vertex chosen among remaining vertices
        double r = (double) random()/RAND_MAX;
        int c = (int) (r * remaining.size());
        m_spec.push_back(remaining[c]);
        // In the next iteration we can not choose the element just picked
        remaining.erase(remaining.begin()+c); 
    }
}

size_t specimen::size(){
    return m_spec.size();
}

// Access elements
int specimen::nucleotide(int spot){
    if(spot < 0 || spot >= m_spec.size()){
        std::cerr << "Trying to access element out of boundaries.\nEXITING PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }
    return m_spec[spot];
}

// Computes "fitness level" (distance)
double specimen::fitness(vertices cities, const std::string& norm){
    double F=0;
    double dx, dy;
    if(norm=="L1") {
        for(size_t i=0; i<m_spec.size() - 1; i++){
           dx = cities.get_x()[m_spec[i+1]-1] - cities.get_x()[m_spec[i]-1];
           dy = cities.get_y()[m_spec[i+1]-1] - cities.get_y()[m_spec[i]-1];
           F += sqrt(dx*dx + dy*dy);
        }
        dx = cities.get_x()[m_spec[0]-1] - cities.get_x()[m_spec[m_spec.size()-1]-1];
        dy = cities.get_y()[m_spec[0]-1] - cities.get_y()[m_spec[m_spec.size()-1]-1];
        F += sqrt(dx*dx + dy*dy);
    }
    else if (norm=="L2") {
        for(size_t i=0; i<m_spec.size() - 1; i++){
           dx = cities.get_x()[m_spec[i+1]-1] - cities.get_x()[m_spec[i]-1];
           dy = cities.get_y()[m_spec[i+1]-1] - cities.get_y()[m_spec[i]-1];
           F += dx*dx + dy*dy;
        }
        dx = cities.get_x()[m_spec[0]-1] - cities.get_x()[m_spec[m_spec.size()-1]-1];
        dy = cities.get_y()[m_spec[0]-1] - cities.get_y()[m_spec[m_spec.size()-1]-1];
        F += dx*dx + dy*dy;
    }
    else {
        std::cerr << "Norm not defined. Accepted arguments are L1 and L2.\n EXIT PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }

    m_fitness = F;
    return F;
}

// MUTATIONS

void specimen::inversion(double p_i){

    // Checking if probability is well defined
    if(p_i < 0 || p_i > 1){
        std::cerr << "Input probability not acceptable.\n EXIT PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }

    double r = (double) random()/RAND_MAX;
    if(r <= p_i){
        // Inversion
        std::vector<int> tmp;
        int P1, P2; // Inversion points
        do { P1 = random() % m_spec.size(); }
        while (P1 == 0); // Position 0 is always left unchanged
        do { P2 = random() % m_spec.size(); }
        while (P2 == 0 || P2 == P1);

        if(P2 < P1) { int a = P2; P2 = P1; P1 = a; } // So that P1 is always the smallest

        for(size_t i=0; i<m_spec.size(); i++){
            if(i<P1 || i>P2) { tmp.push_back(m_spec[i]); }
            else if (P1 <= i && i <= P2) { tmp.push_back(m_spec[P1+P2-i]); }
        }
        m_spec = tmp;
    }
    return;
}

void specimen::translation(double p_t){

    // Checking if probability is well defined
    if(p_t < 0 || p_t > 1){
        std::cerr << "Input probability not acceptable.\n EXIT PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }

    double r = (double) random()/RAND_MAX;
    if(r <= p_t){
        // Translation
        std::vector<int> tmp;
        int P1, P2; // Borders of the sequence to be translated
        int T; // Translation constant

        P1 = 1 + random() % (m_spec.size() - 1);
        do { P2 = 1 + random() % (m_spec.size() - 1); }
        while (abs(P2-P1) == m_spec.size() - 2 || P2 == P1);

        if(P2 <= P1) { int a = P2; P2 = P1; P1 = a; } // So that P1 is always the smallest

        do { T = 1 + random()%(m_spec.size() - 1); }
        while (T >= P1 && T <= P2);

        int shift;
        int diff = P2-P1;
        std::vector<int> index;

        if(T>P2){
            shift=T-P1+1;
            int a = m_spec.size() - shift;
            int b = P2 + shift;
            if(b > m_spec.size()){
                for(int i=a; i <= P2; i++) {index.push_back(i);}
                for(int i=T+1; i<m_spec.size(); i++) {index.push_back(i);}
                for(int i=1; i<P1; i++) {index.push_back(i);}
                for(int i=P2+1; i<=T; i++) {index.push_back(i);}
                for(int i=P1; i<a; i++) {index.push_back(i);}
            }
            else{
                for(int i=1; i<P1; i++) {index.push_back(i);}
                for(int i=P2+1; i<=T; i++) {index.push_back(i);}
                for(int i=P1; i<=P2; i++) {index.push_back(i);}
                for(int i=T+1; i<m_spec.size(); i++) {index.push_back(i);}
            }
        }
        if(T<P1){
            shift=T-P1+1;
            int a = 1 - shift;
            int b = P2 + shift;
            for(int i=1; i <= T; i++) {index.push_back(i);}
            for(int i=P1; i<=P2; i++) {index.push_back(i);}
            for(int i=T+1; i<P1; i++) {index.push_back(i);}
            for(int i=P2+1; i<m_spec.size(); i++) {index.push_back(i);}
        }
        tmp.push_back(m_spec[0]);
        for(int i=0; i<index.size(); i++) { 
            tmp.push_back(index[i]+1);
        }
        m_spec = tmp;
    }
    return;
}

void specimen::switching(double p_s){

    // Checking if probability is well defined
    if(p_s < 0 || p_s > 1){
        std::cerr << "Input probability not acceptable.\n EXIT PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }

    double r = (double)random()/RAND_MAX;
    if(r <= p_s){
        int P1, P2; // Switching points
        P1 = 1 + random() % (m_spec.size() - 1);
        do { P2 = 1 + random() % (m_spec.size() - 1); }
        while (P2 == P1);

        std::vector<int> tmp = m_spec;
        tmp[P1] = m_spec[P2];
        tmp[P2] = m_spec[P1];
        m_spec = tmp;
    }
    return;
}
