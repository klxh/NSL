#include "vertices.h"
#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

vertices::vertices(std::vector<int> index) {
    m_index = index;
}

vertices::vertices(std::vector<int> index, std::vector<double> x, std::vector<double> y) {
    m_index = index;
    m_x = x;
    m_y = y;
}

vertices::vertices() {}

vertices::~vertices() {}


// Initialization

void vertices::init(int Nv, const std::string& geometry, Random& rnd) {
    m_index.clear();
    m_x.clear();
    m_y.clear();

    for(int i=0; i<Nv; i++){
        m_index.push_back(i+1);
    }
    
    if(geometry=="CIRCLE"){
        double x,y,dr;
        // Sampling position of vertex i with uniform probability on unit circle
        for(int i=0; i<Nv; i++){
            do{
                x = 2*(rnd.Rannyu()-.5);
                y = 2*(rnd.Rannyu()-.5);
                dr = sqrt(x*x + y*y);
            } while(dr > 1);
            m_x.push_back(x/dr);
            m_y.push_back(y/dr);
        }
    }
    else if(geometry=="SQUARE") {
        // Sampling position of vertex i with uniform probability in [0,1]^2
        for(int i=0; i<Nv; i++){
            m_x.push_back(rnd.Rannyu());
            m_y.push_back(rnd.Rannyu());
        }
    }
    else {
        std::cerr << "Geometry not defined. Accepted arguments are SQUARE and CIRCLE.\n EXIT PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }
}

// Get elements
std::vector<int> vertices::get_vertices() {
    return m_index;
}

int vertices::get_vertex(int elem) {
    if (elem < 0 || m_index.size() <= (size_t)elem) {
        std::cerr << 
        "vertices::get_vertex ERROR MESSAGE.\nTrying to access element out of bounds. Lower bound = 0. Upper bound = " 
        << m_index.size()-1 << ". EXIT PROGRAM" << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program 
    }
    return m_index[elem];
}

std::vector<double> vertices::get_x() {
    return m_x;
}

std::vector<double> vertices::get_y() {
    return m_y;
}
