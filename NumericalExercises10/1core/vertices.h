#ifndef __VERTICES__
#define __VERTICES__

#include <vector>
#include <string>
#include "random.h"

class vertices{
    private:
    std::vector<int> m_index;
    std::vector<double> m_x;
    std::vector<double> m_y;

    public:
    vertices(std::vector<int> index);
    vertices(std::vector<int> index, std::vector<double> x, std::vector<double> y);
    vertices();
    ~vertices();

    void init(int Nv, const std::string& geometry, Random& rnd);

    std::vector<int> get_vertices();
    int get_vertex(int elem);
    size_t size() { return m_index.size(); }

    std::vector<double> get_x();
    std::vector<double> get_y();

};

#endif
