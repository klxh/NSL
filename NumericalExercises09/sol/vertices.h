#ifndef __VERTICES__
#define __VERTICES__

#include <vector>
#include <string>

class vertices{
    private:
    std::vector<int> m_index;
    std::vector<double> m_x;
    std::vector<double> m_y;

    public:
    vertices(std::vector<int> index);
    vertices();
    ~vertices();

    void init(int Nv, const std::string& geometry);

    std::vector<int> get_vertices();
    int get_vertex(int elem);

    std::vector<double> get_x();
    std::vector<double> get_y();

};

#endif
