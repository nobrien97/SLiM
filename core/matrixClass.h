#include <vector>

template <class T>
class matrix
{
    std::vector<std::vector<T>> m;

 public:
    matrix(unsigned int x, unsigned int y) {
        m.resize(x, std::vector<T>(y, false));
    }
    // If we're given an initialiser value, initialise to that value
    matrix(unsigned int x, unsigned int y, T init) {
        m.resize(x, std::vector<T>(y, false));
        for (int i = 0; i < x; ++i)
        {
            for (int j = 0; j < y; ++j)
            {
                m[i][j] = init;
            }
        }

    }

    class matrix_row
    {
        std::vector<T>& row;
     public:
        matrix_row(std::vector<T>& r) : row(r) {}
        bool& operator[](unsigned int y { return row.at(y);})
    };
    matrix_row& operator[](unsigned int x) { return matrix_row(m.at(x));}

};