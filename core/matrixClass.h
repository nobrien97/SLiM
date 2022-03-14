#include <vector>

template <class T>
class matrix
{
    std::vector<std::vector<T>> m;

 public:
    matrix(unsigned int x, unsigned int y) {
        m.resize(x, std::vector<T>(y, false));
        _y = y;
        _x = x;
    }
    // If we're given an initialiser value, initialise to that value
    matrix(unsigned int x, unsigned int y, T init) {
        m.resize(x, std::vector<T>(y, false));
        _y = y;
        _x = x;
        for (unsigned int i = 0; i < x; ++i)
        {
            for (unsigned int j = 0; j < y; ++j)
            {
                m[i][j] = init;
            }
        }
    }
    std::vector<T> flattenMatrix()
    {
        int x = nCols();
        int y = nRows();
        std::vector<T> result;
        result.reserve(x * y * sizeof(T));

        for (int i = 0; i < x; ++i)
        {
            for (int j = 0; j < y; ++j)
            {
                result.emplace_back(this[i][j]);
            }
        }
        return result; 
    }

    unsigned int nRows() {return _y; }
    unsigned int nCols() {return _x; }


    class matrix_row
    {
        std::vector<T>& row;
     public:
        matrix_row(std::vector<T>& r) : row(r) {}
        T& operator[](unsigned int y) { return row.at(y);}
    };
    matrix_row& operator[](unsigned int x) { return matrix_row(m.at(x));}

 private:
    unsigned int _y;
    unsigned int _x;
};