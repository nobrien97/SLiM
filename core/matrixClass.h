#include <vector>

#define OUT_OF_RANGE 0x0A

template <class T>
class matrix
{
 public:
    matrix(size_t rows, size_t cols) : nRows(rows), nCols(cols), m(rows * cols) {}
    // If we're given an initialiser value, initialise to that value
    matrix(size_t rows, size_t cols, T init) : nRows(rows), nCols(cols), m(rows * cols) 
    {
        for (unsigned int i = 0; i < rows; ++i)
        {
            for (unsigned int j = 0; j < cols; ++j)
            {
                m[i * cols + j] = init;
            }
        }
    }

    std::vector<T> getVec() { return m; }

    T& operator()(size_t i, size_t j) 
    { 
        if ((i * nCols + j) >= m.size())
            throw OUT_OF_RANGE;
        return m[i * nCols + j]; 
    }

    T operator()(size_t i, size_t j) const 
    { 
        if ((i * nCols + j) >= m.size())
            throw OUT_OF_RANGE;
        return m[i * nCols + j]; 
    }

 private:
    unsigned int nCols;
    unsigned int nRows;
    std::vector<T> m;
};