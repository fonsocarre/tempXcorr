#pragma once

#include <vector>
#include <stdlib.h>
#include <cmath>

static inline int ij(int& i, int& j, int& nRows, int& nCols)
{
    return i*nCols + j;
}   

static inline int cij(int i, int j, int& nRows, int& nCols)
{
    return i*nCols + j;
}   
namespace PIV
{
    void xxcorr(const std::vector<double>& mat0,
               const std::vector<double>& mat1,
               const int wNx,
               const int wNy,
               int nRows,
               int nCols,
               std::vector<double>& xdisps,
               std::vector<double>& ydisps,
               std::vector<double>& xcorrvals);
}

template <unsigned int exponent>
inline double intpow(double base)
{
    return intpow<exponent-1>(base) * base;
}

template <>
inline double intpow<0>(double base)
{
    return 1;
}
