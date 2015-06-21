#pragma once

#include <vector>

static inline int ij(int& i, int& j, int& nRows, int& nCols)
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
