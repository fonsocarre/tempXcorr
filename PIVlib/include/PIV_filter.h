#pragma once
#include <vector>
#include <cmath>
#include <math.h>
#define _USE_MATH_DEFINES
#include "fftw3.h"

class PIV_filter
{

public:
    double cutLength;

    virtual void initialise(double cutLength) =0;
    virtual void filter(std::vector<double> in, 
                std::vector<double>& out,
                int nRows, int nCols, double dx,
                int padding) =0;
};
