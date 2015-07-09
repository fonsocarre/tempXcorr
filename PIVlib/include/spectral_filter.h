#pragma once
#include "PIV_filter.h"
#include "PIV_windows.h"
#include <iostream>
#include <cassert>

class spectral_filter: public PIV_filter
{

public:
    void initialise(double cutLength);
    void filter(std::vector<double> in,
                std::vector<double>& out,
                int nRows, int nCols, 
                double dx, int padding);
};
