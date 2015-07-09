#pragma once

#include <vector>
#include <cmath>
#include <math.h>
#include "fftw3.h"
#define _USE_MATH_DEFINES


void getHanningWindow(std::vector<double> initial,
                        std::vector<double>& modified,
                       int nRows, int nCols);

fftw_complex* spectralHanningWindow(int nRows, int nCols);


fftw_complex* zeroPadding(fftw_complex* input, int nRows, int nCols,
                                int padding);
fftw_complex* removePadding(fftw_complex* input, int nRows, int nCols,
                                int padding);
