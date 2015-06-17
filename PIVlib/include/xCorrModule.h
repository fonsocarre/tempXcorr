//
//  stitchutils.h
//  Stitch
//
//  Created by Alfonso del Carre on 03/06/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#pragma once

extern "C"
{
    void xcorr(const double* matL, const double* matR, int* nRows, int* nCols,
               int* overlap, double* xdisp, double* ydisp, bool* withInterpolation);

    void tempxcorr(double* mat0,
                   double* mat1,
                   int* xmin,
                   int* xmax,
                   int* ymin,
                   int* ymax,
                   double* offset,
                   int* expon,
                   int* nRows,
                   int* nCols,
                   int* xdisp,
                   int* ydisp,
                   double* corrVal);
    double calculatexcorr(double* mat0,
                          double* mat1,
                          int* ny,
                          int* nx,
                          int* xdisp,
                          int* expon,
                          double* offset,
                          int* jmin);
}
