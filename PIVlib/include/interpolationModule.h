//
//  interpolationModule.h
//  Stitch
//
//  Created by Alfonso del Carre on 08/06/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//


#pragma once

extern "C"
{
    void interpolate(double* ixinL,
                     double* iyinL,
                     double* ivarin1L,
                     double* ivarin2L,
                     double* ixinR,
                     double* iyinR,
                     double* ivarin1R,
                     double* ivarin2R,
                     int* nRowsin,
                     int* nColsin,
                     double* ixout,
                     double* iyout,
                     double* varout1,
                     double* varout2,
                     int* nRowsout,
                     int* nColsout,
                     double* xdisp,
                     double* ydisp);
}