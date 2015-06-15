//
//  PIV_Stitch.h
//  Stitch
//
//  Created by Alfonso del Carre on 01/06/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>

#include <fstream>
#include <algorithm>    // std::min_element, std::max_element
#include <thread>
#include <mutex>


#include "Set.h"
#include "Picture.h"
#include "Frame.h"
#include "getSettings.h"
#include "xCorrModule.h"
#include "interpolationModule.h"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
typedef CGAL::Interpolation_traits_2<K>               Traits;
typedef K::FT                                         Coord_type;
typedef K::Point_2                                    Point;

// for parallelisation and avoid HDF5 problems
static std::mutex m;

namespace PIV
{
    void parallelStitch(const std::vector<double> displacements,
                        PIV::Set& original,
                        PIV::Set& stitched,
                        int nPics);
    
    void stitchPictureLoop(const std::vector<double> displacements,
                                PIV::Set& original,
                                PIV::Set& stitched,
                                std::vector<int>::iterator begin,
                           std::vector<int>::iterator end,
                           int iThread, int chunkSize);
    
    void stitchPicture(const std::vector<double> displacements,
                       PIV::Set& original,
                       PIV::Set& stitched,
                       int iPic);
    
    std::vector<double> getDisplacements(std::vector<int> order,
                                      PIV::Set& set);
    
    std::vector<double> getFrameDisplacements(const PIV::Frame& frameL,
                                           const PIV::Frame& frameR,
                                              bool withInterpolation);
    
    void naturalNeighbInterpolation
                        (const std::vector<double>& xin,
                         const std::vector<double>& yin,
                         const std::vector<double>& varin,
                         const std::vector<double>& xout,
                         const std::vector<double>& yout,
                         std::vector<double>& varout);
}






















