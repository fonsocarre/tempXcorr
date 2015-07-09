//
//  Frame.h
//  Stitch
//
//  Created by Alfonso del Carre on 29/05/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

namespace PIV
{
    class Frame
    {
 
    public:
        Frame();
        Frame& operator = (Frame&& other); // move constructor
        Frame& operator = (const Frame& other);
        Frame(const Frame& other); // copy constructor
        Frame(const H5File& file, std::string dataset_name,
              int iFrame);
        ~Frame();
        
        void unload(H5File& file, std::string groupName);
        int findWallLocation();
        
        std::string dataset_name;
        hsize_t dims[2];
        
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> vx;
        std::vector<double> vy;
        
        int iFrame;
        int nx;
        int ny;
        int n;
    };
}