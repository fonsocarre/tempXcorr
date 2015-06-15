//
//  Picture.h
//  Stitch
//
//  Created by Alfonso del Carre on 28/05/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#pragma once

#include <vector>
#include <sstream>

#include "Frame.h"

#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

namespace PIV
{
    class Picture
    {
        

    public:
        Picture();
        Picture& operator = (Picture&& other); // move constructor
        Picture(int iPicture, const H5File& file);
        Picture(int iPicture, int nFrames);
        ~Picture();
        
        void unload(H5File& file, std::string groupName);
        
        std::vector<PIV::Frame> frames;
        int nFrames;
        int iPicture;
        double t;
        
        void invertCoordinates();
        void cutFrames(int wallLoc, PIV::Picture& newP);
        void substractFrame(PIV::Frame& avgF);
        void correctFrameCoord(double xdisp, double ydisp);
    };
}
