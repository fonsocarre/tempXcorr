//
//  Data.h
//  Stitch
//
//  Created by Alfonso del Carre on 28/05/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#pragma once

// hack for compiling HDF5 with CGAL
#include <mmintrin.h>

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>

#include "Frame.h"
#include "Picture.h"
#include "getSettings.h"
#include "xCorrModule.h"
#include "PIV_xcorr.h"

#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

namespace PIV
{
    class Set
    {
        std::string fileName;
        H5File file;
        double dt;
        int nx;
        int ny;
        int interrSize;
        
        
        
    public:
        std::map<int, PIV::Picture> container;
        std::set<int> loadedPictures;
        int nPictures;
        
        Set();
        Set(std::string fileName);
        //~Set();
        
        void retrievePicture(int nPicture);
        void removePicture(int nPicture);
        void insertPicture(PIV::Picture& pict);
        Set copyProperties();
        Set copyProperties(std::string fileName);
        void unloadPicture(int iPict);
        void closeFile();
        
        void outputAverageField(PIV::Frame avgFrame, std::string outFile);
        void invertCoordinates();
        H5File copyHdf5File();
        void setHdf5File(H5File file);
        PIV::Frame calculateAvgField();
        
        void tecplotOut(std::string fileName);
        void separateTecplotOutput(std::string fileName,
                                             std::string fileExt);

        void timeXcorr();
        void PIV_like_xcorr();
    };
}
















