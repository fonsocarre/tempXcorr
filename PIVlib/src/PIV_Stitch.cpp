//
//
//  PIV_Stitch.cpp
//  Stitch
//
//  Created by Alfonso del Carre on 02/06/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#include "PIV_Stitch.h"

void PIV::parallelStitch(const std::vector<double> displacements,
                        PIV::Set& original,
                        PIV::Set& stitched,
                        int nPics)
{
    std::vector<std::thread> threads(settings.nThreads);
    const int chunkSize = nPics/settings.nThreads;
    
    std::vector<int> iPics;
    for (int i=0; i<nPics; ++i)
    {
        iPics.push_back(i);
    }
    auto iter = iPics.begin();
    int iThread = 0;
    for (auto it = std::begin(threads); it != std::end(threads)-1; ++it)
    {
        iThread++;
        *it = std::thread(PIV::stitchPictureLoop,
                          std::ref(displacements),
                          std::ref(original),
                          std::ref(stitched),
                          iter, iter+chunkSize,
                          iThread, chunkSize);
        iter += chunkSize;
    }
    iThread++;
    threads.back() = std::thread(PIV::stitchPictureLoop,
                                 std::ref(displacements),
                                 std::ref(original),
                                 std::ref(stitched), iter,
                                 iPics.end(),
                                 iThread, chunkSize);
    
    for(auto&&i: threads)
    {
        i.join();
    }
}

void PIV::stitchPictureLoop(const std::vector<double> displacements,
                            PIV::Set& original,
                            PIV::Set& stitched,
                            std::vector<int>::iterator begin,
                            std::vector<int>::iterator end,
                            int iThread, int chunkSize)
{
    int iPic = 0;
    for (auto it=begin; it!=end; ++it)
    {
        PIV::stitchPicture(displacements, original, stitched, *it);
        if (iThread == 1)
        {
            std::cout << "progress: " << ((iPic++)+1)/(0.01*chunkSize)<< "%"
            << std::endl;
        }
    }
}


void PIV::stitchPicture(const std::vector<double> displacements,
                           PIV::Set& original,
                           PIV::Set& stitched,
                           int iPic)
{
    m.lock();
        original.retrievePicture(iPic);
    m.unlock();
    
    int nCols = original.container[iPic].frames[0].nx;
    int nRows = original.container[iPic].frames[0].ny;
    int nElems= nCols*nRows;
    
    double dx0 = std::abs(original.container[iPic].frames[0].x[1] -
                          original.container[iPic].frames[0].x[0]);
    double dx1 = std::abs(original.container[iPic].frames[1].x[1] -
                          original.container[iPic].frames[1].x[0]);
    
    double dx = 0.5*(dx0 + dx1);
    
    double xminL, xminR;
    double xmaxL, xmaxR;
    double yminL, yminR;
    double ymaxL, ymaxR;
    
    xminL = *std::min_element(original.container[iPic].frames[0].x.begin(),
                              original.container[iPic].frames[0].x.end());
    xmaxL = *std::max_element(original.container[iPic].frames[0].x.begin(),
                              original.container[iPic].frames[0].x.end());
    yminL = *std::min_element(original.container[iPic].frames[0].y.begin(),
                              original.container[iPic].frames[0].y.end());
    ymaxL = *std::max_element(original.container[iPic].frames[0].y.begin(),
                              original.container[iPic].frames[0].y.end());
    xminR = *std::min_element(original.container[iPic].frames[1].x.begin(),
                              original.container[iPic].frames[1].x.end());
    xmaxR = *std::max_element(original.container[iPic].frames[1].x.begin(),
                              original.container[iPic].frames[1].x.end());
    yminR = *std::min_element(original.container[iPic].frames[1].y.begin(),
                              original.container[iPic].frames[1].y.end());
    ymaxR = *std::max_element(original.container[iPic].frames[1].y.begin(),
                              original.container[iPic].frames[1].y.end());
    
    for (auto& x: original.container[iPic].frames[1].x)
    {
        x += displacements[0]*dx - xminR + xmaxL;
    }
    for (auto& y: original.container[iPic].frames[1].y)
    {
        y += displacements[1]*dx - yminR + yminL;
    }
    
    xminR = *std::min_element(original.container[iPic].frames[1].x.begin(),
                              original.container[iPic].frames[1].x.end());
    xmaxR = *std::max_element(original.container[iPic].frames[1].x.begin(),
                              original.container[iPic].frames[1].x.end());
    yminR = *std::min_element(original.container[iPic].frames[1].y.begin(),
                              original.container[iPic].frames[1].y.end());
    ymaxR = *std::max_element(original.container[iPic].frames[1].y.begin(),
                              original.container[iPic].frames[1].y.end());
    
    std::vector<double> xin(2*nElems);
    std::vector<double> yin(2*nElems);
    std::vector<double> vxin(2*nElems);
    std::vector<double> vyin(2*nElems);
    
    for (int i=0; i<nElems; ++i)
    {
        xin[i]  = original.container[iPic].frames[0].x[i];
        yin[i]  = original.container[iPic].frames[0].y[i];
        vxin[i] = original.container[iPic].frames[0].vx[i];
        vyin[i] = original.container[iPic].frames[0].vy[i];
    }
    for (int i=nElems; i<2*nElems; ++i)
    {
        xin[i]  = original.container[iPic].frames[1].x[i-nElems];
        yin[i]  = original.container[iPic].frames[1].y[i-nElems];
        vxin[i] = original.container[iPic].frames[1].vx[i-nElems];
        vyin[i] = original.container[iPic].frames[1].vy[i-nElems];
    }
//    for (int i=0; i<nRows; ++i)
//    {
//        for (int j=0; j<2*nCols; ++j)
//        {
//            if (j<nCols)
//            // first frame
//            {
//                xin[i*2*nCols + j]  = original.container[iPic].frames[0].
//                                                x[i*nCols + j];
//            } else
//            // second frame
//            {
//                xin[i*2*nCols + j]  = original.container[iPic].frames[0].
//                                                x[i*nCols + j];
//            }
//        }
//    }
    
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    
//    xmin = *std::min_element(xin.begin(), xin.end());
//    xmax = *std::max_element(xin.begin(), xin.end());
//    ymin = *std::min_element(yin.begin(), yin.end());
//    ymax = *std::max_element(yin.begin(), yin.end());

    xmin = xminL;
    xmax = xmaxR;
    ymin = fmax(yminL, yminR);
    ymax = fmin(ymaxL, ymaxR);
    
    int newnRows = floor((ymax-ymin)/dx);
    int newnCols = floor((xmax-xmin)/dx);
    
    PIV::Frame newFrame;
    newFrame.nx = newnCols;
    newFrame.ny = newnRows;
    newFrame.n  = newnCols*newnRows;
    newFrame.iFrame = 1;
    newFrame.dataset_name = original.container[iPic].frames[0].dataset_name;
    newFrame.x.resize(newFrame.n);
    newFrame.y.resize(newFrame.n);
    newFrame.vx.resize(newFrame.n);
    newFrame.vy.resize(newFrame.n);
    for (int i=0; i<newnRows; ++i)
    {
        for (int j=0; j<newnCols; ++j)
        {
            newFrame.x[i*newnCols + j] = j*dx + xmin;
        }
    }
    for (int i=0; i<newnRows; ++i)
    {
        for (int j=0; j<newnCols; ++j)
        {
            newFrame.y[i*newnCols + j] = ymax - i*dx;
        }
    }
    
    double xdisp = displacements[0];
    double ydisp = displacements[1];
    interpolate(&(original.container[iPic].frames[0].x[0]),
                &(original.container[iPic].frames[0].y[0]),
                &(original.container[iPic].frames[0].vx[0]),
                &(original.container[iPic].frames[0].vy[0]),
                &(original.container[iPic].frames[1].x[0]),
                &(original.container[iPic].frames[1].y[0]),
                &(original.container[iPic].frames[1].vx[0]),
                &(original.container[iPic].frames[1].vy[0]),
                &(nRows),
                &(nCols),
                &(newFrame.x[0]),
                &(newFrame.y[0]),
                &(newFrame.vx[0]),
                &(newFrame.vy[0]),
                &(newnRows),
                &(newnCols),
                &(xdisp),
                &(ydisp));
    
//    PIV::naturalNeighbInterpolation(
//                                    std::ref(xin),
//                                    std::ref(yin),
//                                    std::ref(vxin),
//                                    std::ref(newFrame.x),
//                                    std::ref(newFrame.y),
//                                    std::ref(newFrame.vx));
//    PIV::naturalNeighbInterpolation(
//                                    std::ref(xin),
//                                    std::ref(yin),
//                                    std::ref(vyin),
//                                    std::ref(newFrame.x),
//                                    std::ref(newFrame.y),
//                                    std::ref(newFrame.vy));
    
    PIV::Picture newPict;
    newPict.frames.push_back(newFrame);
    newPict.nFrames = 1;
    newPict.iPicture = iPic;
    
    m.lock();
        stitched.insertPicture(newPict);
    m.unlock();
    
    m.lock();
        stitched.unloadPicture(iPic);
    m.unlock();
    
    m.lock();
        original.removePicture(iPic);
    m.unlock();
}


void PIV::naturalNeighbInterpolation
                                (const std::vector<double>& xin,
                                 const std::vector<double>& yin,
                                 const std::vector<double>& varin,
                                 const std::vector<double>& xout,
                                 const std::vector<double>& yout,
                                 std::vector<double>& varOut)
{
    int nIn  = static_cast<int> (xin.size());
    int nOut = static_cast<int> (xout.size());
    
    //std::vector<double> varOut;
    varOut.reserve(nOut);
    
    Delaunay_triangulation T;
    std::map<Point, Coord_type, K::Less_xy_2> function_values;
    typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xy_2 > >
                                                            Value_access;
    // input data insertion
    //std::cout << "Inserting points..." << std::endl;
    for (int i=0; i<nIn; ++i)
    {
        K::Point_2 p(xin[i], yin[i]);
        T.insert(p);
        function_values.insert(std::make_pair(p, varin[i]));
    }
    
//     coordinate computation and interpolation
    //std::cout << "Interpolating to " << nOut << " points..." << std::endl;
    for (int i=0; i<nOut; ++i)
    {
        //std::cout << i << std::endl;
        K::Point_2 p(xout[i], yout[i]);
        std::vector<std::pair<Point,Coord_type>> coords;
        
        Coord_type norm = CGAL::natural_neighbor_coordinates_2
                        ( T, p, std::back_inserter(coords)).second;
        
        Coord_type res = CGAL::linear_interpolation(coords.begin(), coords.end(),
                                                    norm, Value_access(function_values));
        varOut.push_back(res);
    }
    //std::cout << " ...DONE" << std::endl;
    
    //return varOut;
}


std::vector<double> PIV::getDisplacements(std::vector<int> order,
                                        PIV::Set& set)
{
    int nPics = settings.nPicsForStitch;
    std::cout << "Getting displacements from " << nPics << " pictures..."
            << std::endl;
    std::vector<std::vector<double>> displacements;
    displacements.resize(nPics);
    for (int iPict=0; iPict<nPics; ++iPict)
    {
        set.retrievePicture(iPict);
        displacements[iPict].resize(2);
        displacements[iPict] = PIV::getFrameDisplacements(
                                      set.container[iPict].frames[order[0]],
                                      set.container[iPict].frames[order[1]],
                                                          false);
        set.removePicture(iPict);
    }
    std::cout << std::endl;
    // find the mode
    // TODO change to min_element and max_element
    int min1 = 100000;
    int max1 =-100000;
    int min2 = 100000;
    int max2 =-100000;
    for (auto& elem: displacements)
    {
        if (elem[0] < min1) {min1=elem[0];}
        if (elem[0] > max1) {max1=elem[0];}
        if (elem[1] < min2) {min2=elem[1];}
        if (elem[1] > max2) {max2=elem[1];}
    }
    std::vector<int> hist1(max1-min1+1);
    std::vector<int> hist2(max2-min2+1);
    for (auto& elem: displacements)
    {
        hist1[elem[0]-min1]++;
        hist2[elem[1]-min2]++;
    }
    int mode1 = static_cast<int>(std::max_element(hist1.begin(), hist1.end())
                                                - hist1.begin()) + min1;
    int mode2 = static_cast<int>(std::max_element(hist2.begin(), hist2.end())
                                                - hist2.begin()) + min2;
    
    
    double xdisp = 0.0;
    double ydisp = 0.0;
    std::vector<double> temp;
    for (int i=0; i<nPics; ++i)
    {
        if (displacements[i][0] == mode1)
        {
            set.retrievePicture(i);
            temp = PIV::getFrameDisplacements(
                                       set.container[i].frames[order[0]],
                                       set.container[i].frames[order[1]],
                                       true);
            xdisp = temp[0];
            set.removePicture(i);
            break;
        }
    }
    for (int i=0; i<nPics; ++i)
    {
        if (displacements[i][1] == mode2)
        {
            set.retrievePicture(i);
            temp = PIV::getFrameDisplacements(
                                              set.container[i].frames[order[0]],
                                              set.container[i].frames[order[1]],
                                              true);
            ydisp = temp[1];
            set.removePicture(i);
            break;
        }
    }
    
    std::cout << "The displacement of the right frame will be " << std::endl;
    std::cout << "x:" << xdisp << " y:" << ydisp << std::endl;
    std::cout << "  in interrogation window lengths" << std::endl;
    
    return std::vector<double>{xdisp, ydisp};
}

std::vector<double> PIV::getFrameDisplacements(const PIV::Frame& frameL,
                                               const PIV::Frame& frameR,
                                               bool withInterpolation)
{
    int rowSize = frameL.nx;
    int colSize = frameL.ny;
    int overlapSize = static_cast<int>(settings.maxXdisp*rowSize);
    
    std::vector<double> vecL(rowSize*colSize);
    std::vector<double> vecR(rowSize*colSize);
    
    for (int i=0; i<rowSize*colSize; ++i)
    {
        vecL[i] = sqrt(pow(frameL.vx[i],2) + pow(frameL.vy[i],2));
        vecR[i] = sqrt(pow(frameR.vx[i],2) + pow(frameR.vy[i],2));
    }
    
    //const int DIM = 2;
    double xdisp;
    double ydisp;
//    xcorr(&(frameL.vx[0]),
//          &(frameL.vx[0]),
//          &(colSize), &(rowSize), &(overlapSize), &(xdisp), &(ydisp));
    xcorr(&(vecL[0]),
          &(vecR[0]),
          &(colSize), &(rowSize), &(overlapSize), &(xdisp), &(ydisp), &(withInterpolation));
    return std::vector<double> {xdisp, ydisp};
}
