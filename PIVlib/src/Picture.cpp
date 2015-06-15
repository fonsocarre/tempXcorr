//
//  Picture.cpp
//  Stitch
//
//  Created by Alfonso del Carre on 28/05/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#include "Picture.h"


PIV::Picture::Picture()
{
    //std::cerr << "Don't use this constructor!" << std::endl;
    //std::cerr << "From PIV::Picture:Picture() @ line "
    //          << __LINE__
    //          << std::endl;
}

PIV::Picture::Picture(int iPicture, const H5File& file)
{
    std::stringstream dataset_name;
    dataset_name << "/" << iPicture << "/";
    Group picture = file.openGroup(dataset_name.str());
    
    Attribute att = picture.openAttribute("nFrames");
    DataType type = PredType::NATIVE_INT;
    att.read(type, &(this->nFrames));

    att = picture.openAttribute("t");
    att.read(PredType::NATIVE_DOUBLE, &(this->t));
    
    
    this->frames.resize(this->nFrames);
    
    for (int iFrame=1; iFrame<=this->nFrames; ++iFrame)
    {
        dataset_name.str("");
        dataset_name.clear();
        dataset_name << "/" << iPicture << "/F" << iFrame << "/";
        this->frames[iFrame-1] = std::move(PIV::Frame(file,
                                     dataset_name.str(),
                                     iFrame));
    }
    
    this->iPicture = iPicture;
}

PIV::Picture::Picture(int iPicture, int nFrames)
{
    this->iPicture = iPicture;
    this->nFrames = nFrames;
    this->frames.resize(nFrames);
}


PIV::Picture::~Picture()
{
    this->frames.clear();
}


PIV::Picture& PIV::Picture::operator= (Picture&& other)
{
    this->nFrames = other.nFrames;
    this->iPicture = other.iPicture;
    
    this->frames = std::move(other.frames);
    
    return *this;
}

void PIV::Picture::unload(H5File& file, std::string groupName)
{
    Group group = file.createGroup(groupName);

    // nFrames
    DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
    Attribute attr = group.createAttribute("nFrames", PredType::NATIVE_INT,
                                           attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(this->nFrames));
    
    // t
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = group.createAttribute("t", PredType::NATIVE_DOUBLE,
                                 attr_dataspace);
    attr.write(PredType::NATIVE_DOUBLE, &(this->t));
    
    
    for (int iFrame=0; iFrame<nFrames; iFrame++)
    {
        std::stringstream name;
        name << groupName << "/F" << iFrame+1 ;
        this->frames[iFrame].unload(file, name.str());
    }
    
}


void PIV::Picture::invertCoordinates()
{
    for (auto& frame: this->frames)
    {
        for (int i=0; i<frame.n; ++i)
        {
            frame.x[i] *= -1;
            frame.vx[i] *= -1;
        }
    }
}

void PIV::Picture::cutFrames(int wallLoc, PIV::Picture& newP)
{
    
    newP.nFrames = this->nFrames;
    newP.iPicture = this->iPicture;
    newP.t = this->t;
    newP.frames.resize(newP.nFrames);
    
    for (int iFrame=0; iFrame<newP.nFrames; ++iFrame)
    {
        newP.frames[iFrame].nx = this->frames[iFrame].nx;
        newP.frames[iFrame].ny = wallLoc;
        newP.frames[iFrame].n = newP.frames[iFrame].nx*newP.frames[iFrame].ny;
        newP.frames[iFrame].iFrame = this->frames[iFrame].iFrame;
        
        newP.frames[iFrame].x.resize(newP.frames[iFrame].n);
        newP.frames[iFrame].y.resize(newP.frames[iFrame].n);
        newP.frames[iFrame].vx.resize(newP.frames[iFrame].n);
        newP.frames[iFrame].vy.resize(newP.frames[iFrame].n);
        
        int nx = newP.frames[iFrame].nx;
        int ny = newP.frames[iFrame].ny;
        
        for (int iRow=0; iRow<ny; ++iRow)
        {
            for (int iCol=0; iCol<nx; ++iCol)
            {
                newP.frames[iFrame].x[iRow*nx + iCol] =
                        this->frames[iFrame].x[iRow*nx + iCol];
                newP.frames[iFrame].y[iRow*nx + iCol] =
                        this->frames[iFrame].y[iRow*nx + iCol];
                newP.frames[iFrame].vx[iRow*nx + iCol] =
                        this->frames[iFrame].vx[iRow*nx + iCol];
                newP.frames[iFrame].vy[iRow*nx + iCol] =
                        this->frames[iFrame].vy[iRow*nx + iCol];
            }
        }
    }
}

void PIV::Picture::substractFrame(PIV::Frame& avgF)
{
    for (auto& frame: this->frames)
    {
        for (int i=0; i<frame.n; ++i)
        {
            frame.vx[i] -= avgF.vx[i];
            frame.vy[i] -= avgF.vy[i];
        }
    }
}

void PIV::Picture::correctFrameCoord(double xdisp, double ydisp)
{
    for (auto& frame: this->frames)
    {
        for (int i=0; i<frame.n; ++i)
        {
            frame.x[i] += xdisp;
            frame.y[i] += ydisp;
        }
    }
}


















