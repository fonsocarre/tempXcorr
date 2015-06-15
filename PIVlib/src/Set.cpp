//
//  Set.cpp
//  Stitch
//
//  Created by Alfonso del Carre on 28/05/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//


#include "Set.h"

PIV::Set::Set()
{
    
}

//PIV::Set::~Set()
//{
//    this->loadedPictures.clear();
//    this->container.clear();
//    this->file.close();
//}

PIV::Set::Set(std::string fileName): file(fileName, H5F_ACC_RDONLY)
{

    this->fileName = fileName;
    std::cout << "Opening file: " << fileName << std::endl;
    
    Group root = this->file.openGroup("/");
    
    Attribute att = root.openAttribute("dt");
    DataType type = PredType::NATIVE_DOUBLE;
    att.read(type, &(this->dt));
    std::cout << "    dt = " << this->dt << std::endl;
    
    att = root.openAttribute("nx");
    type = PredType::NATIVE_INT;
    att.read(type, &(this->nx));
    std::cout << "    nx = " << this->nx << std::endl;
    
    att = root.openAttribute("ny");
    type = PredType::NATIVE_INT;
    att.read(type, &(this->ny));
    std::cout << "    ny = " << this->ny << std::endl;
    
    att = root.openAttribute("interrogationWinSize");
    type = PredType::NATIVE_INT;
    att.read(type, &(this->interrSize));
    std::cout << "    interrogation win size = " << this->interrSize << std::endl;
    
    att = root.openAttribute("nCaptures");
    type = PredType::NATIVE_INT;
    att.read(type, &(this->nPictures));
    std::cout << "    nPictures = " << this->nPictures << std::endl;
}

void PIV::Set::retrievePicture(int nPicture)
{
    this->container[nPicture] = std::move(PIV::Picture(nPicture+1, this->file));
    this->loadedPictures.insert(nPicture);
}

void PIV::Set::removePicture(int nPicture)
{
    auto it = this->container.find(nPicture);
    this->container.erase(it);
    this->loadedPictures.erase(nPicture);
}

PIV::Set PIV::Set::copyProperties()
{
    PIV::Set result;
    
    result.fileName = this->fileName + "new.h5";
    std::cout << "Opening file: " << result.fileName << std::endl;
    result.file = H5Fcreate(result.fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    result.nPictures = this->nPictures;
    result.dt = this->dt;
    result.nx = this->nx;
    result.ny = this->ny;
    result.interrSize = this->interrSize;

    
    // nCaptures
    DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
    Attribute attr = result.file.createAttribute("nCaptures", PredType::NATIVE_INT, 
                                                attr_dataspace);
    attr.write(PredType::NATIVE_INT , &(result.nPictures));
    
    // dt
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = result.file.createAttribute("dt", PredType::NATIVE_DOUBLE, 
                                        attr_dataspace);
    attr.write(PredType::NATIVE_DOUBLE, &(result.dt));
    
    // nx
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = result.file.createAttribute("nx", PredType::NATIVE_INT, attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(result.nx));
    
    // ny
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = result.file.createAttribute("ny", PredType::NATIVE_INT, attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(result.ny));
    
    // interrogationWinSize
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = result.file.createAttribute("interrogationWinSize", PredType::NATIVE_INT,
                                        attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(result.interrSize));
    
    return result;
}

PIV::Set PIV::Set::copyProperties(std::string fileName)
{
    PIV::Set result;
    
    result.fileName = fileName;
    std::cout << "Opening file: " << result.fileName << std::endl;
    result.file = H5Fcreate(result.fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    result.nPictures = this->nPictures;
    result.dt = this->dt;
    result.nx = this->nx;
    result.ny = this->ny;
    result.interrSize = this->interrSize;
    
    
    // nCaptures
    DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
    Attribute attr = result.file.createAttribute("nCaptures", PredType::NATIVE_INT, 
                                                 attr_dataspace);
    attr.write(PredType::NATIVE_INT , &(result.nPictures));
    
    // dt
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = result.file.createAttribute("dt", PredType::NATIVE_DOUBLE, 
                                       attr_dataspace);
    attr.write(PredType::NATIVE_DOUBLE, &(result.dt));
    
    // nx
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = result.file.createAttribute("nx", PredType::NATIVE_INT, attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(result.nx));
    
    // ny
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = result.file.createAttribute("ny", PredType::NATIVE_INT, attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(result.ny));
    
    // interrogationWinSize
    attr_dataspace = DataSpace(H5S_SCALAR);
    attr = result.file.createAttribute("interrogationWinSize", PredType::NATIVE_INT,
                                       attr_dataspace);
    attr.write(PredType::NATIVE_INT, &(result.interrSize));
    
    return result;
}


void PIV::Set::insertPicture(PIV::Picture& pict)
{
    if (this->loadedPictures.find(pict.iPicture) != this->loadedPictures.end())
    {
        std::cerr << "Problem!!! adding twice the same pict to the set!!!"
                  << std::endl
                  << "In Set.cpp @ line " << __LINE__ << std::endl;
    }
    this->container[pict.iPicture-1] = std::move(pict);
    this->loadedPictures.insert(pict.iPicture-1);
}

void PIV::Set::unloadPicture(int iPict)
{
    if (this->loadedPictures.find(iPict) == this->loadedPictures.end())
    // element not in memory
    {
        std::cerr << "Error, tried to unload a picture that does not exist"
                    << std::endl;
        std::cerr << "    in Set::unloadPicture @ line " << __LINE__ 
                    << std::endl;
    }
    
    std::stringstream groupName;
    groupName << "/" << iPict+1 ;
    //this->file.createGroup(groupName.str());
    this->container[iPict].unload(this->file, groupName.str());
}

void PIV::Set::closeFile()
{
    this->file.close();
}

PIV::Frame PIV::Set::calculateAvgField()
{
    PIV::Frame avgFrame;
    
    avgFrame.n  = 0;
    avgFrame.nx = 0;
    avgFrame.ny = 0;
    
    for (int iPict=0; iPict<this->nPictures; ++iPict)
    {
        this->retrievePicture(iPict);
        if (iPict == 0)
        {
            avgFrame.n = this->container[0].frames[0].n;
            avgFrame.nx = this->container[0].frames[0].nx;
            avgFrame.ny = this->container[0].frames[0].ny;
            // Allocation
            avgFrame.x.resize(avgFrame.n);
            avgFrame.y.resize(avgFrame.n);
            avgFrame.vx.resize(avgFrame.n);
            avgFrame.vy.resize(avgFrame.n);
            for (int i=0; i<avgFrame.n; ++i)
            {
                avgFrame.vx[i] = 0.0;
                avgFrame.vy[i] = 0.0;
            }
            avgFrame.x = this->container[0].frames[0].x;
            avgFrame.y = this->container[0].frames[0].y;
        }
        for (int i=0; i<avgFrame.n; ++i)
        {
            avgFrame.vx[i] += this->container[iPict].frames[0].vx[i];
            avgFrame.vy[i] += this->container[iPict].frames[0].vy[i];
        }
        
        this->removePicture(iPict);
        std:: cout << "   progress: " << (iPict)/(this->nPictures*0.01) 
        << std::endl;
    }
    
    for (int i=0; i<avgFrame.n; ++i)
    {
        avgFrame.vx[i] /= this->nPictures;
        avgFrame.vy[i] /= this->nPictures;
    }

    
    return avgFrame;
}

void PIV::Set::outputAverageField(PIV::Frame avgFrame, std::string outFile)
{
    
    //std::cout << "Obtaining average fields..." << std::endl;
    std::ofstream output;
    output.open(outFile);

//    std::vector<double> x;
//    std::vector<double> y;
//    std::vector<double> vx;
//    std::vector<double> vy;
//
//    //std::shared_ptr<PIV::Frame> frame = nullptr;
//
//    int n = 0;
//    int nx = 0;
//    int ny = 0;
//    for (int iPict=0; iPict<this->nPictures; ++iPict)
//    {
//        this->retrievePicture(iPict);
//        if (iPict == 0)
//        {
//            n = this->container[0].frames[0].n;
//            nx = this->container[0].frames[0].nx;
//            ny = this->container[0].frames[0].ny;
//            // Allocation
//            x.resize(n);
//            y.resize(n);
//            vx.resize(n);
//            vy.resize(n);
//            for (int i=0; i<n; ++i)
//            {
//                vx[i] = 0.0;
//                vy[i] = 0.0;
//            }
//            x = this->container[0].frames[0].x;
//            y = this->container[0].frames[0].y;
//        }
//        for (int i=0; i<n; ++i)
//        {
//            vx[i] += this->container[iPict].frames[0].vx[i];
//            vy[i] += this->container[iPict].frames[0].vy[i];
//        }
//
//        this->removePicture(iPict);
//        std:: cout << "   progress: " << (iPict)/(this->nPictures*0.01) 
//                   << std::endl;
//    }
//
//    for (int i=0; i<n; ++i)
//    {
//        vx[i] /= this->nPictures;
//        vy[i] /= this->nPictures;
//    }
//    
    
    output << "VARIABLES = \"X\" \"Y\" \"VX\" \"VY\"" << std::endl;
    int iZone = 0;
     output << "ZONE T=\"" << iZone <<
        "\", I=" << avgFrame.nx;
     output << " ,J=" << avgFrame.ny;
     output << " ,SOLUTIONTIME = " << iZone;
     output << " ,DATAPACKING = POINT" << std::endl;
    
     for (int j=0; j<avgFrame.ny; ++j)
     {
         for (int i=0; i<avgFrame.nx; ++i)
         {
             output << avgFrame.x[j*avgFrame.nx + i]
                     << "  "
                     << avgFrame.y[j*avgFrame.nx + i]
                     << "  "
                     << avgFrame.vx[j*avgFrame.nx + i]
                     << "  "
                     << avgFrame.vy[j*avgFrame.nx + i]
                     << std::endl;
         }
     }
    
    output.close();
}


H5File PIV::Set::copyHdf5File()
{
    return this->file;
}

void PIV::Set::setHdf5File(H5File file)
{
    this->file = file;
}


void PIV::Set::tecplotOut(std::string fileName)
{
    std::ofstream tecFile;
    tecFile.open(fileName);
    
    tecFile << "VARIABLES = \"X\" \"Y\" \"VX\" \"VY\"" << std::endl;
    for (int iZone=0;
         iZone<static_cast<int>(this->loadedPictures.size());
         ++iZone)
         {
             tecFile << "ZONE T=\"" << iZone <<
                "\", I=" << this->container[iZone].frames[0].nx;
             tecFile << " ,J=" << this->container[iZone].frames[0].ny;
             tecFile << " ,SOLUTIONTIME = " << iZone;
             tecFile << " ,DATAPACKING = POINT" << std::endl;
             for (int j=0; j<this->container[iZone].frames[0].nx; ++j)
             {
                 for (int i=0; i<this->container[iZone].frames[0].ny; ++i)
                 {
                     tecFile <<this->container[iZone].frames[0].x
                                [j*this->container[iZone].frames[0].ny + i]
                             << "  "
                             <<this->container[iZone].frames[0].y
                                [j*this->container[iZone].frames[0].ny + i]
                             << "  "
                             <<this->container[iZone].frames[0].vx
                                [j*this->container[iZone].frames[0].ny + i]
                             << "  "
                             <<this->container[iZone].frames[0].vy
                                [j*this->container[iZone].frames[0].ny + i]
                             << "  "
                             << std::endl;
                 }
             }
         }
    
    tecFile.close();
}


void PIV::Set::timeXcorr()
{
    int pict0 = settings.tempXcorrMinPict;
    this->retrievePicture(pict0);
    
    // Ugly parameters
    int xmin = 1;
    int ymin = 0;
    int ymax = 0;
    int expon = 2;
    double offset = 11.6;
    //==================
    
    int nx = this->container[pict0].frames[0].nx;
    int ny = this->container[pict0].frames[0].ny;

    std::vector<int> xdisps(settings.tempXcorrMaxX, 0);
    std::vector<int> ydisps(settings.tempXcorrMaxX, 0);
    std::vector<double> corrVal(settings.tempXcorrMaxX, 0.0);

    for (int iPict=pict0+1;
         iPict<settings.tempXcorrMaxPict;
         ++iPict)
    {
        this->retrievePicture(iPict);
        
        tempXcorr(&(this->container[pict0].frames[0].vx[0]),
                  &(this->container[iPict].frames[0].vx[0]),
                  &(xmin),
                  &(settings.tempXcorrMaxX),
                  &(ymin),
                  &(ymax),
                  &(offset),
                  &(expon),
                  &(ny),
                  &(nx),
                  &(xdisps[0]),
                  &(ydisps[0]),
                  &(corrVal[0]));

        this->removePicture(iPict);
    }
}













