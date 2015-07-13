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

void PIV::Set::separateTecplotOutput(std::string fileName,
                                     std::string fileExt)
{
    std::stringstream ss;
    for (int iZone=0;
         iZone<static_cast<int>(this->loadedPictures.size());
         ++iZone)
     {
        ss.str("");
        ss.clear();
        ss << fileName << "_" << iZone << fileExt;
        std::ofstream tecFile;
        tecFile.open(ss.str());
        
        tecFile << "VARIABLES = \"X\" \"Y\" \"VX\" \"VY\"" << std::endl;
         tecFile << "ZONE T=\"" << 1 <<
            "\", I=" << this->container[iZone].frames[0].nx;
         tecFile << " ,J=" << this->container[iZone].frames[0].ny;
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
        
        tecFile.close();
     }
}

void PIV::Set::tecplotOut(std::string fileName)
{
    std::ofstream tecFile;
    tecFile.open(fileName);
    
    tecFile << "VARIABLES = \"X\" \"Y\" \"VX\" \"VY\"" << std::endl;
    std::cout << std::endl;
    for (int iZone=0;
         iZone<settings.nPics;
         ++iZone)
         {
             this->retrievePicture(iZone);
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
             std::cout << "\r"<< std::flush;
             std::cout << "  pict = " << iZone << std::flush;
             this->removePicture(iZone);
         }
         std::cout << std::endl;
    tecFile.close();
}


void PIV::Set::vorticityTecplotOut(std::string fileName)
{
    std::ofstream tecFile;
    tecFile.open(fileName);
    
    tecFile << "VARIABLES = \"X\" \"Y\" \"VX\" \"VY\" \"VORT\"" << std::endl;
    std::cout << std::endl;
    for (int iZone=0;
         iZone<settings.nPics;
         ++iZone)
         {
             this->retrievePicture(iZone);
             std::vector<double> vort =
                        this->container[iZone].frames[0].calculateVorticity();
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
                             <<vort[j*this->container[iZone].frames[0].ny + i]
                             << std::endl;
                 }
             }
             std::cout << "\r"<< std::flush;
             std::cout << "  pict = " << iZone << std::flush;
             this->removePicture(iZone);
         }
         std::cout << std::endl;
    tecFile.close();
}

void PIV::Set::PIV_like_xcorr(char vxOrVort)
{
    this->retrievePicture(0);
    int nx = this->container[0].frames[0].nx;
    int ny = this->container[0].frames[0].ny;

    double dx = this->container[0].frames[0].x[0] -
                this->container[0].frames[0].x[1];
    this->removePicture(0);
    
    std::ofstream ofile;
    ofile.open("tempXcoor.dat");
    ofile << "VARIABLES = \"x\"  \"y\" \"uc\" \"vc\" \"val\"" << std::endl;

    std::cout << std::endl;
    for (int n=1; n<=settings.maxN; ++n)
    {
        std::vector<double> xdisps(nx*ny, 0.0);
        std::vector<double> ydisps(nx*ny, 0.0);
        std::vector<double> xcorrvals(nx*ny, 0.0);
        for (int iPict=0; iPict<settings.nPics-n; ++iPict)
        {
            std::vector<double> txdisps;
            std::vector<double> tydisps;
            std::vector<double> txcorrvals;

            this->retrievePicture(iPict);
            this->retrievePicture(iPict + n);
            std::cout << "\r" << std::flush;
            std::cout << "  n=" << n << "  pict=" << iPict << std::flush;
            //piv_xcorr(&(this->container[iPict].frames[0].vx[0]),
                      //&(this->container[iPict + n].frames[0].vx[0]),
                      //&(ny),
                      //&(nx),
                      //&(settings.wNx),
                      //&(settings.wNy),
                      //&(txdisps[0]),
                      //&(tydisps[0]),
                      //&(txcorrvals[0]));

            switch (vxOrVort)
            {
                case 'v':
                    PIV::xxcorr(this->container[iPict].frames[0].vx,
                               this->container[iPict + n].frames[0].vx,
                               settings.wNx,
                               settings.wNy,
                               ny,
                               nx,
                               txdisps,
                               tydisps,
                               txcorrvals);
                    break;
                case 'w':
                    PIV::xxcorr(this->container[iPict].frames[0].calculateVorticity(),
                               this->container[iPict + n].frames[0].calculateVorticity(),
                               settings.wNx,
                               settings.wNy,
                               ny,
                               nx,
                               txdisps,
                               tydisps,
                               txcorrvals);
                    break;
            }

            this->removePicture(iPict);
            this->removePicture(iPict + n);

            for (int i=0; i<nx*ny; ++i)
            {
                xdisps[i] += txdisps[i];
                ydisps[i] += tydisps[i];
                xcorrvals[i] += txcorrvals[i];
            }
        }
        for (int i=0; i<nx*ny; ++i)
        {
            xdisps[i] /= (settings.nPics - n);
            ydisps[i] /= (settings.nPics - n);
            xcorrvals[i] /= (settings.nPics - n);
        }
        ofile << "ZONE T=\"n=" << n << "\" I=" << nx << " J=" << ny
              << " DATAPACKING=POINT" << std::endl;
        this->retrievePicture(0);
        for (int i=0; i<nx*ny; ++i)
        {
            ofile << this->container[0].frames[0].x[i] << "   "
                  << this->container[0].frames[0].y[i] << "   "
                  << -xdisps[i]*dx << "   "
                  << ydisps[i]*dx << "   "
                  << xcorrvals[i] << std::endl;
        }
        this->removePicture(0);
    }
    ofile.close();
}

void PIV::Set::timeXcorr()
{
    
    // TODO Ugly parameters
    int expon = 1;
    double offset = 0.0;
    int zero = 0;
    //==================
    this->retrievePicture(0);
    int nx = this->container[0].frames[0].nx;
    int ny = this->container[0].frames[0].ny;
    //std::cout << "nx = " << nx << std::endl;

    // limits determination
    int jmin;
    for (jmin=1; jmin<=nx; ++jmin)
    {
        if (this->container[0].frames[0].x[nx - jmin] > settings.minX) break;
    }
    this->removePicture(0);
    double normRxx = 0.0;
    std::cout << "normRxx = " << normRxx << std::endl;


    double temp;
    int counter;
    std::ofstream ofile;
    ofile.open("tempXcorr.dat");


    std::cout << std::endl;
    ofile << "VARIABLES = \"rDx\" \"y\" \"Rxx\"" << std::endl;
    for (int n=1; n<=settings.maxN; n+=1)
    // n loop for different time steps - ZONE
    {
        ofile << "ZONE T=\"n=" << n << "\" I=" << settings.maxR 
              << " J=" << settings.maxY-settings.minY
              << " DATAPACKING=POINT" << std::endl;
        for (int x2=settings.minY; x2<settings.maxY; ++x2)
        {
            // value for normalisation
            normRxx = 0.0;
            for (int i=0; i<settings.nFramesRxxNorm; ++i)
            {
                this->retrievePicture(i);
                normRxx += calculatexcorr(&(this->container[i].frames[0].vx[0]),
                                          &(this->container[i].frames[0].vx[0]),
                                          &(ny),
                                          &(nx),
                                          &(zero),
                                          &(expon),
                                          &(offset),
                                          &(jmin),
                                          &(x2));
                this->removePicture(i);
            }
            normRxx /= settings.nFramesRxxNorm;

            std::vector<int> rVec(settings.maxR, 0);
            std::vector<double> RxxVec(settings.maxR, 0.0);
            for (int i=0; i<settings.maxR; ++i)
            {
                rVec[i] = i*settings.deltaR;
                RxxVec[i] = 0.0; 
            }

            for (int t=0; t<settings.nPics-n; ++t)
            {
                std::cout << "\r" << std::flush;
                std::cout << "  n=" << n <<" --> pic=" << t << std::flush;
                this->retrievePicture(t);
                this->retrievePicture(t+n);
                for (int i=0; i<settings.maxR; ++i)
                {
                    temp = calculatexcorr(&(this->container[t].frames[0].vx[0]),
                                          &(this->container[t+n].frames[0].vx[0]),
                                          &(ny),
                                          &(nx),
                                          &(rVec[i]),
                                          &(expon),
                                          &(offset),
                                          &(jmin),
                                          &(x2));
                    if (temp > 0.0)
                    {
                        RxxVec[i] += temp/((settings.nPics-n));
                    }
                }
                this->removePicture(t+n);
                this->removePicture(t);
            }
            for (int i=0; i<settings.maxR; ++i)
            {
                ofile << rVec[i] << "    " << x2  << "    "
                      << RxxVec[i]/normRxx << std::endl;
            }
        }
    }
    ofile.close();
    std::cout << std::endl;
}
