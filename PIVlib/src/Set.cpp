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
    if (this->loadedPictures.find(nPicture) == this->loadedPictures.end())
    // element not in memory
    {
        std::cerr << "Error, tried to remove a picture that does not exist"
                    << std::endl;
        std::cerr << "    in Set::removePicture@ line " << __LINE__ 
                    << std::endl;
        std::cerr << "The pict is iPict = " << nPicture << std::endl;
    }
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
    if (this->loadedPictures.find(pict.iPicture-1) != this->loadedPictures.end())
    {
        std::cerr << "Problem!!! adding twice the same pict to the set!!!"
                  << std::endl
                  << "In Set.cpp @ line " << __LINE__ << std::endl;
    }
    //std::cout << "Inserting pict index = " << pict.iPicture-1 << std::endl;
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
        std::cerr << "The pict is iPict = " << iPict << std::endl;
    }
    
    std::stringstream groupName;
    groupName << "/" << iPict+1 ;
    //this->file.createGroup(groupName.str());
    this->container[iPict].unload(this->file, groupName.str());
    this->removePicture(iPict);
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
    
    for (int iPict=0; iPict<settings.nPics; ++iPict)
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
        std:: cout << "   progress: " << (iPict)/(settings.nPics*0.01) 
        << std::endl;
    }
    
    for (int i=0; i<avgFrame.n; ++i)
    {
        avgFrame.vx[i] /= settings.nPics;
        avgFrame.vy[i] /= settings.nPics;
    }

    
    return avgFrame;
}

void PIV::Set::outputAverageField(PIV::Frame avgFrame, std::string outFile)
{
    
    //std::cout << "Obtaining average fields..." << std::endl;
    std::ofstream output;
    output.open(outFile);

    
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
    std::cout << std::endl;
    for (int iZone=0;
         //iZone<static_cast<int>(this->loadedPictures.size());
         iZone<settings.nPics;
         ++iZone)
     {
        this->retrievePicture(iZone);
        ss.str("");
        ss.clear();
        ss << fileName << "_" << iZone << fileExt;
        std::ofstream tecFile;
        tecFile.open(ss.str());
        std::cout << "\r" << std::flush;
        std::cout << "  zone = " << iZone << std::flush;
        
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
        this->removePicture(iZone); 
        tecFile.close();
     }
     std::cout << std::endl;
}

void PIV::Set::tecplotOut(std::string fileName)
{
    std::ofstream tecFile;
    tecFile.open(fileName);
    
    double dt = ((double) 1)/1400.;
    tecFile << "VARIABLES = \"X\" \"Y\" \"VX\" \"VY\"" << std::endl;
    std::cout << std::endl;
    for (int iZone=0;
         iZone<settings.nPics;
         ++iZone)
         {
             this->retrievePicture(iZone);
             tecFile << "ZONE T=\"" << "DATA" <<
                "\", I=" << this->container[iZone].frames[0].nx;
             tecFile << " ,J=" << this->container[iZone].frames[0].ny;
             tecFile << " ,SOLUTIONTIME = " << iZone*dt;
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
    ofile.open("tempXcorr.dat");
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

void PIV::Set::tempRetrievePict(int iPict, 
                                std::vector<double>& vxVec,
                                std::vector<double>& vyVec)
{
    std::stringstream dataset_name;
    dataset_name << iPict+1 << "/F1/";
    hsize_t dims[2];
    auto type = PredType::NATIVE_DOUBLE;

    DataSet dataset;
    DataSpace dataspace;

    //vx
    dataset = file.openDataSet(dataset_name.str() + "vx");
    dataspace = dataset.getSpace();
    //this->type = this->dataset.getTypeClass();
    dataspace.getSimpleExtentDims(dims, NULL);
    vxVec.resize(dims[0]);
    // TODO write type inference instead of just the macro
    dataset.read(&(vxVec[0]), type, dataspace);
    
    //vy
    dataset = file.openDataSet(dataset_name.str() + "vy");
    dataspace = dataset.getSpace();
    //this->type = this->dataset.getTypeClass();
    dataspace.getSimpleExtentDims(dims, NULL);
    vyVec.resize(dims[0]);
    // TODO write type inference instead of just the macro
    dataset.read(&(vyVec[0]), type, dataspace);
}

std::vector<double> PIV::Set::obtainTimeSeries(int tmin,
                                               int tmax,
                                               int index)
{
    std::vector<double> timeSeries(tmax-tmin+1, 0.0);
    std::vector<double> vx;
    std::vector<double> vy;
    for (int t=tmin; t<=tmax; ++t)
    {
        this->tempRetrievePict(t, vx, vy);       
        timeSeries[t] = vx[index];
    }
    return timeSeries;
}

std::vector<std::vector<float>> PIV::Set::preloadData(char variable, 
                                int tmin, 
                                int tmax)
{
    std::cout << "Loading data into memory..." << std::endl;
    std::vector<std::vector<float>> data;
    std::vector<double> tempVecvx;
    std::vector<double> tempVecvy;
    data.resize(tmax-tmin+1);
    this->retrievePicture(0);
    //std::cout << __LINE__ << std::endl;
    int n = this->container[0].frames[0].n;
    //this->removePicture(0);
    //std::cout << __LINE__ << std::endl;
    for (int t=tmin; t<=tmax; ++t)
    {
        data[t].resize(n);
        //std::cout << __LINE__ << std::endl;
        this->tempRetrievePict(t, tempVecvx, tempVecvy);
    
        switch (variable)
        {
            case 'v':
                for (int i=0; i<n; ++i)
                {
                    data[t][i] = static_cast<float>(tempVecvx[i]);
                }
                break;

            case 'w':
                data[t] = this->calculateVorticity(tempVecvx, tempVecvy);
                break;
        }
    }
    std::cout << "  data loaded" << std::endl;

    return data;

}


std::vector<float> PIV::Set::calculateVorticity(std::vector<double>& vx,
                                                 std::vector<double>& vy)
{
    //std::cout << "calculateVorticity: " << __LINE__ << std::endl;
    this->retrievePicture(0);
    std::vector<float> vort(this->container[0].frames[0].x.size(), 0.0);

    //std::cout << "calculateVorticity: " << __LINE__ << std::endl;
    const double dx = std::abs(this->container[0].frames[0].x[1] - 
                               this->container[0].frames[0].x[0]);
    int nRows = this->container[0].frames[0].ny;
    int nCols = this->container[0].frames[0].nx;
    // variable for knowing if coords increase with index (+1) or 
    // decrease (-1)
    const int signx = (this->container[0].frames[0].x[1] -
                        this->container[0].frames[0].x[0]>0.0? 1: -1);
    const int signy = (this->container[0].frames[0].y[this->container[0].frames[0].nx]
                        - this->container[0].frames[0].y[0]>0.0? 1: -1);

    //std::cout << "calculateVorticity: " << __LINE__ << std::endl;
    // this implementation uses the circulation way of calculating the
    // vorticity
    // C = \int \vec{u}\cdot d\vec{l} = \int_A \omega dA
    double elemVort;
    std::vector<int> locations;
    for (int iElem=0; iElem < this->container[0].frames[0].ny-1; ++iElem)
    {
        for (int jElem=0; jElem < this->container[0].frames[0].nx-1; ++jElem)
        {
            //std::cout << "calculateVorticity: " << __LINE__ << std::endl;
            elemVort = 0.0;
            // side 11 (S)
            elemVort += 0.5*(vx[cij(iElem + 1, jElem, nRows, nCols)] +
                             vx[cij(iElem + 1, jElem+1, nRows, nCols)]) *
                                    dx * signx;
            
            //std::cout << "calculateVorticity: " << __LINE__ << std::endl;
            //std::cout << "n = " << nRows*nCols << std::endl;
            //std::cout << "nCols = " << nCols << " nRows = " << nRows << std::endl;
            //std::cout << "cij(iElem + 1, jElem + 1, nRows, nCols) = " <<
                           //cij(iElem + 1, jElem + 1, nRows, nCols) << std::endl; 
            //std::cout << "cij(iElem, jElem + 1, nRows, nCols) = " <<
                           //cij(iElem, jElem + 1, nRows, nCols) << std::endl; 
            //std::cout << "iElem = " << iElem << " jElem = " << jElem << std::endl;
            //std::cout << "size(vy) = " << static_cast<int>(vy.size()) << std::endl;
            //std::cout << "size(vx) = " << static_cast<int>(vx.size()) << std::endl;
            // side 12 (E)
            elemVort += 0.5*(vy[cij(iElem + 1, jElem + 1, nRows, nCols)] +
                             vy[cij(iElem, jElem + 1, nRows, nCols)]) *
                                    dx * signy;

            //std::cout << "calculateVorticity: " << __LINE__ << std::endl;
            // side 13 (N)
            elemVort += -0.5*(vx[cij(iElem, jElem + 1, nRows, nCols)] +
                              vx[cij(iElem, jElem, nRows, nCols)]) *
                                    dx * signx;

            //std::cout << "calculateVorticity: " << __LINE__ << std::endl;
            // side 14 (W)
            elemVort += -0.5*(vy[cij(iElem, jElem, nRows, nCols)] +
                              vy[cij(iElem + 1, jElem, nRows, nCols)]) *
                                    dx * signy;

            //std::cout << "calculateVorticity: " << __LINE__ << std::endl;
            // interpolation to the nodes based on if corner (1, 2, 3, 4),
            // vertex (11, 12, 13, 14), or interior (-1)
            locations = this->container[0].frames[0].determinePositions(iElem, jElem);
            std::vector<double> coeffs(4, 0.25);
            for (int iCoeff=0; iCoeff<4; ++iCoeff)
            {
                switch (locations[iCoeff])
                {
                    case -1:
                        break;
                    case 1:
                    case 2:
                    case 3:
                    case 4:
                        coeffs[iCoeff] = 1.0;
                        break;
                    case 11:
                    case 12:
                    case 13:
                    case 14:
                        coeffs[iCoeff] = 0.5;
                        break;
                    default:
                        std::cerr << "bad input in line" << __LINE__
                                                     << std::endl;
                }
            }

            // inserting the values in the vort matrix
            // i+1,j
            vort[cij(iElem+1, jElem, nRows, nCols)] += coeffs[0]*elemVort;
            // i+1,j+1
            vort[cij(iElem+1, jElem+1, nRows, nCols)] += coeffs[1]*elemVort;
            // i,j+1
            vort[cij(iElem, jElem+1, nRows, nCols)] += coeffs[2]*elemVort;
            // i,j
            vort[cij(iElem, jElem, nRows, nCols)] += coeffs[3]*elemVort;
        }
    }
    this->removePicture(0);
    return vort;
}
void PIV::Set::timeSeriesXcorr(char vOrVort)
{
    int wNy = settings.wNy;
    int wNx = settings.wNx;
    this->retrievePicture(0);
    int nCols = this->container[0].frames[0].nx;
    int nRows = this->container[0].frames[0].ny;

    double dx = this->container[0].frames[0].x[0] -
                this->container[0].frames[0].x[1];
    std::vector<double> x = this->container[0].frames[0].x;
    std::vector<double> y = this->container[0].frames[0].y;
    this->removePicture(0);
    
    std::ofstream ofile;
    ofile.open("timeSeriesXcorr.dat");
    ofile << "VARIABLES = \"x\"  \"y\" \"uc\" \"vc\" \"val\"" << std::endl; 
    std::cout << std::endl;
    double xdisp;
    double ydisp;
    double xcorrVal;
    double rms1;
    double rms2;
    std::vector<double> timeSeries1;
    std::vector<double> timeSeries2;
    int nElemsInXcorr = settings.nPics-settings.wNx;
    std::vector<double> xdisps(nRows*nCols, 0.0);
    std::vector<double> ydisps(nRows*nCols, 0.0);
    std::vector<double> xcorrVals(nRows*nCols, 0.0);
    std::vector<std::vector<float>> storage;
    storage = this->preloadData(vOrVort, 0, nElemsInXcorr + settings.maxN);
    for (int iStep=1; iStep<=settings.maxN; ++iStep)
    {
        ofile << "ZONE T=\"n=" << iStep << "\" I=" << nCols 
              << " J=" << nRows << " DATAPACKING=POINT" << std::endl;
        for (int i=0; i<nRows; ++i)
        {
            int imin = ((i-wNy < 0)? 0:i-wNy);
            int imax = ((i+wNy >= nRows)? nRows-1:i+wNy);
            std::cout << "\r" << std::flush;
            std::cout << "n = " << iStep << "  pos = " << (i*1.0)/(nRows-1)
                      << std::flush;
            for (int j=0; j<nCols; ++j)
            {
                //std::cout << "n = " << iStep << "  pos="
                          //<< static_cast<float>(ij(i,j,nRows,nCols))/(nRows*nCols)
                          //<<  std::flush;
                int jmin = ((j-wNx < 0)? 0:j-wNx);
                int jmax = ((j+wNx >= nCols)? nCols-1:j+wNx);

                xdisp = 0;
                ydisp = 0;
                xcorrVal = 0.0;
                int iimax = 0;
                int jjmax = 0;
                double maxXcorrVal = -1.0e10;

                timeSeries1.clear();
                timeSeries1.resize(nElemsInXcorr+1);
                rms1 = 0.0;
                for (int t=0; t<=nElemsInXcorr; ++t)
                {
                    timeSeries1[t] = storage[t][ij(i,j,nRows,nCols)];
                    rms1 += timeSeries1[t]*timeSeries1[t];
                }
                rms1 = std::sqrt(rms1/(nElemsInXcorr+1));
                // DEBUG
                //if (std::abs(rms1) < 1.0e-6)
                //{
                    //std::cerr << __LINE__ << std::endl;
                //}
                int nnCols = imax-imin+1;
                int nnRows = jmax-jmin+1;
                std::vector<double> xcorrValVec(nnRows*nnCols, 0.0);
                int ii = -1;
                for (int di=imin; di<=imax; ++di)
                {
                    ++ii;
                    int jj = -1;
                    for (int dj=jmin; dj<=jmax; ++dj)
                    {
                        ++jj;
                        //timeSeries1 = this->obtainTimeSeries(0, nElemsInXcorr,
                                                 //ij(i,j,nRows,nCols));
                        //timeSeries2 = this->obtainTimeSeries(0+iStep,
                                                //nElemsInXcorr,
                                                //ij(di,dj,nRows,nCols));
                        timeSeries2.clear();
                        timeSeries2.resize(nElemsInXcorr+1);
                        rms2 = 0.0;
                        xcorrVal = 0.0;
                        for (int t=0; t<=nElemsInXcorr; ++t)
                        {
                            timeSeries2[t] = storage[t+iStep][ij(di,dj,nRows,nCols)];

                            rms2 += timeSeries2[t]*timeSeries2[t];
                            //if (std::abs(rms2) < 1.0e-6)
                            //{
                                //std::cerr << storage[t+iStep][ij(di,dj,nRows,nCols)] << std::endl;
                                //std::cerr << __LINE__ << std::endl;
                            //}
                        }
                        rms2 = std::sqrt(rms2/(nElemsInXcorr+1));
                        //if (std::abs(rms2) < 1.0e-6)
                        //{
                            //std::cerr << __LINE__ << std::endl;
                        //}
                        //double div = rms1*rms2*(nElemsInXcorr+1);
                        //if (std::abs(div) < 1.0e-6)
                        //{
                            //std::cerr << __LINE__ << std::endl;
                        //}
                        for (int t=0; t<=nElemsInXcorr; ++t)
                        {
                            xcorrVal += (timeSeries1[t] * timeSeries2[t])/(rms1*rms2);
                        }
                        xcorrValVec[ij(ii,jj,nnRows,nnCols)] = xcorrVal;
                        if (xcorrVal > maxXcorrVal)
                        {
                            maxXcorrVal = xcorrVal;
                            xdisp = dj;
                            ydisp = di;
                            iimax = ii;
                            jjmax = jj;
                        }
                    }
                }
                //// TODO add subpixel interpolation
                //{
                    //int x1 = xdisp;
                    //int ii1 = iimax;
                    ////if (ii1 == 0)
                    ////{
                        ////std::cerr << "ii == 0!!!" << std::endl;
                        ////int temp;
                        ////std::cin >> temp;
                    ////}
                    //// TODO safety for ii1 = 0 in ii2 and ii3
                    //int x2;
                    //int ii2;
                    //if (ii1 == 0)
                    //{
                        //x2 = xdisp + 2;
                        //ii2 = ii1 + 2;
                    //} else
                    //{
                        //x2 = x1 - 1;
                        //ii2 = ii1 - 1;
                    //}

                    //int x3;
                    //int ii3;
                    //if (x1 == nnRows-1)
                    //{
                        //x3 = xdisp - 2;
                        //ii3 = ii1 - 2;
                    //} else
                    //{
                        //x3 = x1 + 1;
                        //ii3 = ii1 + 1;
                    //}

                    //int jj = jjmax;

                    ////int y1 = ydisp;
                    ////int y2 = (y1==0?ydisp+2:ydisp-1);
                    ////int y3 = (y1==nRows?xdisp-2:ydisp+1);

                    //double y1 = (xcorrValVec[ij(ii1, jj, nnRows, nnCols)]);
                    //double y2 = (xcorrValVec[ij(ii2, jj, nnRows, nnCols)]);
                    //double y3 = (xcorrValVec[ij(ii3, jj, nnRows, nnCols)]);

                    //double A = std::log(y3) - std::log(y1);
                    //double B = std::log(y2) - std::log(y1);
                    //xdisp = (B*(x1*x1 - x3*x3) - A*(x1*x1 - x2*x2))/
                                //(2*A*(x2-x1) - 2*B*(x3-x1));
                //}
                xdisps[ij(i,j, nRows, nCols)] = j-xdisp;
                ydisps[ij(i,j, nRows, nCols)] = i-ydisp;
                xcorrVals[ij(i,j, nRows, nCols)] = maxXcorrVal;
            }
        }
        for (int i=0; i<nRows*nCols; ++i)
        {
            ofile << x[i] << "  "
                  << y[i] << "  "
                  << xdisps[i]/iStep << "  "
                  << ydisps[i]/iStep << "  "
                  << xcorrVals[i] << std::endl;
        }
        ofile << std::flush;
    }
    std::cout << std::endl;
    ofile.close();
    storage.clear();
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
