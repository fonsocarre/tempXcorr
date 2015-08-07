#include "PIV_functions.h"

void outputAvgVelocityField(PIV::Set& set)
{
    std::cout << "Getting avg velocity field..." << std::endl;
    PIV::Frame avgFrame = set.calculateAvgField();
    std::ofstream ofile;
    ofile.open(settings.avgVelFile);
    ofile << avgFrame.n << std::endl;
    for (int i=0; i<avgFrame.n; ++i)
    {
        ofile << avgFrame.vx[i] << std::endl;
    }
    for (int i=0; i<avgFrame.n; ++i)
    {
        ofile << avgFrame.vy[i] << std::endl; 
    }
    ofile.close();
}

std::vector<std::vector<double>> readAvgVelocityField()
{
    std::vector<std::vector<double>> avgField;
    int n;
    std::ifstream ifile;
    ifile.open(settings.avgVelFile);
    ifile >> n;
    avgField.resize(2);
    avgField[0].resize(n);
    avgField[1].resize(n);
    for (int i=0; i<n; ++i)
    {
        ifile >> avgField[0][i];
    }
    for (int i=0; i<n; ++i)
    {
        ifile >> avgField[1][i];
    }
    ifile.close();
    return avgField;
}

void filterAbsVelocity(PIV::Set& set)
{
    // Input has to be the fluc component!
    PIV::Set filteredSet;
    filteredSet = set.copyProperties(settings.inputFile + ".abs.filtered.h5");
    set.retrievePicture(0);
    int n = set.container[0].frames[0].nx * 
            set.container[0].frames[0].ny;
    std::vector<double> out(set.container[0].frames[0].n, 0.0);
    spectral_filter filter;
    filter.initialise(settings.cutLength);
    set.removePicture(0);

    std::vector<std::vector<double>> avgField = readAvgVelocityField();
    std::vector<double> tvx;
    std::vector<double> tvy;
    for (int iframe=0; iframe<settings.nPics; ++iframe)
    {
        tvx.clear();
        tvy.clear();
        std::cout << "\r" << std::flush;
        std::cout << "  p: " << iframe << std::flush;
        set.retrievePicture(iframe);
        filter.filter(set.container[iframe].frames[0].vx,
                      tvx,
                      set.container[iframe].frames[0].ny,
                      set.container[iframe].frames[0].nx,
                      (double)0.4707, settings.padding);
        filter.filter(set.container[iframe].frames[0].vy,
                      tvy,
                      set.container[iframe].frames[0].ny,
                      set.container[iframe].frames[0].nx,
                      (double)0.4707, settings.padding);



        //std::cout << __LINE__ << std::endl << std::flush;
        filteredSet.insertPicture(set.container[iframe]);
        //std::cout << __LINE__ << std::endl << std::flush;
        for (int i=0; i<n; ++i)
        {
            //std::cout << tvx[i] << ", " << tvy[i] << std::endl;
            //filteredSet.container[iframe].frames[0].vx[i] = tvx[i] + avgField[0][i];
            //filteredSet.container[iframe].frames[0].vy[i] = tvy[i] + avgField[1][i];
            filteredSet.container[iframe].frames[0].vx[i] = tvx[i] + avgField[0][i];
            filteredSet.container[iframe].frames[0].vy[i] = tvy[i] + avgField[1][i];
        }
        //std::cout << __LINE__ << std::endl << std::flush;
        if (iframe == 0)
        {
            std::ofstream file;
            file.open("filtertest.dat");
            file << "VARIABLES = \"x\"  \"y\" \"avgU\" \"avgV\" \"u'\" \"v'\" \"u\" \"v\" " << std::endl;

            file << "ZONE T=\"0 \" I=" 
                <<  filteredSet.container[0].frames[0].nx 
                << " J=" << filteredSet.container[0].frames[0].ny
                << " DATAPACKING=POINT" << std::endl;
            for (int i=0; i<n; ++i)
            {
                file << filteredSet.container[0].frames[0].x[i] << "  ";
                file << filteredSet.container[0].frames[0].y[i] << "  ";
                file << avgField[0][i] << "  ";
                file << avgField[1][i] << "  ";
                file << tvx[i] << "  ";
                file << tvy[i] << "  ";
                file << filteredSet.container[0].frames[0].vx[i] << "  ";
                file << filteredSet.container[0].frames[0].vy[i] << "  ";
                file << std::endl;
            }
        }
        //std::cout << __LINE__ << std::endl << std::flush;
        filteredSet.unloadPicture(iframe);

        //std::cout << __LINE__ << std::endl << std::flush;
        set.removePicture(iframe);
        //std::cout << __LINE__ << std::endl;
        //std::cout << filteredSet.loadedPictures[0] << std::endl;
        //filteredSet.removePicture(iframe);
    }
    std::cout << std::endl;
    filteredSet.closeFile();
}


void removeAvg(PIV::Set& set)
{
    PIV::Set flucSet = set.copyProperties(settings.inputFile + ".fluc.h5");
    set.retrievePicture(0);
    int n = set.container[0].frames[0].nx * 
            set.container[0].frames[0].ny;
    PIV::Frame avgFrame = set.calculateAvgField();
    for (int iPict=0; iPict<settings.nPics; ++iPict)
    {
        set.retrievePicture(iPict);

        for (int i=0; i<n; ++i)
        {
            set.container[iPict].frames[0].vx[i] -= avgFrame.vx[i];
            set.container[iPict].frames[0].vy[i] -= avgFrame.vy[i];
        }

        flucSet.insertPicture(set.container[iPict]);
        set.removePicture(iPict);
        flucSet.unloadPicture(iPict);
    }
    flucSet.closeFile();
}

void invertSet(PIV::Set& set)
{
    //========================INVERT COORDIN=========================
    PIV::Set invertedSet = set.copyProperties(settings.inputFile + 
                                                ".inverted.h5");
    int nPict = invertedSet.nPictures;
    
    for (int iPict=0; iPict<settings.nPics; ++iPict)
    {
        //std::cout << __LINE__ << std::flush << std::endl;
        set.retrievePicture(iPict);
        //std::cout << __LINE__ << std::flush << std::endl;
        invertedSet.insertPicture(set.container[iPict]);
        //std::cout << __LINE__ << std::flush << std::endl;
        invertedSet.container[iPict].invertCoordinates();
        //std::cout << __LINE__ << std::flush << std::endl;
        invertedSet.unloadPicture(iPict);
        //std::cout << __LINE__ << std::flush << std::endl;
        //invertedSet.removePicture(iPict);
        //std::cout << __LINE__ << std::flush << std::endl;
        std::cout << "\r" << std::flush;
        std::cout << "  progress: " << (iPict+1)/(nPict*0.01) << " %"
                    << std::flush; 
    }
    std::cout << std::endl;
    
    invertedSet.closeFile();
}
void trimSet(PIV::Set& set)
{
    PIV::Frame avgFrame = set.calculateAvgField();
    int wallLoc = avgFrame.findWallLocation();
    
    PIV::Set trimmedSet = set.copyProperties(settings.inputFile + 
                                                ".trimmed.h5");
     //for all picts
    //int nPict = trimmedSet.nPictures;
     //only for a few
    int nPict = settings.nPics;
    
    
    for (int iPict=0; iPict<nPict; ++iPict)
    {
        set.retrievePicture(iPict);
        PIV::Picture tempPict;
        set.container[iPict].cutFrames(wallLoc, tempPict);
        trimmedSet.insertPicture(tempPict);
        trimmedSet.unloadPicture(iPict);
        //trimmedSet.removePicture(iPict);
        std::cout << " " << std::flush;
        std::cout << "\r" << std::flush;
        std::cout << "  progress: " << (iPict+1)/(nPict*0.01) << " %";
    }
    std::cout << std::endl;
}
