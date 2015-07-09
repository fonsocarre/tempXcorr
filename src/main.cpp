#include "Set.h"
#include "Picture.h"
#include "Frame.h"
#include "getSettings.h"
#include "spectral_filter.h"

#include <iostream>

int main()
{
    PIV::Set set(settings.inputFile);

    //int nPics = 547;
    //for (int i=0; i<nPics; ++i)
    //{
        //set.retrievePicture(i);
    //}



    // FILTERING
    //PIV::Set filteredSet;
    //filteredSet = set.copyProperties(settings.inputFile + ".filtered.h5");
    //set.retrievePicture(0);
    //int n = set.container[0].frames[0].nx * 
            //set.container[0].frames[0].ny;
    //std::vector<double> out(set.container[0].frames[0].n, 0.0);
    //spectral_filter filter;
    //filter.initialise(settings.cutLength);
    
    //std::ofstream file;
    //file.open("filtertest.dat");
    //file << "VARIABLES = \"x\"  \"y\" \"u\" \"v\" \"fu\"  " << std::endl;
    //std::cout << std::endl;
    //for (int iframe=0; iframe<settings.nPics; ++iframe)
    //{
        //std::cout << "\r" << std::flush;
        //std::cout << "  p: " << iframe << std::flush;
        //set.retrievePicture(iframe);
        //filter.filter(set.container[iframe].frames[0].vx,
                      //out,
                      //set.container[iframe].frames[0].ny,
                      //set.container[iframe].frames[0].nx,
                      //(double)0.4707, settings.padding);
        //filteredSet.insertPicture(set.container[iframe]);
        //filteredSet.unloadPicture(iframe);
        //file << "ZONE T=\"" << iframe << "\" I=" 
                //<< set.container[iframe].frames[0].nx 
                //<< " J=" << set.container[iframe].frames[0].ny
              //<< " DATAPACKING=POINT" << std::endl;
        //for (int i=0; i<n; ++i)
        //{
            //file << set.container[iframe].frames[0].x[i] << "  ";
            //file << set.container[iframe].frames[0].y[i] << "  ";
            //file << set.container[iframe].frames[0].vx[i] << "  ";
            //file << set.container[iframe].frames[0].vy[i] << "  ";
            //file << out[i] << std::endl;
        //}
        //set.removePicture(iframe);
        //filteredSet.removePicture(iframe);
    //}
    //std::cout << std::endl;
    //filteredSet.closeFile();
    //file.close();

    //set.separateTecplotOutput("tmp/tecplotOutSaw", ".dat");
    //set.timeXcorr();
    set.PIV_like_xcorr();
    //set.tecplotOut("completetec_saw.dat");

    //for (int i=0; i<nPics; ++i)
    //{
        //set.removePicture(i);
    //}

    set.closeFile();
    return 0;
}
