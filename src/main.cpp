#include "Set.h"
#include "Picture.h"
#include "Frame.h"
#include "getSettings.h"

#include <iostream>

int main()
{
    PIV::Set set(settings.inputFile);

    //int nPics = 5472;
    //for (int i=0; i<nPics; ++i)
    //{
        //set.retrievePicture(i);
    //}

    //set.separateTecplotOutput("tmp/tecplotOutSaw", ".dat");
    //set.timeXcorr();
    //set.PIV_like_xcorr();
    set.tecplotOut("completetec_saw.dat");

    //for (int i=0; i<nPics; ++i)
    //{
        //set.removePicture(i);
    //}
    set.closeFile();


    return 0;
}
