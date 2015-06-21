#include "Set.h"
#include "Picture.h"
#include "Frame.h"
#include "getSettings.h"

#include <iostream>

int main()
{
    PIV::Set set(settings.inputFile);

    //int nPics = 200;
    //for (int i=0; i<nPics; ++i)
    //{
        //set.retrievePicture(i);
    //}

    //set.separateTecplotOutput("tmp/tecplotOutTest", ".dat");
    //set.timeXcorr();
    set.PIV_like_xcorr();

    //for (int i=0; i<nPics; ++i)
    //{
        //set.removePicture(i);
    //}
    set.closeFile();


    return 0;
}
