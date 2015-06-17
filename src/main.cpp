#include "Set.h"
#include "Picture.h"
#include "Frame.h"
#include "getSettings.h"

#include <iostream>

int main()
{
    PIV::Set set(settings.inputFile);

    set.timeXcorr();

    set.closeFile();


    return 0;
}
