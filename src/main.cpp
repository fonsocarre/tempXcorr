#include "Set.h"
#include "Picture.h"
#include "Frame.h"
#include "getSettings.h"
#include "spectral_filter.h"
#include "PIV_Stitch.h"
#include "PIV_Functions.h"

#include <iostream>
/*
 order: stitch -> invert -> trim ->subs avg -> coord
 
 */

int main()
{
    PIV::Set set(settings.inputFile);

    //int nPics = 547;
    //for (int i=0; i<nPics; ++i)
    //{
        //set.retrievePicture(i);
    //}

    //outputAvgVelocityField(set);
    //set.outputAverageField(set.calculateAvgField(), "output_omd_saw_avg.dat");
    //filterAbsVelocity(set);
    set.timeSeriesXcorr('v');

    //set.tecplotOut("completetec_omd_saw.dat");
    //removeAvg(set);


    //invertSet(set);
    //trimSet(set);


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
//    for (int iframe=0; iframe<settings.nPics; ++iframe)
//    {
//        std::cout << "\r" << std::flush;
//        std::cout << "  p: " << iframe << std::flush;
//        set.retrievePicture(iframe);
//        filter.filter(set.container[iframe].frames[0].vx,
//                      out,
//                      set.container[iframe].frames[0].ny,
//                      set.container[iframe].frames[0].nx,
//                      (double)0.4707, settings.padding);
//        filteredSet.insertPicture(set.container[iframe]);
//        filteredSet.unloadPicture(iframe);
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
//        set.removePicture(iframe);
//        filteredSet.removePicture(iframe);
//    }
//    std::cout << std::endl;
//    filteredSet.closeFile();
    //file.close();

    //set.separateTecplotOutput("/Volumes/SSD/Animation/tecplot_cyl_separated/cyl_output", ".dat");
    // 'v' for vx, 'w' for vorticity
    //set.PIV_like_xcorr('v');
    //set.timeSeriesXcorr('v');
    //set.tecplotOut("completetec.dat");
    //set.vorticityTecplotOut("vorticityOut.dat");

    //========================INVERT COORDIN=========================
    //PIV::Set invertedSet = set.copyProperties(settings.inputFile + 
                                                //".inverted.h5");
    //int nPict = invertedSet.nPictures;
    
    //for (int iPict=0; iPict<nPict; ++iPict)
    //{
        //std::cout << __LINE__ << std::flush << std::endl;
        //set.retrievePicture(iPict);
        //std::cout << __LINE__ << std::flush << std::endl;
        //invertedSet.insertPicture(set.container[iPict]);
        //std::cout << __LINE__ << std::flush << std::endl;
        //invertedSet.container[iPict].invertCoordinates();
        //std::cout << __LINE__ << std::flush << std::endl;
        //invertedSet.unloadPicture(iPict);
        //std::cout << __LINE__ << std::flush << std::endl;
        ////invertedSet.removePicture(iPict);
        //std::cout << __LINE__ << std::flush << std::endl;
        //std::cout << "\r" << std::flush;
        //std::cout << "  progress: " << (iPict+1)/(nPict*0.01) << " %"
                    //<< std::flush; 
    //}
    //std::cout << std::endl;
    
    //invertedSet.closeFile();
//    //========================END INVERT=============================
//    // =======================FIND WALL=============================
    //PIV::Frame avgFrame = set.calculateAvgField();
    //int wallLoc = avgFrame.findWallLocation();
    
    //PIV::Set trimmedSet = set.copyProperties(settings.inputFile + 
                                                //".trimmed.h5");
     ////for all picts
    //int nPict = trimmedSet.nPictures;
     ////only for a few
    ////int nPict = settings.nPics;
    
    
    //for (int iPict=0; iPict<nPict; ++iPict)
    //{
        //set.retrievePicture(iPict);
        //PIV::Picture tempPict;
        //set.container[iPict].cutFrames(wallLoc, tempPict);
        //trimmedSet.insertPicture(tempPict);
        //trimmedSet.unloadPicture(iPict);
        ////trimmedSet.removePicture(iPict);
        //std::cout << " " << std::flush;
        //std::cout << "\r" << std::flush;
        //std::cout << "  progress: " << (iPict+1)/(nPict*0.01) << " %";
    //}
    //std::cout << std::endl;
//    // =======================END FIND WALL=========================
    //for (int i=0; i<nPics; ++i)
    //{
        //set.removePicture(i);
    //}
    //

    //========================CORRECT COORD=========================
        //PIV::Set correcSet = set.copyProperties(settings.inputFile + 
                                                    //".correc.h5");
        //int nPict = correcSet.nPictures;
        
        //for (int iPict=0; iPict<nPict; ++iPict)
        //{
            //set.retrievePicture(iPict);
            //correcSet.insertPicture(set.container[iPict]);
            //correcSet.container[iPict].correctFrameCoord(settings.xdisp,
                                                         //settings.ydisp);
            //correcSet.unloadPicture(iPict);
            ////correcSet.removePicture(iPict);
            //set.removePicture(iPict);
            //std::cout << "  progress: " << (iPict+1)/(nPict*0.01) << " %"
                        //<< std::endl; 
        //}
        
        //correcSet.closeFile();
    //========================END CORRECT COORD=============================

     //=======================STITCH=============================
    //int nPics = settings.nPics;
    //std::vector<double> displacements;
    
    //displacements = PIV::getDisplacements(std::vector<int>{0,1},
                                          //set);
    
    //PIV::Set stitchedSet = set.copyProperties();
    //PIV::parallelStitch(displacements, set, stitchedSet, nPics);
        
        //// DEBUG
        //// Tecplot output
    ////    
        //std::cout << "Printing to tecplot file..." << std::endl << std::flush;
        //std::ofstream tecFile;
        //tecFile.open("tecFIle.dat");
        
        //tecFile << "VARIABLES = \"X\" \"Y\" \"VX\" \"VY\"" << std::endl;
        //for (int iZone=0;
             //iZone<settings.nPics;
             //++iZone)
             //{
                 //stitchedSet.retrievePicture(iZone);
                 //tecFile << "ZONE T=\"" << iZone <<
                    //"\", I=" << stitchedSet.container[iZone].frames[0].nx;
                 //tecFile << " ,J=" << stitchedSet.container[iZone].frames[0].ny;
                 //tecFile << " ,SOLUTIONTIME = " << iZone;
                 //tecFile << " ,DATAPACKING = POINT" << std::endl;
                 //for (int j=0; j<stitchedSet.container[iZone].frames[0].nx; ++j)
                 //{
                     //for (int i=0; i<stitchedSet.container[iZone].frames[0].ny; ++i)
                     //{
                         //tecFile <<stitchedSet.container[iZone].frames[0].x
                                    //[j*stitchedSet.container[iZone].frames[0].ny + i]
                                 //<< "  "
                                 //<<stitchedSet.container[iZone].frames[0].y
                                    //[j*stitchedSet.container[iZone].frames[0].ny + i]
                                 //<< "  "
                                 //<<stitchedSet.container[iZone].frames[0].vx
                                    //[j*stitchedSet.container[iZone].frames[0].ny + i]
                                 //<< "  "
                                 //<<stitchedSet.container[iZone].frames[0].vy
                                    //[j*stitchedSet.container[iZone].frames[0].ny + i]
                                 //<< "  "
                                 //<< std::endl;
                     //}
                 //}
                 //stitchedSet.removePicture(iZone);
             //}
        
        //tecFile.close();
    
    //stitchedSet.closeFile();
    // =======================END STITCH=============================
    set.closeFile();
    return 0;
}
