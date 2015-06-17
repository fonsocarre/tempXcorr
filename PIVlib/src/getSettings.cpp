//
//  getSettings.cpp
//  Stitch
//
//  Created by Alfonso del Carre on 01/06/2015.
//  Copyright (c) 2015 Alfonso del Carre. All rights reserved.
//

#include "getSettings.h"

Settings settings = getSettings(
    "../settings.ini");


Settings getSettings(const std::string settingsFile)
{
    std::cout << "===========================================" << std::endl;
    std::cout << "    Reading settings file: " << settingsFile << std::endl;
    std::cout << "===========================================" << std::endl;
    
    Settings output;
    boost::property_tree::ptree pt;
    try
    {
        boost::property_tree::ini_parser::read_ini(settingsFile, pt);
    }
    catch (std::exception &ex)
    {
        std::cerr << "Error reading settings file!!!" << std::endl;
    }
    
        output.inputFile = boost::lexical_cast<std::string>(pt.get<std::string>
                                                    ("Settings.inputFile"));
    //    output.circOutFile = boost::lexical_cast<std::string>(pt.get<std::string>
    //                                                ("Files.CirculationOutputFile"));
        output.nPicsForStitch = boost::lexical_cast<int>(pt.get<std::string>
                                                     ("Settings.nPicsForStitch"));
        output.nPics = boost::lexical_cast<int>(pt.get<std::string>
                                                     ("Settings.nPics"));
        output.nThreads = boost::lexical_cast<int>(pt.get<std::string>
                                                ("Settings.nThreads"));
        output.maxXdisp = boost::lexical_cast<double>(pt.get<std::string>
                                                        ("Settings.maxXdisp"));
//        output.stitchXcorrOffset = boost::lexical_cast<double>(pt.get<std::string>
//                                                    ("DVM.stitchXcorrOffset"));
        output.xdisp = boost::lexical_cast<double>(pt.get<std::string>
                                                    ("Settings.xdisp"));
        output.ydisp = boost::lexical_cast<double>(pt.get<std::string>
                                                ("Settings.ydisp"));
    //    output.separationAngle = boost::lexical_cast<double>(pt.get<std::string>
    //                                                ("DVM.separationAngle"));
    //    output.separationAngle *= (constants::PI/180.);
        output.tempXcorrMinPict = boost::lexical_cast<int>(pt.get<std::string>
                                                    ("Settings.tempXcorrMinPict"));
        output.tempXcorrMaxPict = boost::lexical_cast<int>(pt.get<std::string>
                                                    ("Settings.tempXcorrMaxPict"));
        output.tempXcorrMaxX = boost::lexical_cast<int>(pt.get<std::string>
                                                    ("Settings.tempXcorrMaxX"));
        output.tempXcorrn = boost::lexical_cast<int>(pt.get<std::string>
                                                    ("Settings.tempXcorrn"));
        output.maxN = boost::lexical_cast<int>(pt.get<std::string>
                                                    ("Settings.maxN"));
        output.maxR = boost::lexical_cast<int>(pt.get<std::string>
                                                    ("Settings.maxR"));
        output.deltaR = boost::lexical_cast<int>(pt.get<std::string>
                                                    ("Settings.deltaR"));
        output.minX = boost::lexical_cast<double>(pt.get<std::string>
                                                    ("Settings.minX"));
    //    output.dxsdys.resize(constants::DIM);
    //    output.dxsdys[0] = boost::lexical_cast<double>(pt.get<std::string>
    //                                                ("DVM.dxs"));
    //    output.dxsdys[1] = boost::lexical_cast<double>(pt.get<std::string>
    //                                                ("DVM.dys"));
    //    output.xNear = boost::lexical_cast<double>(pt.get<std::string>
    //                                                   ("DVM.xNear"));
    //    output.xFar = boost::lexical_cast<double>(pt.get<std::string>
    //                                                   ("DVM.xFar"));
    //
    //    std::cout << "    Settings file read." << std::endl;
    //    std::cout << "*********" << std::endl << std::endl;
    //    std::cout << "  Output File: " << output.outputFile << std::endl;
    //    std::cout << "  Physical time step of the simulation: "
    //                << output.dt*output.a/output.u << std::endl;
    //    std::cout << "  Non-dimensional time step of the simulation: "
    //                << output.dt << std::endl;
    return output;
}

