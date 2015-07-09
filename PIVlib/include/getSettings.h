/*
 © Alfonso del Carre 2015
 Imperial College London
 ===================================
 getSettings.h
 Settings file parser functions
 ===================================
 */
#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <exception>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

struct Settings {
    double maxXdisp;
    int nPicsForStitch;
    int nPics;
    int nThreads;
    std::string inputFile;
    double stitchXcorrOffset;
    double xdisp;
    double ydisp;

    int tempXcorrMinPict;
    int tempXcorrMaxPict;
    int tempXcorrMaxX;
    int tempXcorrn;

    int maxN;
    int maxR;
    int deltaR;
    double minX;
    int nFramesRxxNorm;
    int minY;
    int maxY;
    int wNx;
    int wNy;

    double cutLength;
    int padding;
};

Settings getSettings(const std::string settingsFile);

extern Settings settings;
