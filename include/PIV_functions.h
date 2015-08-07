#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Set.h"
#include "Picture.h"
#include "getSettings.h"
#include "spectral_filter.h"
#include "PIV_Stitch.h"


void outputAvgVelocityField(PIV::Set& set);
std::vector<std::vector<double>> readAvgVelocityField();
void filterAbsVelocity(PIV::Set& set);
void removeAvg(PIV::Set& set);
void invertSet(PIV::Set& set);
void trimSet(PIV::Set& set);
