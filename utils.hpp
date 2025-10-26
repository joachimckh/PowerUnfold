#pragma once

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMinuit.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

#include "settings.hpp"

inline double GetE0() {
  // return 0.031+0.008*gCent;
  return 0.015 + 0.008 * gCent; // Ja IV 0.15
};
inline double GetAlpha() {
  // return 5.372+TMath::Exp(-0.052*gCent);
  return 60.0 + TMath::Exp(-0.052 * gCent); // Ja IV 60.0
};
