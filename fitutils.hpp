#pragma once

#include "TMath.h"
#include "settings.hpp"
#include "utils.hpp"
#include "v2u.hpp"
#include <iostream>

inline void fcn(int &npar, double *gin, double &chisq, double *par, int iflag) {
  double eps0 = par[0];
  double alpha = par[1];
  double kapa = par[2];
  double kapapr = par[3];

  // double eps0 = v/kapa;
  // double alpha = sigm/kapa;
  double kapa2, kapa4, kapa6, kapa8;
  if (FittingV2) {
    kapa2 = kapa;
    kapa4 = kapa;
    kapa6 = kapa;
    kapa8 = kapa;
  } else { // for C2 fitting, it's not just simple kapa
    kapa2 = kapa * kapa;
    kapa4 = -TMath::Power(kapa, 4.);
    kapa6 = 4 * TMath::Power(kapa, 6);
    kapa8 = -33 * TMath::Power(kapa, 8);
  };

  double v22ep = kapa2 * eps22EP(alpha, eps0) +
                 kapapr * kapa2 * eps22EP(alpha, eps0) * eps22EP(alpha, eps0) *
                     eps22EP(alpha, eps0);
  double v24ep = kapa4 * eps24EP(alpha, eps0) +
                 kapapr * kapa4 * eps24EP(alpha, eps0) * eps24EP(alpha, eps0) *
                     eps24EP(alpha, eps0);
  double v26ep = kapa6 * eps26EP(alpha, eps0) +
                 kapapr * kapa6 * eps26EP(alpha, eps0) * eps26EP(alpha, eps0) *
                     eps26EP(alpha, eps0);
  double v28ep = kapa8 * eps28EP(alpha, eps0) +
                 kapapr * kapa8 * eps28EP(alpha, eps0) * eps28EP(alpha, eps0) *
                     eps28EP(alpha, eps0);

  double d1 = (v22ep - v2i) / v2ei;
  double d2 = (v24ep - v4i) / v4ei;
  double d3 = (v26ep - v6i) / v6ei;
  double d4 = (v28ep - v8i) / v8ei;

  chisq = d1 * d1 + d2 * d2 + d3 * d3 + d4 * d4;
}

inline void meanvEP246(double v2, double v2e, double v4, double v4e, double v6,
                       double v6e, double v8, double v8e, double &eps0,
                       double &eps0Err, double &alpha, double &alphaErr,
                       double &kapa, double &kapaErr, double &kapapr,
                       double &kapaprErr, double KappaInp) {

  v2i = v2;
  v2ei = v2e;
  v4i = v4;
  v4ei = v4e;
  v6i = v6;
  v6ei = v6e;
  v8i = v8;
  v8ei = v8e;

  std::cout << "v2=" << v2 << "  v2e=" << v2e << std::endl;
  std::cout << "v4=" << v4 << "  v4e=" << v4e << std::endl;
  std::cout << "v6=" << v6 << "  v6e=" << v6e << std::endl;
  std::cout << "v8=" << v8 << "  v8e=" << v8e << std::endl;
  //     minuit->mnexcm("MIGRAD", arglist ,1,ierflg);
  if (v2 > v4) {
    std::cout << "STOP" << endl;
  }

  // initialize TMinuit with a maximum of 3 params
  TMinuit *minuit = new TMinuit(4);
  minuit->SetFCN(fcn);

  double arglist[10];
  int ierflg = 0;
  // Now ready for minimization step
  arglist[0] = 2000; // 2000 standard
  arglist[1] = 1;

  double step[4] = {0.01, 0.01, 0.01, 0.01};
  double lowLim[4] = {0, -500, 0, 0};
  double upLim[4] = {1, 1500, 1, 1}; // kpr can be 1

  double vstart[4] = {GetE0(), GetAlpha(), 0.25,
                      0.1}; // 0.4 initially for k_n (3rd) 0.3 for default value
  if (FixKappaPrime)
    vstart[3] = 0;

  minuit->mnparm(0, "eps0", vstart[0], step[0], lowLim[0], upLim[0], ierflg);
  minuit->mnparm(1, "alpha", vstart[1], step[1], lowLim[1], upLim[1], ierflg);
  minuit->mnparm(2, "kapa", vstart[2], step[2], lowLim[2], upLim[2], ierflg);
  minuit->mnparm(3, "kapapr", vstart[3], step[3], lowLim[3], upLim[3], ierflg);
  if (FixKappaPrime)
    minuit->FixParameter(3);
  minuit->SetPrintLevel(0);

  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  int nfits = 0;
  TString status = minuit->fCstatu.Data();

  while ((!status.Contains("CONVERGED")) && (nfits < 10)) {
    double scalef = (nfits + 1) / 10;
    double vstartn[4] = {GetE0() * (1 - scalef - 0.5),
                         GetAlpha() * (1 - scalef + 0.2),
                         0.4 * (1 - scalef - 0.5), 0};
    minuit->mnparm(0, "eps0", vstartn[0], step[0], lowLim[0], upLim[0], ierflg);
    minuit->mnparm(1, "alpha", vstartn[1], step[1], lowLim[1], upLim[1],
                   ierflg);
    minuit->mnparm(2, "kapa", vstartn[2], step[2], lowLim[2], upLim[2], ierflg);
    minuit->mnparm(3, "kapapr", vstartn[3], step[3], lowLim[3], upLim[3],
                   ierflg);
    if (FixKappaPrime)
      minuit->FixParameter(3);

    minuit->SetPrintLevel(0);

    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    status = minuit->fCstatu.Data();
    nfits++;
  }

  //  minuit->Migrad();
  minuit->GetParameter(0, eps0, eps0Err);
  minuit->GetParameter(1, alpha, alphaErr);
  minuit->GetParameter(2, kapa, kapaErr);
  minuit->GetParameter(3, kapapr, kapaprErr);

  std::cout << " eps0 = " << eps0 << " +/- " << eps0Err << std::endl;
  std::cout << " alpha = " << alpha << " +/- " << alphaErr << std::endl;
  std::cout << " kapa = " << kapa << " +/- " << kapaErr << std::endl;
  std::cout << " kapapr = " << kapapr << " +/- " << kapaprErr << std::endl;

  // Print results
  double amin, edm, errdef;
  int nvpar, nparx, icstat;
}
