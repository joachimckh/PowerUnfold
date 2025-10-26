#pragma once
#include "TMath.h"
#include "settings.hpp"



double eps22EP(double alpha, double eps0) {
  double hypergf1 = ROOT::Math::hyperg(1.5, 1., alpha + 2., eps0 * eps0);

  // Double_t hypergf1 = gsl_sf_hyperg_2F1(1.5, 1., alpha+2., eps0*eps0);

  double qc2 = 1. - alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1;
  double eps22 = FittingV2 ? (TMath::Sqrt(qc2)) : (qc2);

  return eps22;
}

double eps24EP(double alpha, double eps0) {

  double hypergf1 = ROOT::Math::hyperg(1.5, 1., alpha + 2., eps0 * eps0);
  double hypergf2 = ROOT::Math::hyperg(2.5, 2., alpha + 3., eps0 * eps0);

  // Double_t hypergf1 = gsl_sf_hyperg_2F1(1.5, 1., alpha+2., eps0*eps0);
  // Double_t hypergf2 = gsl_sf_hyperg_2F1(2.5, 2., alpha+3., eps0*eps0);

  double qc4 =
      1. - 2. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 +
      2. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 * alpha /
          (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 -
      alpha / (alpha + 2.) * (1. - eps0 * eps0) * (1. - eps0 * eps0) * hypergf2;
  double eps24 = FittingV2 ? (TMath::Power(qc4, 1. / 4.)) : qc4;

  return eps24;
}

double eps26EP(double alpha, double eps0) {

  double hypergf1 = ROOT::Math::hyperg(1.5, 1., alpha + 2., eps0 * eps0);
  double hypergf2 = ROOT::Math::hyperg(2.5, 2., alpha + 3., eps0 * eps0);
  double hypergf3 = ROOT::Math::hyperg(3.5, 3., alpha + 4., eps0 * eps0);

  /*
  Double_t hypergf1 = gsl_sf_hyperg_2F1(1.5, 1., alpha+2., eps0*eps0);
  Double_t hypergf2 = gsl_sf_hyperg_2F1(2.5, 2., alpha+3., eps0*eps0);
  Double_t hypergf3 = gsl_sf_hyperg_2F1(3.5, 3., alpha+4., eps0*eps0);
  */

  double qc6 = 1. +
               9. / 2. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                   alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 -
               3. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                   alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                   alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 +
               3. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                   (3. / 4. * alpha / (alpha + 2.) * (1. - eps0 * eps0) *
                        (1. - eps0 * eps0) * hypergf2 -
                    1.) -
               3. / 2. * alpha / (alpha + 2.) * (1. - eps0 * eps0) *
                   (1. - eps0 * eps0) * hypergf2 -
               1. / 4. * alpha / (alpha + 3.) * (1. - eps0 * eps0) *
                   (1. - eps0 * eps0) * (1. - eps0 * eps0) * hypergf3;

  double eps26 = FittingV2 ? (TMath::Power(qc6, 1. / 6.)) : qc6;

  return eps26;
}

double eps28EP(double alpha, double eps0) {

  double hypergf1 = ROOT::Math::hyperg(1.5, 1., alpha + 2., eps0 * eps0);
  double hypergf2 = ROOT::Math::hyperg(2.5, 2., alpha + 3., eps0 * eps0);
  double hypergf3 = ROOT::Math::hyperg(3.5, 3., alpha + 4., eps0 * eps0);
  double hypergf4 = ROOT::Math::hyperg(4.5, 4., alpha + 5., eps0 * eps0);

  /*
  Double_t hypergf1 = gsl_sf_hyperg_2F1(1.5, 1., alpha+2., eps0*eps0);
  Double_t hypergf2 = gsl_sf_hyperg_2F1(2.5, 2., alpha+3., eps0*eps0);
  Double_t hypergf3 = gsl_sf_hyperg_2F1(3.5, 3., alpha+4., eps0*eps0);
  Double_t hypergf4 = gsl_sf_hyperg_2F1(4.5, 4., alpha+5., eps0*eps0);
  */

  double term2 = -288. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                 alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 * alpha /
                 (alpha + 1.) * (1. - eps0 * eps0) * hypergf1;
  double term3 = 144. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                 alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 * alpha /
                 (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 * alpha /
                 (alpha + 1.) * (1. - eps0 * eps0) * hypergf1;
  double term4 = -66. * alpha / (alpha + 2.) * (1. - eps0 * eps0) *
                 (1. - eps0 * eps0) * hypergf2;
  double term5 = 18. * alpha / (alpha + 2.) * (1. - eps0 * eps0) *
                 (1. - eps0 * eps0) * hypergf2 * alpha / (alpha + 2.) *
                 (1. - eps0 * eps0) * (1. - eps0 * eps0) * hypergf2;
  double term6 = -24. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                 alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                 (-11. + 6. * alpha / (alpha + 2.) * (1. - eps0 * eps0) *
                             (1. - eps0 * eps0) * hypergf2);
  double term7 = -12. * alpha / (alpha + 3.) * (1. - eps0 * eps0) *
                 (1. - eps0 * eps0) * (1. - eps0 * eps0) * hypergf3;
  double term8 = 4. * alpha / (alpha + 1.) * (1. - eps0 * eps0) * hypergf1 *
                 (-33. +
                  42. * alpha / (alpha + 2.) * (1. - eps0 * eps0) *
                      (1. - eps0 * eps0) * hypergf2 +
                  4. * alpha / (alpha + 3.) * (1. - eps0 * eps0) *
                      (1. - eps0 * eps0) * (1. - eps0 * eps0) * hypergf3);
  double term9 = -alpha / (alpha + 4.) * (1. - eps0 * eps0) *
                 (1. - eps0 * eps0) * (1. - eps0 * eps0) * (1. - eps0 * eps0) *
                 hypergf4;

  double qc8 =
      (33. + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9) /
      33.;

  double eps28 = FittingV2 ? (TMath::Power(qc8, 1. / 8.)) : qc8;

  return eps28;
}
