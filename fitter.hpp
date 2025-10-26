#pragma once

#include "fitutils.hpp"
#include "settings.hpp"
#include <vector>

R__LOAD_LIBRARY(libMathMore);

void getParamEP(int optDrawEP = 1)
{

  gSystem->Load("libMathCore.so");
  gSystem->Load("libMathMore.so"); 
  gSystem->AddIncludePath("-I$ROOTSYS/include"); 


  ifstream inFile;
  inFile.open(Form("%sout.txt",FittingV2?"v2":"c2"));


  double cent, v2sp, v2sper, v24, v24er, v26, v26er, v28, v28er;
  vector<double> vcent, vv22, vv22er, vv24, vv24er, vv26, vv26er, vv28, vv28er;
  
  while (!inFile.eof()) {
    inFile>>cent>>v2sp>>v2sper>>v24>>v24er>>v26>>v26er>>v28>>v28er;
    vcent.push_back(cent);
    vv22.push_back(v2sp);
    vv22er.push_back(v2sper);
    vv24.push_back(v24);
    vv24er.push_back(v24er);
    vv26.push_back(v26);
    vv26er.push_back(v26er);
    vv28.push_back(v28);
    vv28er.push_back(v28er);
  }
  gCent = cent; //storing for initial fit parameters

  inFile.close();

    
    auto nPoints = vcent.size();
    
    std::cout<<"Number of points = "<<nPoints<<std::endl;
  
    TGraphErrors* grEps0 = new TGraphErrors(nPoints-1);
    grEps0->SetName("grEps0");
    TGraphErrors* grAlpha = new TGraphErrors(nPoints-1);
    grAlpha->SetName("grAlpha");
    TGraphErrors* grKappa = new TGraphErrors(nPoints-1);
    grKappa->SetName("grKappa");
    TGraphErrors* grKappaPr = new TGraphErrors(nPoints-1);
    grKappaPr->SetName("grKappaPr");
   
    
    
    TGraphErrors* grv22 = new TGraphErrors(nPoints-1);
    TGraphErrors* grv24 = new TGraphErrors(nPoints-1);
    TGraphErrors* grv26 = new TGraphErrors(nPoints-1);
    TGraphErrors* grv28 = new TGraphErrors(nPoints-1);
        
    TGraphErrors* grv22EP = new TGraphErrors(nPoints-1);
    TGraphErrors* grv24EP = new TGraphErrors(nPoints-1);
    TGraphErrors* grv26EP = new TGraphErrors(nPoints-1);
    TGraphErrors* grv28EP = new TGraphErrors(nPoints-1);
    
    
    TGraphErrors* grv22Rat = new TGraphErrors(nPoints-1);
    grv22Rat->SetName("grv22Rat");
    TGraphErrors* grv24Rat = new TGraphErrors(nPoints-1);
    grv24Rat->SetName("grv24Rat");
    TGraphErrors* grv26Rat = new TGraphErrors(nPoints-1);
    grv26Rat->SetName("grv26Rat");
    TGraphErrors* grv28Rat = new TGraphErrors(nPoints-1);
    grv28Rat->SetName("grv28Rat");

    //testing fix kap
    double kapInp[8] = {0.2033558221352963,0.18424869966769836,0.13678503767135297,0.18338891964171322,0.20838626240128438,0.22669365529172347,0.2228178854885996,0.19562231905346555};
    ///
    for (int i = 0; i < nPoints-1; i++){
      std::cout<<"     "<<std::endl;
      std::cout<<"******"<<std::endl;
      std::cout<<"!!!!!!!!! Fitting centrality "<<vcent[i]<<std::endl;
      std::cout<<"******"<<std::endl;

      double Eps0, Eps0Err, Alpha, AlphaErr, Kappa, KappaErr, Kappapr, KappaprErr;
        meanvEP246(vv22[i], vv22er[i], vv24[i], vv24er[i], vv26[i], vv26er[i], vv28[i], vv28er[i], Eps0, Eps0Err, Alpha, AlphaErr, Kappa, KappaErr, Kappapr, KappaprErr,kapInp[i]);


        grEps0->SetPoint(i, vcent[i], Eps0);
        grEps0->SetPointError(i, 0, Eps0Err);

        grAlpha->SetPoint(i, vcent[i], Alpha);
        grAlpha->SetPointError(i, 0, AlphaErr);

        grKappa->SetPoint(i, vcent[i], Kappa);
        grKappa->SetPointError(i, 0, KappaErr);
        
        grKappaPr->SetPoint(i, vcent[i], Kappapr);
        grKappaPr->SetPointError(i, 0, KappaprErr);

        
        
        if (optDrawEP){

	  double kapa2, kapa4, kapa6, kapa8;
	  if(FittingV2) { 
	    kapa2 = Kappa;
	    kapa4 = Kappa;
	    kapa6 = Kappa;
	    kapa8 = Kappa;
	  } else { //for C2 fitting, it's not just simple kapa
	    kapa2 = Kappa*Kappa;
	    kapa4 = -TMath::Power(Kappa, 4.);
	    kapa6 = 4*TMath::Power(Kappa, 6);
	    kapa8 = -33*TMath::Power(Kappa, 8);
	  };
        
            double V22ep = kapa2*eps22EP(Alpha, Eps0) + kapa2*Kappapr*eps22EP(Alpha, Eps0)*eps22EP(Alpha, Eps0)*eps22EP(Alpha, Eps0);
            double V24ep = kapa4*eps24EP(Alpha, Eps0) + kapa4*Kappapr*eps24EP(Alpha, Eps0)*eps24EP(Alpha, Eps0)*eps24EP(Alpha, Eps0);
            double V26ep = kapa6*eps26EP(Alpha, Eps0) + kapa6*Kappapr*eps26EP(Alpha, Eps0)*eps26EP(Alpha, Eps0)*eps26EP(Alpha, Eps0);
            double V28ep = kapa8*eps28EP(Alpha, Eps0) + kapa8*Kappapr*eps28EP(Alpha, Eps0)*eps28EP(Alpha, Eps0)*eps28EP(Alpha, Eps0);
        
        
            double V22epErr = KappaErr*eps22EP(AlphaErr, Eps0Err) + KappaErr*KappaprErr*eps22EP(AlphaErr, Eps0Err)*eps22EP(AlphaErr, Eps0Err)*eps22EP(AlphaErr, Eps0Err);
            double V24epErr = KappaErr*eps24EP(AlphaErr, Eps0Err) + KappaErr*KappaprErr*eps24EP(AlphaErr, Eps0Err)*eps24EP(AlphaErr, Eps0Err)*eps24EP(AlphaErr, Eps0Err);
            double V26epErr = KappaErr*eps26EP(AlphaErr, Eps0Err) + KappaErr*KappaprErr*eps26EP(AlphaErr, Eps0Err)*eps26EP(AlphaErr, Eps0Err)*eps26EP(AlphaErr, Eps0Err);
            double V28epErr = KappaErr*eps28EP(AlphaErr, Eps0Err) + KappaErr*KappaprErr*eps28EP(AlphaErr, Eps0Err)*eps28EP(AlphaErr, Eps0Err)*eps28EP(AlphaErr, Eps0Err);
        
      
            grv22->SetPoint(i, vcent[i], vv22[i]);
            grv22->SetPointError(i, 0, vv22er[i]);
      
            grv24->SetPoint(i, vcent[i], vv24[i]);
            grv24->SetPointError(i, 0, vv24er[i]);
      
            grv26->SetPoint(i, vcent[i], vv26[i]);
            grv26->SetPointError(i, 0, vv26er[i]);
        
            grv28->SetPoint(i, vcent[i], vv28[i]);
            grv28->SetPointError(i, 0, vv28er[i]);
        
      
            grv22EP->SetPoint(i, vcent[i], V22ep);
            grv22EP->SetPointError(i, 0, V22epErr);
      
            grv24EP->SetPoint(i, vcent[i], V24ep);
            grv24EP->SetPointError(i, 0, V24epErr);
      
            grv26EP->SetPoint(i, vcent[i], V26ep);
            grv26EP->SetPointError(i, 0, V26epErr);
        
            grv28EP->SetPoint(i, vcent[i], V28ep);
            grv28EP->SetPointError(i, 0, V28epErr);
        
        
        
            double ratV22 = V22ep/vv22[i];
            double errRatV22 = TMath::Sqrt(TMath::Abs(vv22er[i]*vv22er[i] - V22epErr*V22epErr));
        
            double ratV24 = V24ep/vv24[i];
            double errRatV24 = TMath::Sqrt(TMath::Abs(vv24er[i]*vv24er[i] - V24epErr*V24epErr));
        
            double ratV26 = V26ep/vv26[i];
            double errRatV26 = TMath::Sqrt(TMath::Abs(vv26er[i]*vv26er[i] - V26epErr*V26epErr));
        
            double ratV28 = V28ep/vv28[i];
            double errRatV28 = TMath::Sqrt(TMath::Abs(vv28er[i]*vv28er[i] - V28epErr*V28epErr));
        
        
            grv22Rat->SetPoint(i, vcent[i], ratV22);
            grv22Rat->SetPointError(i, 0, errRatV22);
        
            grv24Rat->SetPoint(i, vcent[i], ratV24);
            grv24Rat->SetPointError(i, 0, errRatV24);
        
            grv26Rat->SetPoint(i, vcent[i], ratV26);
            grv26Rat->SetPointError(i, 0, errRatV26);
        
            grv28Rat->SetPoint(i, vcent[i], ratV28);
            grv28Rat->SetPointError(i, 0, errRatV28);
        
        
        }

  }


  TCanvas* cEps0 = new TCanvas("cEps0", "cEps0");
  cEps0->cd();

  TH1D* hdumeps = new TH1D("hdumeps", "; centrality percentile; #epsilon_{0}", 100, 0, 100);
  hdumeps->SetMaximum(1);
  hdumeps->SetMinimum(0);
  hdumeps->Draw();

  grEps0->SetMarkerStyle(20);
  grEps0->SetMarkerColor(1);
  grEps0->SetLineColor(1);
  grEps0->Draw("sameP");


  TCanvas* cAlph = new TCanvas("cAlph", "cAlph");
  cAlph->cd();

  TH1D* hdumalp = new TH1D("hdumalp", "; centrality percentile; #alpha", 100, 0, 100);
  hdumalp->SetMaximum(180);
  hdumalp->SetMinimum(0);
  hdumalp->Draw();

  grAlpha->SetMarkerStyle(20);
  grAlpha->SetMarkerColor(1);
  grAlpha->SetLineColor(1);
  grAlpha->Draw("sameP");
  


  TCanvas* cKap = new TCanvas("cKap", "cKap");
  cKap->cd();

  TH1D* hdumkap = new TH1D("hdumkap", "; centrality percentile; #Kappa", 100, 0, 100);
  hdumkap->SetMaximum(1);
  hdumkap->SetMinimum(0);
  hdumkap->Draw();

  grKappa->SetMarkerStyle(20);
  grKappa->SetMarkerColor(1);
  grKappa->SetLineColor(1);
  grKappa->Draw("sameP");
    
    
    
    TCanvas* cKappr = new TCanvas("cKappr", "cKappr");
    cKappr->cd();
    
    TH1D* hdumkappr = new TH1D("hdumkappr", "; centrality percentile; #Kappa^{'}", 100, 0, 100);
    hdumkappr->SetMaximum(1);
    hdumkappr->SetMinimum(0);
    hdumkappr->Draw();
    
    grKappaPr->SetMarkerStyle(20);
    grKappaPr->SetMarkerColor(1);
    grKappaPr->SetLineColor(1);
    grKappaPr->Draw("sameP");
    
    
    
   if (optDrawEP){
       
       TCanvas* cv2 = new TCanvas("cv2", "cv2");
       cv2->cd();
       
       TH1D* hdumv2 = new TH1D("hdumv2", "; centrality percentile; v_{2}", 100, 0, 100);
       hdumv2->SetMaximum(0.15);
       hdumv2->SetMinimum(0);
       hdumv2->Draw();
       
       grv22->SetLineColor(1);
       grv22->SetMarkerColor(1);
       grv22->SetMarkerStyle(20);
       grv22->Draw("Psame");
       
       grv22EP->SetLineColor(2);
       grv22EP->SetMarkerColor(2);
       grv22EP->SetMarkerStyle(24);
       grv22EP->Draw("Psame");
       
       
       grv24->SetLineColor(1);
       grv24->SetMarkerColor(1);
       grv24->SetMarkerStyle(21);
       grv24->Draw("Psame");
       
       grv24EP->SetLineColor(2);
       grv24EP->SetMarkerColor(2);
       grv24EP->SetMarkerStyle(25);
       grv24EP->Draw("Psame");
       
       
       grv26->SetLineColor(kGreen+2);
       grv26->SetMarkerColor(kGreen+2);
       grv26->SetMarkerStyle(22);
       grv26->Draw("Psame");
       
       grv26EP->SetLineColor(4);
       grv26EP->SetMarkerColor(4);
       grv26EP->SetMarkerStyle(26);
       grv26EP->Draw("Psame");
       
       
       grv28->SetLineColor(kMagenta+2);
       grv28->SetMarkerColor(kMagenta+2);
       grv28->SetMarkerStyle(23);
       grv28->Draw("Psame");
       
       grv28EP->SetLineColor(kCyan+2);
       grv28EP->SetMarkerColor(kCyan+2);
       grv28EP->SetMarkerStyle(32);
       grv28EP->Draw("Psame");
       
       
       TLegend* lg = new TLegend(0.12, 0.53, 0.28, 0.88);
       lg->SetFillColor(0);
       lg->SetBorderSize(0);
       lg->AddEntry(grv22, "v_{2}{2, |#Delta#eta|>1} (Data)", "LP");
       lg->AddEntry(grv24, "v_{2}{4} (Data)", "LP");
       lg->AddEntry(grv26, "v_{2}{6} (Data)", "LP");
       lg->AddEntry(grv28, "v_{2}{8} (Data)", "LP");
       lg->AddEntry(grv22EP, "v_{2}{2, |#Delta#eta|>1} (Fit)", "LP");
       lg->AddEntry(grv24EP, "v_{2}{4} (Fit)", "LP");
       lg->AddEntry(grv26EP, "v_{2}{6} (Fit)", "LP");
       lg->AddEntry(grv28EP, "v_{2}{8} (Fit)", "LP");
       lg->Draw();
       
       //cv2->SaveAs("v2_comp_data_fit.png");
       
       
       
       TCanvas* cv2r = new TCanvas("cv2r", "cv2r");
       cv2r->SetGridy();
       cv2r->cd();
       
       TH1D* hdumv2r = new TH1D("hdumv2r", "; centrality percentile; v_{2} (fit/data)", 100, 0, 100);
       hdumv2r->SetMaximum(1.03);
       hdumv2r->SetMinimum(0.95);
       //hdumv2r->SetMaximum(1.03);
       //hdumv2r->SetMinimum(0.96);
       hdumv2r->Draw();
       
       grv22Rat->SetLineColor(1);
       grv22Rat->SetMarkerColor(1);
       grv22Rat->SetMarkerStyle(20);
       grv22Rat->Draw("Psame");
       
       grv24Rat->SetLineColor(2);
       grv24Rat->SetMarkerColor(2);
       grv24Rat->SetMarkerStyle(25);
       grv24Rat->Draw("Psame");
       
       grv26Rat->SetLineColor(4);
       grv26Rat->SetMarkerColor(4);
       grv26Rat->SetMarkerStyle(26);
       grv26Rat->Draw("Psame");
       
       grv28Rat->SetLineColor(kMagenta+2);
       grv28Rat->SetMarkerColor(kMagenta);
       grv28Rat->SetMarkerStyle(32);
       grv28Rat->Draw("Psame");
       
       
       TLegend* lgr = new TLegend(0.3, 0.68, 0.5, 0.88);
       lgr->SetFillColor(0);
       lgr->SetBorderSize(0);
       lgr->AddEntry(grv22Rat, "v_{2}{2, |#Delta#eta|>1}", "LP");
       lgr->AddEntry(grv24Rat, "v_{2}{4}", "LP");
       lgr->AddEntry(grv26Rat, "v_{2}{6}", "LP");
       lgr->AddEntry(grv28Rat, "v_{2}{8}", "LP");
       lgr->Draw();
       
       //cv2r->SaveAs("v2_comp_data_fit_rat.png");
       
    }
    

    //TFile* out = new TFile("param_fitEP_run2_JacopoStatCubicKpr01.root", "RECREATE");
    TFile* out = new TFile("ParameterOutput.root", "RECREATE");
    grEps0->Write();
    grAlpha->Write();
    grKappa->Write();
    grKappaPr->Write();
    if (optDrawEP){
        grv22Rat->Write();
        grv24Rat->Write();
        grv26Rat->Write();
        grv28Rat->Write();
    }
    out->Close();
    
    

}


void fitter() {
  getParamEP();
};
