#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "classes/DelphesClasses.h"
#include "wj_jetAnalyser.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// select define for getting output about Gen or Reco
//#define Gen
#define Reco
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

  auto inf = std::string{argv[1]};
  auto outf = std::string{argv[2]};

  // check cpu time (start)
  std::clock_t c_start = std::clock();

  //read input file
  auto tfiles = TFile::Open(inf.c_str(), "READ");
  auto trees = (TTree*) tfiles->Get("Events");

  UInt_t b_nGenPart;
  Int_t b_GenPart_status[10000];
  Int_t b_GenPart_pdgId[10000];
  Float_t b_GenPart_eta[10000];
  Float_t b_GenPart_phi[10000];
  Float_t b_GenPart_pt[10000];

  UInt_t b_nV0GenPart;
  Int_t b_V0GenPart_pdgId[10000];
  Float_t b_V0GenPart_eta[10000];
  Float_t b_V0GenPart_phi[10000];
  Float_t b_V0GenPart_pt[10000];

  UInt_t b_nJet;
  Float_t b_Jet_eta[10000];
  Float_t b_Jet_phi[10000];
  Float_t b_Jet_pt[10000];

  UInt_t b_nGenJet;
  Int_t b_GenJet_partonFlavour[10000];
  Float_t b_GenJet_eta[10000];
  Float_t b_GenJet_phi[10000];
  Float_t b_GenJet_pt[10000];

  UInt_t b_nGenJetAK8;
  Int_t b_GenJetAK8_partonFlavour[10000];
  Float_t b_GenJetAK8_eta[10000];
  Float_t b_GenJetAK8_phi[10000];
  Float_t b_GenJetAK8_pt[10000];

  UInt_t b_nKshort;
  Float_t b_Kshort_x[10000];
  Float_t b_Kshort_y[10000];
  Float_t b_Kshort_z[10000];
  Float_t b_Kshort_eta[10000];
  Float_t b_Kshort_phi[10000];
  Float_t b_Kshort_pt[10000];
  Float_t b_Kshort_mass[10000];

  UInt_t b_nmeson;
  Int_t b_meson_pdgId[10000];
  Float_t b_meson_x[10000];
  Float_t b_meson_y[10000];
  Float_t b_meson_z[10000];
  Float_t b_meson_eta[10000];
  Float_t b_meson_phi[10000];
  Float_t b_meson_pt[10000];
  Float_t b_meson_mass[10000];
  Float_t b_meson_chi2[10000];
  Float_t b_meson_dca[10000];
  Float_t b_meson_lxy[10000];
  Float_t b_meson_angleXY[10000];

  UInt_t b_ncmeson;
  Int_t b_cmeson_pdgId[10000];
  Float_t b_cmeson_x[10000];
  Float_t b_cmeson_y[10000];
  Float_t b_cmeson_z[10000];
  Float_t b_cmeson_eta[10000];
  Float_t b_cmeson_phi[10000];
  Float_t b_cmeson_pt[10000];
  Float_t b_cmeson_mass[10000];
  Float_t b_cmeson_chi2[10000];
  Float_t b_cmeson_dca[10000];
  Float_t b_cmeson_lxy[10000];
  Float_t b_cmeson_angleXY[10000];

  //GenPart
  trees->SetBranchAddress("nGenPart",&b_nGenPart);
  trees->SetBranchAddress("GenPart_status",b_GenPart_status);
  trees->SetBranchAddress("GenPart_pdgId",b_GenPart_pdgId);
  trees->SetBranchAddress("GenPart_eta", b_GenPart_eta);
  trees->SetBranchAddress("GenPart_phi", b_GenPart_phi);
  trees->SetBranchAddress("GenPart_pt", b_GenPart_pt);

  trees->SetBranchAddress("nV0GenPart",&b_nV0GenPart);
  trees->SetBranchAddress("V0GenPart_pdgId",b_V0GenPart_pdgId);
  trees->SetBranchAddress("V0GenPart_eta", b_V0GenPart_eta);
  trees->SetBranchAddress("V0GenPart_phi", b_V0GenPart_phi);
  trees->SetBranchAddress("V0GenPart_pt", b_V0GenPart_pt);


  //Reco & Gen jet
  trees->SetBranchAddress("nJet", &b_nJet);
  trees->SetBranchAddress("Jet_eta", b_Jet_eta);
  trees->SetBranchAddress("Jet_phi", b_Jet_phi);
  trees->SetBranchAddress("Jet_pt", b_Jet_pt);

  trees->SetBranchAddress("nGenJet", &b_nGenJet);
  trees->SetBranchAddress("GenJet_partonFlavour", b_GenJet_partonFlavour);
  trees->SetBranchAddress("GenJet_eta", b_GenJet_eta);
  trees->SetBranchAddress("GenJet_phi", b_GenJet_phi);
  trees->SetBranchAddress("GenJet_pt", b_GenJet_pt);

  trees->SetBranchAddress("nGenJetAK8", &b_nGenJetAK8);
  trees->SetBranchAddress("GenJetAK8_partonFlavour", b_GenJetAK8_partonFlavour);
  trees->SetBranchAddress("GenJetAK8_eta", b_GenJetAK8_eta);
  trees->SetBranchAddress("GenJetAK8_phi", b_GenJetAK8_phi);
  trees->SetBranchAddress("GenJetAK8_pt", b_GenJetAK8_pt);

  //Meson from several producer
  trees->SetBranchAddress("nKshort", &b_nKshort);
  trees->SetBranchAddress("Kshort_eta", b_Kshort_eta);
  trees->SetBranchAddress("Kshort_phi", b_Kshort_phi);
  trees->SetBranchAddress("Kshort_x", b_Kshort_x);
  trees->SetBranchAddress("Kshort_y", b_Kshort_y);
  trees->SetBranchAddress("Kshort_z", b_Kshort_z);
  trees->SetBranchAddress("Kshort_pt", b_Kshort_pt);
  trees->SetBranchAddress("Kshort_mass", b_Kshort_mass);

  trees->SetBranchAddress("nmeson", &b_nmeson);
  trees->SetBranchAddress("meson_pdgId", b_meson_pdgId);
  trees->SetBranchAddress("meson_eta", b_meson_eta);
  trees->SetBranchAddress("meson_phi", b_meson_phi);
  trees->SetBranchAddress("meson_x", b_meson_x);
  trees->SetBranchAddress("meson_y", b_meson_y);
  trees->SetBranchAddress("meson_z", b_meson_z);
  trees->SetBranchAddress("meson_pt", b_meson_pt);
  trees->SetBranchAddress("meson_mass", b_meson_mass);
  trees->SetBranchAddress("meson_chi2", b_meson_chi2);
  trees->SetBranchAddress("meson_dca", b_meson_dca);
  trees->SetBranchAddress("meson_lxy", b_meson_lxy);
  trees->SetBranchAddress("meson_angleXY", b_meson_angleXY);

  trees->SetBranchAddress("ncmeson", &b_ncmeson);
  trees->SetBranchAddress("cmeson_pdgId", b_cmeson_pdgId);
  trees->SetBranchAddress("cmeson_eta", b_cmeson_eta);
  trees->SetBranchAddress("cmeson_phi", b_cmeson_phi);
  trees->SetBranchAddress("cmeson_x", b_cmeson_x);
  trees->SetBranchAddress("cmeson_y", b_cmeson_y);
  trees->SetBranchAddress("cmeson_z", b_cmeson_z);
  trees->SetBranchAddress("cmeson_pt", b_cmeson_pt);
  trees->SetBranchAddress("cmeson_mass", b_cmeson_mass);
  trees->SetBranchAddress("cmeson_chi2", b_cmeson_chi2);
  trees->SetBranchAddress("cmeson_dca", b_cmeson_dca);
  trees->SetBranchAddress("cmeson_lxy", b_cmeson_lxy);
  trees->SetBranchAddress("cmeson_angleXY", b_cmeson_angleXY);

  //make output file
  auto out = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("tsw", "tsw");

  //make Branches

  //make Histograms
  TString cutflow_title = "cutflow" + inf;
  TH1F * cutflow = new TH1F("cutflow", cutflow_title, 7,-1,6); // -1 : all events, 0 : events after lepton selection, 1~ : events after step

  // global parameters
  int nS = 0; int nSM = 0;
  int nB = 0; int nBM = 0;

  //nSB 
  //nVC

  std::vector<TH1F*> hList;
  TH1F * hC = new TH1F("C", "C Meson collection;Mass [GeV];", 50, 0.43, 0.57); hList.push_back(hC);
  TH1F * hCc = new TH1F("Cc", "C Meson collection;Mass [GeV];", 50, 0.43, 0.57); hList.push_back(hCc);
  TH1F * hV = new TH1F("V", "V0Meson collection;Mass [GeV];", 50, 0.43, 0.57); hList.push_back(hV);
  TH1F * hM = new TH1F("M", "Meson collection;Mass [GeV];", 50, 0.43, 0.57);  hList.push_back(hM);
  TH1F * hMc = new TH1F("Mc", "Meson collection;Mass [GeV];", 50, 0.43, 0.57);  hList.push_back(hMc);

  TH1F * hC2D = new TH1F("C2D", "C Meson collection;Decay Length 2D [cm];", 50, 0., 10); hList.push_back(hC2D);
  TH1F * hV2D = new TH1F("V2D", "V0Meson collection;Decay Length 2D [cm];", 50, 0., 10); hList.push_back(hV2D);
  TH1F * hM2D = new TH1F("M2D", "Meson collection;Decay Length 2D [cm];", 50, 0., 10); hList.push_back(hM2D);

  TH1F * hC3D = new TH1F("C3D", "C Meson collection;Decay Length 3D [cm];", 50, 0., 25); hList.push_back(hC3D);
  TH1F * hV3D = new TH1F("V3D", "V0Meson collection;Decay Length 3D [cm];", 50, 0., 25); hList.push_back(hV3D);
  TH1F * hM3D = new TH1F("M3D", "Meson collection;Decay Length 3D [cm];", 50, 0., 25); hList.push_back(hM3D);

  TH1F * hCx = new TH1F("Cx", "C Meson collection;x;", 50, 0., 1.); hList.push_back(hCx);
  TH1F * hVx = new TH1F("Vx", "V0Meson collection;x;", 50, 0., 1.); hList.push_back(hVx);
  TH1F * hMx = new TH1F("Mx", "Meson collection;x;", 50, 0., 1.); hList.push_back(hMx);

  TH1F * hCj = new TH1F("Cj", "C Meson collection;j;", 50, 0., 1.); hList.push_back(hCj);
  TH1F * hVj = new TH1F("Vj", "V0Meson collection;j;", 50, 0., 1.); hList.push_back(hVj);
  TH1F * hMj = new TH1F("Mj", "Meson collection;j;", 50, 0., 1.); hList.push_back(hMj);
////
  TH1F * hMCs = new TH1F("MCs" , "gen s-quark & gen s-jet matching check", 3, 0,2); hList.push_back(hMCs);
  TH1F * hMCb = new TH1F("MCb" , "gen b-quark & gen b-jet matching check", 3, 0,2); hList.push_back(hMCb);
  //TH1F * hMCsAK8 = new TH1F("MCsAK8" , "gen s-quark & gen s-jet(AK8) matching check", 3, 0,2); hList.push_back(hMCsAK8);
  //TH1F * hMCbAK8 = new TH1F("MCbAK8" , "gen b-quark & gen b-jet(AK8) matching check", 3, 0,2); hList.push_back(hMCbAK8);

  TH1F * hGx = new TH1F("Gx", "x_KS as Gen KS vs Gen Jet", 100, 0, 1); hList.push_back(hGx);
  TH1F * hGSx = new TH1F("GSx", "x_KS as Gen KS vs Gen S Jet", 100, 0, 1); hList.push_back(hGSx);
  TH1F * hGBx = new TH1F("GBx", "x_KS as Gen KS vs Gen B Jet", 100, 0, 1); hList.push_back(hGBx);

  TH1F * hRx = new TH1F("Rx", "x_KS as Reco KS vs Reco Jet", 100, 0, 1); hList.push_back(hRx);

  //Event Loop Start!
  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%1000 == 0 ) std::cout << "event check    iev    ----> " << iev << std::endl;
    trees->GetEntry(iev);
    nS = 0;
    nB = 0;
    int q = -1;
    int qb = -1;
    bool isS = false;
    //Find s/b quark from Gen Info.
    for(auto i=0; i<b_nGenPart; ++i){
      if (std::abs(b_GenPart_status[i] - 25) < 5 && b_GenPart_pdgId[i] == 3) {
	q = i;
	nS += 1; 
	isS = true;
      }
      if (std::abs(b_GenPart_status[i] - 25) < 5 && b_GenPart_pdgId[i] == -3) {
        qb = i; 
        nS += 1; 
      }
      if (std::abs(b_GenPart_status[i] - 25) < 5 && b_GenPart_pdgId[i] == 5) {
        q = i; 
        nB += 1; 
      }
      if (std::abs(b_GenPart_status[i] - 25) < 5 && b_GenPart_pdgId[i] == -5) {
        qb = i; 
        nB += 1; 
      }
    }
    if (q == -1 ) continue;
#ifdef Gen
    //Gen Particle & Gen Jet matching check
    int qgjS = -1;
    int qbgjS = -1;
    int qgjB = -1;
    int qbgjB = -1;
    int qgj = -1;
    int qbgj = -1;
    for(auto j=0; j<b_nGenJet; ++j){
      auto dr = DeltaR(b_GenJet_eta[j] - b_GenPart_eta[q], DeltaPhi(b_GenJet_phi[j], b_GenPart_phi[q]));
      auto drb = DeltaR(b_GenJet_eta[j] - b_GenPart_eta[qb], DeltaPhi(b_GenJet_phi[j], b_GenPart_phi[qb]));
      bool isMat = false;
      if(b_GenJet_partonFlavour[j] == 3) {
        if(dr < 0.5) {
	  isMat = true;
	  hMCs->Fill(isMat);
	  qgjS = j;
	  qgj = j;
	}
      }
      if(b_GenJet_partonFlavour[j] == -3) {
        if(drb < 0.5) {
          isMat = true;
          hMCs->Fill(isMat);
          qbgjS = j;
          qbgj = j;
        }
      }
      if(b_GenJet_partonFlavour[j] == 5) {
        if(dr < 0.5) {
          isMat = true;
          hMCb->Fill(isMat);
          qgjB = j;
          qgj = j;
        }
      }
      if(b_GenJet_partonFlavour[j] == -5) {
        if(drb < 0.5) {
          isMat = true;
          hMCs->Fill(isMat);
          qbgjB = j;
          qbgj = j;
        }
      }
    }
    if(qgj == -1 || qbgj == -1) continue;
    int sKS = -1;
    int bKS = -1;
    std::vector<int> genKsInJet;
    std::vector<int> genKsInS;
    std::vector<int> genKsInB;
    //Find KS in jets
    for(auto i=0; i<b_nV0GenPart; ++i){
      if(abs(b_V0GenPart_pdgId[i]) != 310) continue;
      auto drSKS = DeltaR(b_GenJet_eta[qgjS] - b_V0GenPart_eta[i], DeltaPhi(b_GenJet_phi[qgjS], b_V0GenPart_phi[i]));
      auto drbSKS = DeltaR(b_GenJet_eta[qbgjS] - b_V0GenPart_eta[i], DeltaPhi(b_GenJet_phi[qbgjS], b_V0GenPart_phi[i]));
      auto drBKS = DeltaR(b_GenJet_eta[qgjB] - b_V0GenPart_eta[i], DeltaPhi(b_GenJet_phi[qgjB], b_V0GenPart_phi[i]));
      auto drbBKS = DeltaR(b_GenJet_eta[qbgjB] - b_V0GenPart_eta[i], DeltaPhi(b_GenJet_phi[qbgjB], b_V0GenPart_phi[i]));

      if(drSKS < 0.3 && qgjS != -1) {
        if(b_V0GenPart_pt[i]/b_GenJet_pt[qgjS] >= 0.15){
          genKsInJet.push_back(i);
          sKS = i;
          genKsInS.push_back(sKS);
          hGx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[qgjS]);
	  hGSx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[qgjS]);
        }
      }
      if(drbSKS < 0.3 && qbgjS != -1) {
        if(b_V0GenPart_pt[i]/b_GenJet_pt[qbgjS] >= 0.15){
          genKsInJet.push_back(i);
          sKS = i;
          genKsInS.push_back(sKS);
          hGx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[qbgjS]);
          hGSx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[qbgjS]);
        }
      }
      if(drBKS < 0.3 && qgjB != -1) {
        if(b_V0GenPart_pt[i]/b_GenJet_pt[qgjB] >= 0.15){
          genKsInJet.push_back(i);
          bKS = i;
          genKsInB.push_back(bKS);
          hGx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[qgjB]);
          hGBx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[qgjB]);
        }
      }
      if(drbBKS < 0.3 && qbgjB != -1) {
        if(b_V0GenPart_pt[i]/b_GenJet_pt[qbgjB] >= 0.15){
          genKsInJet.push_back(i);
          bKS = i;
          genKsInB.push_back(bKS);
          hGx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[qbgjB]);
          hGBx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[qbgjB]);
        }
      }
    }
#endif

#ifdef Reco
    int qj = -1;
    int qbj = -1;
    for(auto j=0; j<b_nJet;++j){
      if(DeltaR(b_Jet_eta[j] - b_GenPart_eta[q], DeltaPhi(b_Jet_phi[j], b_GenPart_phi[q])) < 0.5) qj = j;
      else if(DeltaR(b_Jet_eta[j] - b_GenPart_eta[qb], DeltaPhi(b_Jet_phi[j], b_GenPart_phi[qb])) < 0.5) qbj = j;
    }
    if (qj == -1 || qbj == -1) continue; 

    std::vector<int> match;
    //KS from V0Producer
    for(auto k=0; k<b_nKshort; ++k){
      auto dr = DeltaR(b_Kshort_eta[k] - b_Jet_eta[qj], DeltaPhi(b_Kshort_phi[k],b_Jet_phi[qj]));
      auto drb = DeltaR(b_Kshort_eta[k] - b_Jet_eta[qbj], DeltaPhi(b_Kshort_phi[k],b_Jet_phi[qbj]));
      if(dr < 0.3){
	if((b_Kshort_pt[k]/b_Jet_pt[qj]) < 0.15) continue;
        if(isS){
	  nSM += 1;
	  match.push_back(k);
	  hV->Fill(b_Kshort_mass[k]);
	  if (fabs(b_Kshort_mass[k] - 0.50) < 0.02){
	    hV2D->Fill(vL2D(b_Kshort_x[k], b_Kshort_y[k]));
	    hVx->Fill(b_Kshort_pt[k]/b_Jet_pt[qj]);
	    hVj->Fill(dr);
	  }
	}
	else{
	  nBM += 1;
	  hV->Fill(b_Kshort_mass[k]);
          if (fabs(b_Kshort_mass[k] - 0.50) < 0.02){
            hV2D->Fill(vL2D(b_Kshort_x[k], b_Kshort_y[k]));
            hVx->Fill(b_Kshort_pt[k]/b_Jet_pt[qj]);
            hVj->Fill(dr);
          }
	}
      }
      else if(drb < 0.3){
	if((b_Kshort_pt[k]/b_Jet_pt[qbj]) < 0.15) continue;
	if(isS){
	  nSM +=1;
	  match.push_back(k);
          hV->Fill(b_Kshort_mass[k]);
          if (fabs(b_Kshort_mass[k] - 0.50) < 0.02){
            hV2D->Fill(vL2D(b_Kshort_x[k], b_Kshort_y[k]));
	    hV3D->Fill(vL3D(b_Kshort_x[k], b_Kshort_y[k], b_Kshort_z[k]));
            hVx->Fill((b_Kshort_pt[k]/b_Jet_pt[qbj]));
            hVj->Fill(drb);
          }
	}
	else{
          nBM += 1;
          hV->Fill(b_Kshort_mass[k]);
          if (fabs(b_Kshort_mass[k] - 0.50) < 0.02){
            hV2D->Fill(vL2D(b_Kshort_x[k], b_Kshort_y[k]));
	    hV3D->Fill(vL3D(b_Kshort_x[k], b_Kshort_y[k], b_Kshort_z[k]));
            hVx->Fill((b_Kshort_pt[k]/b_Jet_pt[qbj]));
            hVj->Fill(drb);
          }
	}
      }
    }

    // KS from MesonProducer
    for (auto m=0; m<b_nmeson; ++m){
      if(b_meson_pdgId[m] != 310) continue;
      auto dr = DeltaR(b_meson_eta[m] - b_Jet_eta[qj], DeltaPhi(b_meson_phi[m],b_Jet_phi[qj]));
      auto drb = DeltaR(b_meson_eta[m] - b_Jet_eta[qbj], DeltaPhi(b_meson_phi[m],b_Jet_phi[qbj]));
      if(dr<0.3){
	if((b_meson_pt[m]/b_Jet_pt[qj]) < 0.15) continue;
	if(isS){
	  nSM += 1;
	  match.push_back(m);
	  hM->Fill(b_meson_mass[m]);
	  if(fabs(b_meson_chi2[m]) < 1.0 && b_meson_dca[m] < 0.2 && b_meson_lxy[m] > 0.25 && b_meson_angleXY[m] > 0){
	    hMc->Fill(b_meson_mass[m]);
	    if(fabs(b_meson_mass[m] - 0.50) < 0.02){
	      hM2D->Fill(vL2D(b_meson_x[m], b_meson_y[m]));
	      hMx->Fill(b_meson_pt[m]/b_Jet_pt[qj]);
	      hMj->Fill(dr);
	    }
	  } 
	}
	else{
	  nBM += 1;
	  hM->Fill(b_meson_mass[m]);
	  if(fabs(b_meson_chi2[m]) < 1.0 && b_meson_dca[m] < 0.2 && b_meson_lxy[m] > 0.25 && b_meson_angleXY[m] > 0){
            hMc->Fill(b_meson_mass[m]);
            if(fabs(b_meson_mass[m] - 0.50) < 0.02){
              hM2D->Fill(vL2D(b_meson_x[m], b_meson_y[m]));
              hMx->Fill(b_meson_pt[m]/b_Jet_pt[qj]);
              hMj->Fill(dr);
            }	 
	  }
	}
      }
      else if(drb < 0.3){
        if((b_meson_pt[m]/b_Jet_pt[qbj]) < 0.15) continue;
        if(isS){
          nSM += 1;
          match.push_back(m);
          hM->Fill(b_meson_mass[m]);
          if(fabs(b_meson_chi2[m]) < 1.0 && b_meson_dca[m] < 0.2 && b_meson_lxy[m] > 0.25 && b_meson_angleXY[m] > 0){
            hMc->Fill(b_meson_mass[m]);
            if(fabs(b_meson_mass[m] - 0.50) < 0.02){
              hM2D->Fill(vL2D(b_meson_x[m], b_meson_y[m]));
	      hM3D->Fill(vL3D(b_meson_x[m], b_meson_y[m], b_meson_z[m]));
              hMx->Fill(b_meson_pt[m]/b_Jet_pt[qbj]);
              hMj->Fill(drb);
            }
          }  
        }
        else{
          nBM += 1;
          hM->Fill(b_meson_mass[m]);
          if(fabs(b_meson_chi2[m]) < 1.0 && b_meson_dca[m] < 0.2 && b_meson_lxy[m] > 0.25 && b_meson_angleXY[m] > 0){
            hMc->Fill(b_meson_mass[m]);
            if(fabs(b_meson_mass[m] - 0.50) < 0.02){
              hM2D->Fill(vL2D(b_meson_x[m], b_meson_y[m]));
              hM3D->Fill(vL3D(b_meson_x[m], b_meson_y[m], b_meson_z[m]));
              hMx->Fill(b_meson_pt[m]/b_Jet_pt[qbj]);
              hMj->Fill(drb);
            }
          }
        }
      }
    }

    // KS from CMesonProducer
    for (auto c=0; c<b_ncmeson; ++c){
      if(b_cmeson_angleXY[c] < 0) continue;
      if(b_cmeson_pdgId[c] == 310){
        auto dr = DeltaR(b_cmeson_eta[c] - b_Jet_eta[qj], DeltaPhi(b_cmeson_phi[c],b_Jet_phi[qj]));
        auto drb = DeltaR(b_cmeson_eta[c] - b_Jet_eta[qbj], DeltaPhi(b_cmeson_phi[c],b_Jet_phi[qbj]));
        if(dr < 0.3){
          if((b_cmeson_pt[c]/b_Jet_pt[qj]) < 0.15) continue;
          hC->Fill(b_cmeson_mass[c]);
          if(fabs(b_cmeson_chi2[c]) < 1.0 && b_cmeson_dca[c] < 0.2 && b_cmeson_lxy[c] > 0.25 && b_cmeson_angleXY[c] > 0){
            hCc->Fill(b_cmeson_mass[c]);
            if(fabs(b_cmeson_mass[c] - 0.50) < 0.02){
              hC2D->Fill(vL2D(b_cmeson_x[c], b_cmeson_y[c]));
              hC3D->Fill(vL3D(b_cmeson_x[c], b_cmeson_y[c], b_cmeson_z[c]));
              hCx->Fill(b_cmeson_pt[c]/b_Jet_pt[qj]);
              hCj->Fill(dr);
            }
          }
        }
        else if(drb < 0.3){
          if((b_cmeson_pt[c]/b_Jet_pt[qbj]) < 0.15) continue;
          hC->Fill(b_cmeson_mass[c]);
          if(fabs(b_cmeson_chi2[c]) < 1.0 && b_cmeson_dca[c] < 0.2 && b_cmeson_lxy[c] > 0.25 && b_cmeson_angleXY[c] > 0){
            hCc->Fill(b_cmeson_mass[c]);
            if(fabs(b_cmeson_mass[c] - 0.50) < 0.02){
              hC2D->Fill(vL2D(b_cmeson_x[c], b_cmeson_y[c]));
              hC3D->Fill(vL3D(b_cmeson_x[c], b_cmeson_y[c], b_cmeson_z[c]));
              hCx->Fill(b_cmeson_pt[c]/b_Jet_pt[qbj]);
              hCj->Fill(drb);
            }
          }   
        }
      }
    } 
#endif  
  }  

  tfiles->Close();
  outtr->Write();

  for(auto h=0; h<hList.size(); ++h){
    hList[h]->Write();
  }

  out->Close();

  //check cpu time (end)
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  std::cout << "CPU time used(sec): " << time_elapsed_ms/1000 << " sec\n";
  return 0;
}
