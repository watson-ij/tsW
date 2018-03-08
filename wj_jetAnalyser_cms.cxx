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
// select define for getting output about EvnetSelection / Gen or Reco
#define evSel
#define Gen
//#define Reco
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
  Int_t b_GenPart_status[1000];
  Int_t b_GenPart_pdgId[1000];
  Float_t b_GenPart_eta[1000];
  Float_t b_GenPart_phi[1000];
  Float_t b_GenPart_pt[1000];

  UInt_t b_nV0GenPart;
  Int_t b_V0GenPart_pdgId[1000];
  Float_t b_V0GenPart_eta[1000];
  Float_t b_V0GenPart_phi[1000];
  Float_t b_V0GenPart_pt[1000];
  Float_t b_V0GenPart_mass[1000];

  UInt_t b_nJet;
  Int_t b_Jet_partonFlavour[1000];
  Int_t b_Jet_jetId[1000];
  Float_t b_Jet_eta[1000];
  Float_t b_Jet_phi[1000];
  Float_t b_Jet_pt[1000];
  Float_t b_Jet_mass[1000];

  UInt_t b_nGenJet;
  Int_t b_GenJet_partonFlavour[1000];
  Float_t b_GenJet_eta[1000];
  Float_t b_GenJet_phi[1000];
  Float_t b_GenJet_pt[1000];
  Float_t b_GenJet_mass[1000];

  UInt_t b_nGenJetAK8;
  Int_t b_GenJetAK8_partonFlavour[1000];
  Float_t b_GenJetAK8_eta[1000];
  Float_t b_GenJetAK8_phi[1000];
  Float_t b_GenJetAK8_pt[1000];

  UInt_t b_nKshort;
  Float_t b_Kshort_x[1000];
  Float_t b_Kshort_y[1000];
  Float_t b_Kshort_z[1000];
  Float_t b_Kshort_eta[1000];
  Float_t b_Kshort_phi[1000];
  Float_t b_Kshort_pt[1000];
  Float_t b_Kshort_mass[1000];

  UInt_t b_nmeson;
  Int_t b_meson_pdgId[1000];
  Float_t b_meson_x[1000];
  Float_t b_meson_y[1000];
  Float_t b_meson_z[1000];
  Float_t b_meson_eta[1000];
  Float_t b_meson_phi[1000];
  Float_t b_meson_pt[1000];
  Float_t b_meson_mass[1000];
  Float_t b_meson_chi2[1000];
  Float_t b_meson_dca[1000];
  Float_t b_meson_lxy[1000];
  Float_t b_meson_angleXY[1000];

  UInt_t b_ncmeson;
  Int_t b_cmeson_pdgId[1000];
  Float_t b_cmeson_x[1000];
  Float_t b_cmeson_y[1000];
  Float_t b_cmeson_z[1000];
  Float_t b_cmeson_eta[1000];
  Float_t b_cmeson_phi[1000];
  Float_t b_cmeson_pt[1000];
  Float_t b_cmeson_mass[1000];
  Float_t b_cmeson_chi2[1000];
  Float_t b_cmeson_dca[1000];
  Float_t b_cmeson_lxy[1000];
  Float_t b_cmeson_angleXY[1000];

  Int_t b_PV_npvs;
  Float_t b_PV_ndof;
  Float_t b_PV_x;
  Float_t b_PV_y;
  Float_t b_PV_z;

  Float_t b_MET_pt;

  UInt_t b_nMuon;
  Int_t b_Muon_pdgId[1000];
  Int_t b_Muon_charge[1000];
  Float_t b_Muon_pt[1000];
  Float_t b_Muon_eta[1000];
  Float_t b_Muon_phi[1000];
  Float_t b_Muon_mass[1000];
  Float_t b_Muon_pfRelIso04_all[1000];
  Bool_t b_Muon_tightId[1000];

  UInt_t b_nElectron;
  Int_t b_Electron_pdgId[1000];
  Int_t b_Electron_charge[1000];
  Int_t b_Electron_cutBased[1000];
  Float_t b_Electron_pt[1000];
  Float_t b_Electron_eta[1000];
  Float_t b_Electron_phi[1000];
  Float_t b_Electron_mass[1000];
  Float_t b_Electron_deltaEtaSC[1000];

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
  trees->SetBranchAddress("V0GenPart_mass", b_V0GenPart_mass);

  //Reco & Gen jet
  trees->SetBranchAddress("nJet", &b_nJet);
  trees->SetBranchAddress("Jet_partonFlavour", b_Jet_partonFlavour);
  trees->SetBranchAddress("Jet_jetId", b_Jet_jetId);
  trees->SetBranchAddress("Jet_eta", b_Jet_eta);
  trees->SetBranchAddress("Jet_phi", b_Jet_phi);
  trees->SetBranchAddress("Jet_pt", b_Jet_pt);
  trees->SetBranchAddress("Jet_mass", b_Jet_mass);

  trees->SetBranchAddress("nGenJet", &b_nGenJet);
  trees->SetBranchAddress("GenJet_partonFlavour", b_GenJet_partonFlavour);
  trees->SetBranchAddress("GenJet_eta", b_GenJet_eta);
  trees->SetBranchAddress("GenJet_phi", b_GenJet_phi);
  trees->SetBranchAddress("GenJet_pt", b_GenJet_pt);
  trees->SetBranchAddress("GenJet_mass", b_GenJet_mass);

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

  //PV
  trees->SetBranchAddress("PV_npvs", &b_PV_npvs);
  trees->SetBranchAddress("PV_ndof", &b_PV_ndof);
  trees->SetBranchAddress("PV_x", &b_PV_x);
  trees->SetBranchAddress("PV_y", &b_PV_y);
  trees->SetBranchAddress("PV_z", &b_PV_z);

  //MET
  trees->SetBranchAddress("MET_pt", &b_MET_pt);

  //Leptons
  trees->SetBranchAddress("nMuon", &b_nMuon);
  trees->SetBranchAddress("Muon_pdgId", b_Muon_pdgId);
  trees->SetBranchAddress("Muon_charge", b_Muon_charge);
  trees->SetBranchAddress("Muon_pt", b_Muon_pt);
  trees->SetBranchAddress("Muon_eta", b_Muon_eta);
  trees->SetBranchAddress("Muon_phi", b_Muon_phi);
  trees->SetBranchAddress("Muon_mass", b_Muon_mass);
  trees->SetBranchAddress("Muon_pfRelIso04_all", b_Muon_pfRelIso04_all);
  trees->SetBranchAddress("Muon_tightId", b_Muon_tightId);

  trees->SetBranchAddress("nElectron", &b_nElectron);
  trees->SetBranchAddress("Electron_pdgId", b_Electron_pdgId);
  trees->SetBranchAddress("Electron_charge", b_Electron_charge);
  trees->SetBranchAddress("Electron_cutBased", b_Electron_cutBased);
  trees->SetBranchAddress("Electron_pt", b_Electron_pt);
  trees->SetBranchAddress("Electron_eta", b_Electron_eta);
  trees->SetBranchAddress("Electron_phi", b_Electron_phi);
  trees->SetBranchAddress("Electron_mass", b_Electron_mass);
  trees->SetBranchAddress("Electron_deltaEtaSC", b_Electron_deltaEtaSC);

  //make output file
  auto out = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("tsw", "tsw");

  //make Branches
  // global parameters
  int nS = 0; int nSM = 0;
  int nB = 0; int nBM = 0;

  //make Histograms
  std::vector<TH1F*> hList;
  TString cutflow_title = "cutflow" + inf;
  TH1F * cutflow = new TH1F("cutflow", cutflow_title, 10,-2,8); // -2 : all events, -1 : PV cut, 0 : events after lepton selection, 1~5 : events after step, 6 : events after quark selection, 7 : events after quark-jet matching
  hList.push_back(cutflow);

  TH1F * hDLCH = new TH1F("dilepCh", "dilepton channel after event selection", 10, 20, 30);hList.push_back(hDLCH);

  TH1F * hKSn = new TH1F("KSn", "number of KS in gen jet", 100, 0, 10);hList.push_back(hKSn);
  TH1F * hKSnS = new TH1F("KSnS", "number of KS in gen s jet", 100, 0, 10);hList.push_back(hKSnS);
  TH1F * hKSnB = new TH1F("KSnB", "number of KS in gen b jet", 100, 0, 10);hList.push_back(hKSnB);

  TH1F * hLn = new TH1F("Ln", "number of Lamb in gen jet", 100, 0, 10);hList.push_back(hLn);
  TH1F * hLnS = new TH1F("LnS", "number of Lamb in gen s jet", 100, 0, 10);hList.push_back(hLnS);
  TH1F * hLnB = new TH1F("LnB", "number of Lamb in gen b jet", 100, 0, 10);hList.push_back(hLnB);

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

  TH1F * hCx = new TH1F("Cx", "C Meson collection;x;", 100, 0., 1.); hList.push_back(hCx);
  TH1F * hVx = new TH1F("Vx", "V0Meson collection;x;", 100, 0., 1.); hList.push_back(hVx);
  TH1F * hMx = new TH1F("Mx", "Meson collection;x;", 100, 0., 1.); hList.push_back(hMx);

  TH1F * hCj = new TH1F("Cj", "C Meson collection;j;", 100, 0., 1.); hList.push_back(hCj);
  TH1F * hVj = new TH1F("Vj", "V0Meson collection;j;", 100, 0., 1.); hList.push_back(hVj);
  TH1F * hMj = new TH1F("Mj", "Meson collection;j;", 100, 0., 1.); hList.push_back(hMj);
////
  TH1F * hMCs = new TH1F("MCs" , "gen s-quark & gen s-jet matching check", 3, 0,2); hList.push_back(hMCs);
  TH1F * hMCb = new TH1F("MCb" , "gen b-quark & gen b-jet matching check", 3, 0,2); hList.push_back(hMCb);

  TH1F * hGx = new TH1F("Gx", "Gen KS collection;x", 100, 0, 1); hList.push_back(hGx);
  TH1F * hGSx = new TH1F("GSx", "Gen KS collection from S;x", 100, 0, 1); hList.push_back(hGSx);
  TH1F * hGBx = new TH1F("GBx", "Gen KS collection from B;x", 100, 0, 1); hList.push_back(hGBx);

  TH1F * hGj = new TH1F("Gj", "Gen KS collection;j", 100, 0, 1); hList.push_back(hGj);
  TH1F * hGSj = new TH1F("GSj", "Gen KS collection from S;j", 100, 0, 1); hList.push_back(hGSj);
  TH1F * hGBj = new TH1F("GBj", "Gen KS collection from B;j", 100, 0, 1); hList.push_back(hGBj);

  TH1F * hGd = new TH1F("Gd", "Gen KS collection;d", 1000, 0, 20); hList.push_back(hGd);
  TH1F * hGSd = new TH1F("GSd", "Gen KS collection from S;d", 1000, 0, 20); hList.push_back(hGSd);
  TH1F * hGBd = new TH1F("GBd", "Gen KS collection from B;d", 1000, 0, 20); hList.push_back(hGBd);
  TH1F * hGdL = new TH1F("GdL", "Gen KS collection;log(d)", 1000, -10, 5); hList.push_back(hGdL);
  TH1F * hGSdL = new TH1F("GSdL", "Gen KS collection from S;log(d)", 1000, -10, 5); hList.push_back(hGSdL);
  TH1F * hGBdL = new TH1F("GBdL", "Gen KS collection from B;log(d)", 1000, -10, 5); hList.push_back(hGBdL);

  TH1F * hLGx = new TH1F("LGx", "Gen Lamb collection;x", 100, 0, 1); hList.push_back(hGx);
  TH1F * hLGSx = new TH1F("LGSx", "Gen Lamb collection from S;x", 100, 0, 1); hList.push_back(hGSx);
  TH1F * hLGBx = new TH1F("LGBx", "Gen Lamb collection from B;x", 100, 0, 1); hList.push_back(hGBx);

  TH1F * hLGj = new TH1F("LGj", "Gen Lamb collection;j", 100, 0, 1); hList.push_back(hLGj);
  TH1F * hLGSj = new TH1F("LGSj", "Gen Lamb collection from S;j", 100, 0, 1); hList.push_back(hLGSj);
  TH1F * hLGBj = new TH1F("LGBj", "Gen Lamb collection from B;j", 100, 0, 1); hList.push_back(hLGBj);

  TH1F * hLGd = new TH1F("LGd", "Gen Lamb collection;d", 1000, 0, 20); hList.push_back(hLGd);
  TH1F * hLGSd = new TH1F("LGSd", "Gen Lamb collection from S;d", 1000, 0, 20); hList.push_back(hLGSd);
  TH1F * hLGBd = new TH1F("LGBd", "Gen Lamb collection from B;d", 1000, 0, 20); hList.push_back(hLGBd);
  TH1F * hLGdL = new TH1F("LGdL", "Gen Lamb collection;log(d)", 1000, -10, 5); hList.push_back(hLGdL);
  TH1F * hLGSdL = new TH1F("LGSdL", "Gen Lamb collection from S;log(d)", 1000, -10, 5); hList.push_back(hLGSdL);
  TH1F * hLGBdL = new TH1F("LGBdL", "Gen Lamb collection from B;log(d)", 1000, -10, 5); hList.push_back(hLGBdL);



  //Event Loop Start!
  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%1000 == 0 ) std::cout << "event check    iev    ----> " << iev << std::endl;
    trees->GetEntry(iev);

    cutflow->Fill(-2);
#ifdef evSel
    //Object selection
    std::vector<struct Lepton> recolep;
    for(auto i=0; i<b_nMuon; ++i){
      if(!b_Muon_tightId[i]) continue;
      if(b_Muon_pt[i] < 20) continue;
      if(fabs(b_Muon_eta[i]) > 2.4) continue;
      if(b_Muon_pfRelIso04_all[i] > 0.15) continue;

      struct Lepton mu;
      mu.tlv.SetPtEtaPhiM(b_Muon_pt[i], b_Muon_eta[i], b_Muon_phi[i], b_Muon_mass[i]);
      mu.charge = b_Muon_charge[i];
      mu.pdgid = b_Muon_pdgId[i];

      recolep.push_back(mu);
    }
    for(auto i=0; i<b_nElectron; ++i){
      if(b_Electron_pt[i] < 20) continue;
      if(fabs(b_Electron_eta[i]) > 2.4) continue;
      if(b_Electron_cutBased[i] < 3) continue;
      float el_scEta = b_Electron_deltaEtaSC[i] + b_Electron_eta[i];
      if(fabs(el_scEta) > 1.4442 && fabs(el_scEta) < 1.566) continue;

      struct Lepton elec;
      elec.tlv.SetPtEtaPhiM(b_Electron_pt[i], b_Electron_eta[i], b_Electron_phi[i], b_Electron_mass[i]);
      elec.charge = b_Electron_charge[i];
      elec.pdgid = b_Electron_pdgId[i];

      recolep.push_back(elec);
    }

    if(fabs(b_PV_z) >= 24.) continue;
    if(b_PV_npvs == 0 ) continue;
    if(b_PV_ndof < 4) continue;

    cutflow->Fill(-1);

    if (recolep.size() < 2 ) continue;
    //pick highest pt lepton pair
    sort(recolep.begin(), recolep.end(), [](struct Lepton a, struct Lepton b){return a.tlv.Pt() > b.tlv.Pt();});
    recolep.erase(recolep.begin()+2,recolep.end());

    cutflow->Fill(0);

    auto dilepton = recolep[0].tlv + recolep[1].tlv;
    int dilepton_ch = abs(recolep[0].pdgid) + abs(recolep[1].pdgid); // 22 -> ee , 24 -> emu , 26 -> mumu

    // step1
    if(dilepton.M() < 20. || recolep[0].charge * recolep[1].charge > 0 ) continue;
    cutflow->Fill(1);

    // step2
    if ( (dilepton_ch != 24) && ( dilepton.M() > 76. && dilepton.M() < 106.) ) continue;
    cutflow->Fill(2);

    // step3
    if ( (dilepton_ch != 24) && (b_MET_pt < 40.) ) continue;
    cutflow->Fill(3);

    // step4
    std::vector<int> selectedJets;
    for(auto i=0; i<b_nJet; ++i){
      if(b_Jet_pt[i] < 30) continue;
      if(fabs(b_Jet_eta[i]) > 2.4) continue;
      if(b_Jet_jetId[i] < 1) continue;
      bool hasOverLap = false;
      TLorentzVector jet;
      jet.SetPtEtaPhiM(b_Jet_pt[i], b_Jet_eta[i], b_Jet_phi[i], b_Jet_mass[i]);
      for (auto lep : recolep) { if (jet.DeltaR(lep.tlv) < 0.4) hasOverLap = true; }
      if(hasOverLap) continue;
      selectedJets.push_back(i);
    }
    if (selectedJets.size() < 2) continue;
    cutflow->Fill(4);

    hDLCH->Fill(dilepton_ch);
#endif

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

    cutflow->Fill(6);

#ifdef Gen
    //Gen Particle & Gen Jet matching check
    std::vector<int> genJet;

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

    cutflow->Fill(7);

    genJet.push_back(qgj);
    genJet.push_back(qbgj);
    for(auto gj : genJet){
      std::vector<int> KsInJet;
      std::vector<int> LambInJet;
      //Find KS in GenJets
      for(auto i=0; i <b_nV0GenPart;++i){
        if(abs(b_V0GenPart_pdgId[i]) != 310) continue;
        auto dr = DeltaR(b_GenJet_eta[gj] - b_V0GenPart_eta[i], DeltaPhi(b_GenJet_phi[gj], b_V0GenPart_phi[i]));

        TLorentzVector GenKS;
        GenKS.SetPtEtaPhiM(b_V0GenPart_pt[i], b_V0GenPart_eta[i], b_V0GenPart_phi[i], b_V0GenPart_mass[i]);
        std::vector<Double_t> PVtoKS = {GenKS.X() - b_PV_x, GenKS.Y() - b_PV_y, GenKS.Z() - b_PV_z};
        auto magMomentum = sqrt(GenKS.Px()*GenKS.Px() + GenKS.Py()*GenKS.Py() + GenKS.Pz()*GenKS.Pz());
        std::vector<Double_t> unitMomentum = { GenKS.Px()/magMomentum, GenKS.Py()/magMomentum, GenKS.Pz()/magMomentum };
        std::vector<Double_t> Cross = cross3D(PVtoKS, unitMomentum);
        auto d = sqrt(Cross[0]*Cross[0] + Cross[1]*Cross[1] + Cross[2]*Cross[2]);

	if(dr < 0.3) {
	  KsInJet.push_back(i);
          if(b_V0GenPart_pt[i] / b_GenJet_pt[gj] < 0.15) continue;
          hGx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[gj]);
          hGj->Fill(dr);
          hGd->Fill(d);
          hGdL->Fill(log(d));
	  if(abs(b_GenJet_partonFlavour[gj]) == 3) {
            hGSx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[gj]);
            hGSj->Fill(dr);
            hGSd->Fill(d);
            hGSdL->Fill(log(d));
	  }
          if(abs(b_GenJet_partonFlavour[gj]) == 5) {
            hGBx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[gj]);
            hGBj->Fill(dr);
            hGBd->Fill(d);
            hGBdL->Fill(log(d));
          }
	}
      }
      hKSn->Fill(KsInJet.size());
      if(abs(b_GenJet_partonFlavour[gj]) == 3) {
	hKSnS->Fill(KsInJet.size());
      }
      if(abs(b_GenJet_partonFlavour[gj]) == 5) {
	hKSnB->Fill(KsInJet.size());
      }
      // Find Lambda in GenJets
      for(auto i=0; i <b_nV0GenPart;++i){
        if(abs(b_V0GenPart_pdgId[i]) != 3122) continue;
        auto dr = DeltaR(b_GenJet_eta[gj] - b_V0GenPart_eta[i], DeltaPhi(b_GenJet_phi[gj], b_V0GenPart_phi[i]));

        TLorentzVector GenLamb;
        GenLamb.SetPtEtaPhiM(b_V0GenPart_pt[i], b_V0GenPart_eta[i], b_V0GenPart_phi[i], b_V0GenPart_mass[i]);
        std::vector<Double_t> PVtoLamb = {GenLamb.X() - b_PV_x, GenLamb.Y() - b_PV_y, GenLamb.Z() - b_PV_z};
        auto magMomentum = sqrt(GenLamb.Px()*GenLamb.Px() + GenLamb.Py()*GenLamb.Py() + GenLamb.Pz()*GenLamb.Pz());
        std::vector<Double_t> unitMomentum = { GenLamb.Px()/magMomentum, GenLamb.Py()/magMomentum, GenLamb.Pz()/magMomentum };
        std::vector<Double_t> Cross = cross3D(PVtoLamb, unitMomentum);
        auto d = sqrt(Cross[0]*Cross[0] + Cross[1]*Cross[1] + Cross[2]*Cross[2]);

        if(dr < 0.3) {
          LambInJet.push_back(i);
          if(b_V0GenPart_pt[i] / b_GenJet_pt[gj] < 0.15) continue;
          hLGx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[gj]);
          hLGj->Fill(dr);
          hLGd->Fill(d);
          hLGdL->Fill(log(d));
          if(abs(b_GenJet_partonFlavour[gj]) == 3) {
            hLGSx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[gj]);
            hLGSj->Fill(dr);
            hLGSd->Fill(d);
            hLGSdL->Fill(log(d));
          }
          if(abs(b_GenJet_partonFlavour[gj]) == 5) {
            hLGBx->Fill(b_V0GenPart_pt[i]/b_GenJet_pt[gj]);
            hLGBj->Fill(dr);
            hLGBd->Fill(d);
            hLGBdL->Fill(log(d));
          }
        }
      }
      hLn->Fill(LambInJet.size());
      if(abs(b_GenJet_partonFlavour[gj]) == 3) {
        hLnS->Fill(LambInJet.size());
      }
      if(abs(b_GenJet_partonFlavour[gj]) == 5) {
        hLnB->Fill(LambInJet.size());
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

    cutflow->Fill(7);

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
