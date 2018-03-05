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

  UInt_t b_nJet;
  Float_t b_Jet_eta[10000];
  Float_t b_Jet_phi[10000];
  Float_t b_Jet_pt[10000];

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

//  trees->SetBranchStatus("*", true);
  trees->SetBranchAddress("nGenPart",&b_nGenPart);
  trees->SetBranchAddress("GenPart_status",b_GenPart_status);
  trees->SetBranchAddress("GenPart_pdgId",b_GenPart_pdgId);
  trees->SetBranchAddress("GenPart_eta", b_GenPart_eta);
  trees->SetBranchAddress("GenPart_phi", b_GenPart_phi);

  trees->SetBranchAddress("nJet", &b_nJet);
  trees->SetBranchAddress("Jet_eta", b_Jet_eta);
  trees->SetBranchAddress("Jet_phi", b_Jet_phi);
  trees->SetBranchAddress("Jet_pt", b_Jet_pt);

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

  //Event Loop Start!
  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%1000 == 0 ) std::cout << "event check    iev    ----> " << iev << std::endl;
    trees->GetEntry(iev);
    nS = 0;
    nB = 0;
    int q = -1;
    int qb = -1;
    bool isS = false;
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
    int qj = -1;
    int qbj = -1;
    for(auto j=0; j<b_nJet;++j){
      if(DeltaR(b_Jet_eta[j] - b_GenPart_eta[q], DeltaPhi(b_Jet_phi[j], b_GenPart_phi[q])) < 0.5) qj = j;
      else if(DeltaR(b_Jet_eta[j] - b_GenPart_eta[qb], DeltaPhi(b_Jet_phi[j], b_GenPart_phi[qb])) < 0.5) qbj = j;
    }
    if (qj == -1 || qbj == -1) continue; 
    //nSB
    std::vector<int> match;
    //auto vL2D = sqrt(b_Kshort_x*b_Kshort_x + b_Kshort_y*b_Kshort_y);
    //auto vL3D = sqrt(b_Kshort_x*b_Kshort_x + b_Kshort_y*b_Kshort_y + b_Kshort_z*b_Kshort_z);
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
