#include <stdio.h>
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
#include "wj_analysis.h"

int main(int argc, char* argv[])
{
  auto inf = std::string{argv[1]};
  auto outf = std::string{argv[2]};

  // check cpu time (start)
  std::clock_t c_start = std::clock();

  //read input file
  auto tfiles = TFile::Open(inf.c_str(), "READ");
  auto trees = (TTree*) tfiles->Get("Delphes");
  trees->SetBranchStatus("*", true);
  trees->SetBranchAddress("GenJet",      &gen_jets);
  trees->SetBranchAddress("Particle",    &particles);
  trees->SetBranchAddress("Electron",    &electrons);
  trees->SetBranchAddress("Muon",        &muons);
  trees->SetBranchAddress("Jet",         &jets);
  trees->SetBranchAddress("MissingET",   &missingET);
  trees->SetBranchAddress("Track",       &tracks);

  //make out file
  auto out = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("tsw", "tsw");

  BranchI(nJets);
  BranchO(matched);
  BranchI(nMatchedJets);
  BranchF(dr); BranchF(drTrue);

  BranchVF(massTrackPair);
  BranchVF(massPion1); BranchVF(massPionTrue1); BranchVF(massPion2); BranchVF(massPionTrue2);

  TH2F * histo_M_dr = new TH2F("hist_M_dr", "hist_M_dr" ,1000, 0, 5, 1000, 0,1);

  //Event Loop Start!
  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%1000 == 0 ) std::cout << "event check    iev    ----> " << iev << std::endl;
    trees->GetEntry(iev);
    float dRCut = 0.5;
    float mass_KS = 0.49761;
    //find s or b
    std::vector<const GenParticle*> genJets;
    for (unsigned i = 0; i < particles->GetEntries(); ++ i){
      auto p = (const GenParticle*) particles->At(i);
      if (p->Status > 30) continue;
      if (abs(p->PID) != 6) continue;
      auto top = getLast(particles,p);
      auto jet = (const GenParticle*) particles->At(top->D2);
      genJets.push_back(jet);
    }
    //find jet match to s or b    
    std::vector<Jet*> jetCand;
    nJets = jets->GetEntries();
    for(size_t i = 0; i < nJets; ++i){
      matched=false;
      nMatchedJets = -1;
      dr = -1;
      drTrue = -1;
      massTrackPair.clear();
      massPion1.clear(); massPionTrue1.clear(); massPion2.clear(); massPionTrue2.clear();

      //jet matching
      auto jet = (Jet*) jets->At(i);
      for (auto& p : genJets) {
        float dR = DeltaR(jet->Eta - p->Eta, DeltaPhi(jet->Phi, p->Phi));
        if (dR < dRCut) {
          jetCand.push_back(jet);
          matched = true;
        }
      }

      if(i==nJets-1) nMatchedJets = jetCand.size();

      if(!matched) continue;
      auto nDau = jet->Constituents.GetEntries();
      //select first track
      for(size_t j=0; j < nDau; ++j){
        auto dau_1 = jet->Constituents.At(j);
        auto dau_trk_1 = dynamic_cast<Track*>(dau_1);
        if(!dau_trk_1) continue;
        else if(abs(dau_trk_1->PID) == 11 || abs(dau_trk_1->PID) == 13) continue;
        TLorentzVector dau_trk_tlv_1;
        dau_trk_tlv_1.SetPtEtaPhiM(dau_trk_1->PT, dau_trk_1->Eta, dau_trk_1->Phi, 0.13957);
        //select second track
        for(auto k= j+1; k < nDau; ++k){
          auto dau_2 = jet->Constituents.At(k);
          auto dau_trk_2 = dynamic_cast<Track*>(dau_2);
          if(!dau_trk_2) continue;
          else if(abs(dau_trk_2->PID) == 11 || abs(dau_trk_2->PID) == 13) continue;
          TLorentzVector dau_trk_tlv_2;
          dau_trk_tlv_2.SetPtEtaPhiM(dau_trk_2->PT, dau_trk_2->Eta, dau_trk_2->Phi, 0.13957);
          auto pair = dau_trk_tlv_1 + dau_trk_tlv_2;
          float mass_pair = pair.M();
          //float mass_pair = dau_trk_tlv_1.M()*cosh(dau_trk_tlv_1.Eta()) + dau_trk_tlv_2.M()*cosh(dau_trk_tlv_2.Eta()); // v = tanh(eta), c = 1 => gamma(v) = cosh(eta)
          massTrackPair.push_back(mass_pair);
          //pair matching
          if ( fabs(mass_pair - mass_KS) < 0.1 ) {
            massPion1.push_back(mass_pair);
            if(abs(dau_trk_1->PID) == 211 & abs(dau_trk_2->PID) == 211) massPionTrue1.push_back(mass_pair);
          }
          if (dau_trk_tlv_1.DeltaR(dau_trk_tlv_2) < fabs(dr)) {
            dr = dau_trk_tlv_1.DeltaR(dau_trk_tlv_2);
            massPion2.push_back(mass_pair);
          }
          if ( abs(dau_trk_1->PID) == 211 & abs(dau_trk_2->PID) == 211){
            if (dau_trk_tlv_1.DeltaR(dau_trk_tlv_2) < fabs(drTrue)) {
              drTrue = dau_trk_tlv_1.DeltaR(dau_trk_tlv_2);
              massPionTrue2.push_back(mass_pair);
            }
          }
        }
      }
      if(massPion2.size() >0) histo_M_dr->Fill(massPion2[massPion2.size()-1],dr);
      outtr->Fill();
    }
  }
  tfiles->Close();
  outtr->Write();
  histo_M_dr->Write();
  out->Close();
  //check cpu time (end)
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  std::cout << "CPU time used(sec): " << time_elapsed_ms/1000 << " sec\n";
  return 0;
}
