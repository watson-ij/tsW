#include <stdio.h>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "wj_analysis.h"

int main(int argc, char* argv[])
{
  auto inf = std::string{argv[1]};
  auto outf = std::string{argv[2]};

  // check cpu time (start)
  std::clock_t c_start = std::clock();

  //read input file
  TClonesArray *gen_jets = 0, *particles = 0;
  TClonesArray *electrons = 0, *muons = 0, *jets = 0, *missingET = 0;

  auto tfiles = TFile::Open(inf.c_str(), "READ");
  auto trees = (TTree*) tfiles->Get("Delphes");
  trees->SetBranchStatus("*", true);
  trees->SetBranchAddress("GenJet",      &gen_jets);
  trees->SetBranchAddress("Particle",    &particles);
  trees->SetBranchAddress("Electron",    &electrons);
  trees->SetBranchAddress("Muon",        &muons);
  trees->SetBranchAddress("Jet",         &jets);
  trees->SetBranchAddress("MissingET",   &missingET);

  //make out file
  auto out = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("tsW", "tsW");

  defBranchFucns();
  BranchF(dilepton_mass);

  BranchI(dilepton_ch); BranchI(step);
  BranchO(step1); BranchO(step2); BranchO(step3); BranchO(step4); BranchO(step5);

  BranchO(gen_step0);
  //BranchLV(jet_tlv); BranchLV(Ks_tlv);
  BranchVF(jet_energy); BranchVF(jet_pt); BranchVI(jet_pid);
  BranchVF(Ks_energy); BranchVF(Ks_pt); BranchVI(nKs);

  //Event Loop Start!
  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%1000 == 0 ) std::cout << "event check    iev    ----> " << iev << std::endl;
    trees->GetEntry(iev);

    //init values
    dilepton_mass = -99;
    dilepton_ch = 0; step = 0;
    step1 = false; step2 = false; step3 = false; step4 = false; step5 = false;

    gen_step0 = false;
    //jet_tlv = TLorentzVector(); Ks_tlv = TLorentzVector();
    jet_energy.clear(); jet_pt.clear(); jet_pid.clear();
    Ks_energy.clear(); Ks_pt.clear(); nKs.clear();

    // event selection (dilepton channel)
    std::vector<struct lepton> recolep;
    std::vector<Jet*> selectedJets;

    //muon selection
    for (unsigned i = 0; i < muons->GetEntries(); ++i){
      auto mu = (Muon*) muons->At(i);
      if (abs(mu->Eta) > 2.4) continue;
      if (mu->PT < 20.) continue;
      TLorentzVector mu_tlv;
      mu_tlv = mu->P4();

      struct lepton selmuons;
      selmuons.tlv = mu_tlv;
      selmuons.charge = mu->Charge;
      selmuons.pdgid = 13;

      recolep.push_back(selmuons);
    }
    //electron selection
    for (unsigned j = 0; j < electrons->GetEntries(); ++j){
      auto elec = (Electron*) electrons->At(j);
      if (abs(elec->Eta) > 2.4) continue;
      if (elec->PT < 20.) continue;
      TLorentzVector elec_tlv;
      elec_tlv = elec->P4();

      struct lepton selelecs;
      selelecs.tlv = elec_tlv;
      selelecs.charge = elec->Charge;
      selelecs.pdgid = 11;

      recolep.push_back(selelecs);
    }

    if (recolep.size() >= 2 ) {
      sort(recolep.begin(), recolep.end(), [](struct lepton a, struct lepton b){return a.tlv.Pt() > b.tlv.Pt();});
      recolep.erase(recolep.begin()+2,recolep.end());

      auto dilepton = recolep[0].tlv + recolep[1].tlv;
      dilepton_ch = recolep[0].pdgid + recolep[1].pdgid; // 22 -> ee , 24 -> emu , 26 -> mumu
      dilepton_mass = dilepton.M();

      // step 1
      if(dilepton.M() > 20. && recolep[0].charge * recolep[1].charge < 0 ) {
        step1 = true;
        step++;
      }

      // step2
      if ( (dilepton_ch == 24) || ( dilepton.M() < 76. || dilepton.M() > 106.) ) {
        step2 = true;
        if (step == 1) step++;
      }

      // step3
      auto met = ((MissingET*) missingET->At(0))->MET;
      if ( (dilepton_ch == 24) || (met > 40.) ) {
        step3 = true;
        if (step == 2) step++;
      }

      // step4  
      for ( unsigned k = 0; k < jets->GetEntries(); ++k){
        auto jet = (Jet*) jets->At(k);
        TLorentzVector jet_tlv;
        jet_tlv = jet->P4();
        bool hasOverLap = false;
        if (jet->PT < 30.) continue;
        if (abs(jet->Eta) > 2.4) continue;
        for( unsigned kk=0; kk < recolep.size(); ++kk){
	  if ( jet_tlv.DeltaR(recolep[kk].tlv) < 0.4) hasOverLap = true;
        }
        if (hasOverLap) continue;
        selectedJets.push_back(jet);
      }

      if (selectedJets.size() > 1) {
        step4 = true;
        if ( step == 3) step++;
      }

      // step5
      std::vector<Jet*> selectedBJets;
      for (auto jet : selectedJets){
        if (jet->BTag) selectedBJets.push_back(jet);
      }
      
      if (selectedBJets.size() > 0){
	    step5 = true;
	    if ( step == 4) step++;
      }
    }

    // collect gen objects
    std::vector<GenParticle> genTops;
    std::vector<GenParticle> genJets;
    std::vector<struct lepton> genLeps;
    for (unsigned i = 0; i < particles->GetEntries(); ++ i){
      auto p = (const GenParticle*) particles->At(i);
      if (p->Status > 30) continue;
      if (abs(p->PID) != 6) continue;

      auto top = getLast(particles,iev,p);
      auto jet = (const GenParticle*) particles->At(top->D2);
      auto Wboson = ( const GenParticle*) particles->At(top->D1);
      genTops.push_back(*top);
      genJets.push_back(*jet);

      auto lastBoson = getLast(particles, iev, Wboson);
      auto wDau1 = (const GenParticle*) particles->At(lastBoson->D1);
      auto wDau2 = (const GenParticle*) particles->At(lastBoson->D2);
      struct lepton lep;
      lep.tlv = wDau1->P4();
      lep.charge = wDau1->Charge;
      lep.pdgid = wDau1->PID;
      genLeps.push_back(lep);

      // what about GenJet???
      auto jet_tlv = jet->P4();
      jet_energy.push_back(jet_tlv.E());
      jet_pt.push_back(jet_tlv.Pt());
      jet_pid.push_back(jet->PID);

      // collect jet constitues
      std::vector<GenParticle> jetConstitues;
      std::vector<int> pdgList = {211, 321, 2212, 310, 3122, 11, 13};
      for (unsigned i2 = 0; i2 < particles->GetEntries(); ++ i2){
        auto p2 = (const GenParticle*) particles->At(i2);

        if ( !(std::find(pdgList.begin(), pdgList.end(), abs(p2->PID)) != pdgList.end()) ) continue;
	    if ( jet_tlv.DeltaR(p2->P4()) > 0.5) continue;
        jetConstitues.push_back(*p2);
      }

      std::vector<GenParticle> KsInjet;
      std::vector<GenParticle> lambdaInjet;
      std::vector<GenParticle> pionInjet;
      for (auto c : jetConstitues){
        if (abs(c.PID) == 310)  KsInjet.push_back(c);
        if (abs(c.PID) == 3122) lambdaInjet.push_back(c);
        if (abs(c.PID) == 211)  pionInjet.push_back(c);
      }

      // save highest pT Ks only
      //sort(KsInjet.begin(), KsInjet.end(), [](std::vector<GenParticle> a, std::vector<GenParticle> b){return a.PT > b.PT;});
      nKs.push_back(KsInjet.size());
      if (KsInjet.size() >0){
        Ks_energy.push_back(KsInjet[0].E);
        Ks_pt.push_back(KsInjet[0].PT);
      }

    }

    // accepted leptons
    std::vector<struct lepton> selectedLeps;
    for (auto genLep : genLeps){
      if (abs(genLep.tlv.Eta()) > 2.4) continue;
      if (genLep.tlv.Pt() < 20.) continue;
      selectedLeps.push_back(genLep);
    }
    if (selectedLeps.size() >= 2) gen_step0 = true;

    outtr->Fill();  
  }

  //br->Write();
  tfiles->Close();
  outtr->Write();


  Int_t st;
  outtr->SetBranchAddress("step", &st);
  TH1F * cutflow = new TH1F("cutflow", "cutflow", 6,0,6);
  for ( size_t i = 0; i < outtr->GetEntries(); ++i){
    outtr->GetEntry(i);
    for( size_t j = 0; j < st; j++){
      cutflow->Fill(j);
    }
  }

  cutflow->Write();
  out->Close();

  //check cpu time (end)
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  std::cout << "CPU time used(sec): " << time_elapsed_ms/1000 << " sec\n";
  return 0;
}

