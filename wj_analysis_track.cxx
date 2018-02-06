#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <list>
#include <ctime>
#include <stdio.h>

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

// initialization functions
const GenParticle* getLast(TClonesArray *particles, size_t iev, const GenParticle* p);
Double_t getdr(const GenParticle* p1, const GenParticle* p2);
std::vector<const GenParticle*> getMlist(TClonesArray * particles, size_t iev, const GenParticle* p);
void histo1D(TClonesArray * jets, TClonesArray * gen_jets, TClonesArray * particles, TClonesArray * electrons, TClonesArray * muons, TFile * out, TTree* trees, int jetpid);
void histo_numLepton(TClonesArray * electrons, TClonesArray * muons, TFile * out, TTree * trees);

// struct for selected leptons
struct lepton {
  TLorentzVector tlv;
  int charge;
  int pdgid;
};

int main(int argc, char* argv[])
{
  auto inf = std::string{argv[1]};
  auto outf = std::string{argv[2]};

  // check cpu time (start)
  std::clock_t c_start = std::clock();

  auto tfiles = TFile::Open(inf.c_str(), "READ");
  auto trees = (TTree*) tfiles->Get("Delphes");
  trees->SetBranchStatus("*", true);

  TClonesArray *jets = 0;
  trees->SetBranchAddress("Jet", &jets);
  TClonesArray *gen_jets = 0;
  trees->SetBranchAddress("GenJet", &gen_jets);
  TClonesArray *particles = 0;
  trees->SetBranchAddress("Particle", &particles);

  TClonesArray *electrons = 0;
  trees->SetBranchAddress("Electron", &electrons);
  TClonesArray *muons = 0;
  trees->SetBranchAddress("Muon", &muons);

  TClonesArray *tracks = 0;
  trees->SetBranchAddress("Track", &tracks);

  TClonesArray *missing = 0;
  trees->SetBranchAddress("MissingET", &missing);

  auto out = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("tsW", "tsW");

// not complete yet

#define Branch_(type, name, suffix) type name = 0; outtr->Branch(#name, &name, #name "/" #suffix);
#define BranchI(name) Branch_(Int_t, name, I)
#define BranchF(name) Branch_(Float_t, name, F)
#define BranchO(name) Branch_(Bool_t, name, O)
#define BranchA_(type, name, size, suffix) type name[size] = {0.}; outtr->Branch(#name, &name, #name"["#size"]/"#suffix);
#define BranchAI(name, size) BranchA_(Int_t, name, size, I);
#define BranchAF(name, size) BranchA_(Float_t, name, size, F);
#define BranchAO(name, size) BranchA_(Bool_t, name, size, O);
#define BranchVF(name) std::vector<float> name; outtr->Branch(#name, "vector<float>", &name);
#define BranchVI(name) std::vector<int> name; outtr->Branch(#name, "vector<int>", &name);
#define BranchVO(name) std::vector<bool> name; outtr->Branch(#name, "vector<bool>", &name);

#define BranchP_(type, branchname, name, suffix) type name = 0; TBranch *branchname =  outtr->Branch(#name, &name, #name "/" #suffix);
#define BranchPI(branchname,name) BranchP_(Int_t, branchname,name, I);
#define BranchPF(branchname,name) BranchP_(Float_t,branchname, name, F);
#define BranchPO(branchname,name) BranchP_(Bool_t,branchname, name, O);
#define BranchPA_(type, branchname, name, size, suffix) type name[size] = {0.}; TBranch *branchname = outtr->Branch(#name, &name, #name"["#size"]/"#suffix);
#define BranchPAI(branchname, name, size) BranchA_(Int_t, name, size, I);
#define BranchPAF(branchname, name, size) BranchA_(Float_t, name, size, F);
#define BranchPAO(branchname, name, size) BranchA_(Bool_t, name, size, O);
#define BranchPVF(branchname, name) std::vector<float> name; TBranch *branchname = outtr->Branch(#name, "vector<float>", &name);
#define BranchPVI(branchname, name) std::vector<int> name; TBranch *branchname = outtr->Branch(#name, "vector<int>", &name);
#define BranchPVO(branchname, name) std::vector<bool> name; TBranch *branchname = outtr->Branch(#name, "vector<bool>", &name);


  BranchI(nParticles);
  BranchF(dilepton_mass);
  BranchI(nJets);

  BranchI(nLepton);
  BranchI(dilepton_ch);

  BranchF(pion_track_X);
  BranchF(pion_track_Y);
  BranchF(pion_track_Z);
  BranchF(pion_track_R);
  BranchF(pion_track_path);
  BranchF(pion_track_Xd);
  BranchF(pion_track_Yd);
  BranchF(pion_track_Zd);
  BranchF(pion_track_Rd);
  BranchF(pion_track_Xout);
  BranchF(pion_track_Yout);
  BranchF(pion_track_Zout);
  BranchF(pion_track_Rout);

  BranchF(track_X);
  BranchF(track_Y);
  BranchF(track_Z);
  BranchF(track_R);
  BranchF(track_path);
  BranchF(track_Xd);
  BranchF(track_Yd);
  BranchF(track_Zd);
  BranchF(track_Rd);
  BranchF(track_Xout);
  BranchF(track_Yout);
  BranchF(track_Zout);
  BranchF(track_Rout);

  //histo_numLepton(electrons,muons,out,trees);

  TString cutflow_title = "cutflow" + inf;
  TH1F * cutflow = new TH1F("cutflow", cutflow_title, 6, 0, 6);
  TH2F * eff_err = new TH2F("eff_err", "efficiency error", 6, 0, 6, 100, 0.0, 0.004);
  auto hTrack = new TH2F("Track", "Track R vs Track Z", 1500,0,1500,7000,-3500,3500);
  auto hPion_Track = new TH2F("Pion_Track", "Pion Track R vs Z",1500,0,1500,7000,-3500,3500);
  auto hPion_Track_KS = new TH2F("Pion_Track_with_KS", "Pion Track taking KS as mother R vs Z",1500,0,1500,7000,-3500,3500);
  auto hPion_Track_KS_3d = new TH3F("Pion_Track_3d", "Pion Track R vs Z 3D", 300,-1500,1500,300,-1500,1500,700,-3500,3500);
  Double_t Nev = trees->GetEntries();

  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    trees->GetEntry(iev);

    cutflow->Fill(0);

    nParticles = 0;
    dilepton_mass = 0;
    nLepton = 0;
    dilepton_ch = 0;
    nJets = 0;

    pion_track_X = 0;
    pion_track_Y = 0;
    pion_track_Z = 0;
    pion_track_R = 0;
    pion_track_path = 0;
    pion_track_Xd = 0;
    pion_track_Yd = 0;
    pion_track_Zd = 0;
    pion_track_Rd = 0;
    pion_track_Xout = 0;
    pion_track_Yout = 0;
    pion_track_Zout = 0;
    pion_track_Rout = 0;

    track_X = 0;
    track_Y = 0;
    track_Z = 0;
    track_R = 0;
    track_path = 0;
    track_Xd = 0;
    track_Yd = 0;
    track_Zd = 0;
    track_Rd = 0;
    track_Xout = 0;
    track_Yout = 0;
    track_Zout = 0;
    track_Rout = 0;

    if (iev%10000 == 0 ) std::cout << "event check    iev    ----> " << iev << std::endl;
    nParticles = particles->GetEntries();

    // event selection (dilepton channel)
    std::vector<struct lepton> recolep;
    std::vector<Jet*> selectJets;
    std::vector<Jet*> selectBJets;

    recolep.clear();
    selectJets.clear();
    selectBJets.clear();

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

    nLepton = recolep.size(); 
    if (nLepton < 2 ) {
      continue;
    }
    if(recolep.size() > 2){
      sort(recolep.begin(), recolep.end(), [](struct lepton a, struct lepton b){return a.tlv.Pt() > b.tlv.Pt();});
      recolep.erase(recolep.begin()+2,recolep.end());
    }
    auto dilepton = recolep[0].tlv + recolep[1].tlv;
    dilepton_ch = recolep[0].pdgid + recolep[1].pdgid; // 22 -> ee , 24 -> emu , 26 -> mumu

    // step 1
    if(dilepton.M() < 20. || recolep[0].charge * recolep[1].charge > 0 ) {
      continue;
    }
    cutflow->Fill(1);

    // step2
    if ( (dilepton_ch != 24) && ( dilepton.M() > 76. && dilepton.M() < 106.) ) {
      continue;
    }
    cutflow->Fill(2);

    // step3
    auto miss = (MissingET*) missing->At(0);
    if ( (dilepton_ch != 24) && (miss->MET < 40.) ) {
      continue;
    }
    cutflow->Fill(3);

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
      selectJets.push_back(jet);
      nJets = nJets + 1;
    }
   if (selectJets.size() < 2) {
     continue;
   }
   cutflow->Fill(4);

    // step5
/*
    for (unsigned l = 0; l < selectJets.size(); ++l){
      if (selectJets[l]->BTag) {
        selectBJets.push_back(selectJets[l]);
      }
    }
    if (selectBJets.size() < 1) {
      continue;
    }
    cutflow->Fill(5);
*/
    dilepton_mass = dilepton.M();


    //for (size_t a = 0; a < particles->GetEntries(); ++a){
    for (size_t a = 0; a < tracks->GetEntries(); ++a){
      auto trackp = (Track*) tracks->At(a);
      auto gen_pa = (GenParticle*) trackp->Particle.GetObject();

      track_X = trackp->X;
      track_Y = trackp->Y;
      track_Z = trackp->Z;
      track_R = sqrt(trackp->X*trackp->X + trackp->Y*trackp->Y);
      track_path = trackp->L;
      track_Xd = trackp->Xd;
      track_Yd = trackp->Yd;
      track_Zd = trackp->Zd;
      track_Rd = sqrt(trackp->Xd*trackp->Xd + trackp->Yd*trackp->Yd);
      track_Xout = trackp->XOuter;
      track_Yout = trackp->YOuter;
      track_Zout = trackp->ZOuter;
      track_Rout = sqrt(trackp->XOuter*trackp->XOuter + trackp->YOuter*trackp->YOuter);

      hTrack->Fill(sqrt(trackp->X*trackp->X + trackp->Y*trackp->Y), trackp->Z);

      if(abs(trackp->PID) == 211) hPion_Track->Fill(sqrt(trackp->X*trackp->X + trackp->Y*trackp->Y), trackp->Z);

      if(!(gen_pa->M1 == -1)){
	auto gen_pa_M1 = (GenParticle*) particles->At(gen_pa->M1);
        if(abs(gen_pa_M1->PID) == 310 && abs(trackp->PID) == 211){
          std::cout << "track check : " << trackp->PID << " track's mother check : " << gen_pa_M1->PID << std::endl;
          pion_track_X = trackp->X;
          pion_track_Y = trackp->Y;
          pion_track_Z = trackp->Z;
          pion_track_R = sqrt(trackp->X*trackp->X + trackp->Y*trackp->Y);
          pion_track_path = trackp->L;
          pion_track_Xd = trackp->Xd;
          pion_track_Yd = trackp->Yd;
          pion_track_Zd = trackp->Zd;
          pion_track_Rd = sqrt(trackp->Xd*trackp->Xd + trackp->Yd*trackp->Yd);
          pion_track_Xout = trackp->XOuter;
          pion_track_Yout = trackp->YOuter;
          pion_track_Zout = trackp->ZOuter;
          pion_track_Rout = sqrt(trackp->XOuter*trackp->XOuter + trackp->YOuter*trackp->YOuter);

          hPion_Track_KS->Fill(sqrt(trackp->X*trackp->X + trackp->Y*trackp->Y), trackp->Z);
  	  hPion_Track_KS_3d->Fill(trackp->X, trackp->Y, trackp->Z);
        }
      }
    }
    outtr->Fill();  
  }
  Double_t P_no_step = cutflow->GetBinContent(1);
  Double_t P_step1 = cutflow->GetBinContent(2);
  Double_t P_step2 = cutflow->GetBinContent(3);
  Double_t P_step3 = cutflow->GetBinContent(4);
  Double_t P_step4 = cutflow->GetBinContent(5);
  //Double_t P_step5 = cutflow->GetBinContent(6);
  //eff_err->Fill(0.,(P_no_step*(Nev - P_no_step))/(Nev*Nev*Nev));
  eff_err->Fill(1.,sqrt((P_step1*(Nev - P_step1))/(Nev*Nev*Nev)));
  eff_err->Fill(2.,sqrt((P_step2*(Nev - P_step2))/(Nev*Nev*Nev)));
  eff_err->Fill(3.,sqrt((P_step3*(Nev - P_step3))/(Nev*Nev*Nev)));
  eff_err->Fill(4.,sqrt((P_step4*(Nev - P_step4))/(Nev*Nev*Nev)));

  std::cout << "check eff err for step 1 : " << sqrt((P_step1*(Nev - P_step1))/(Nev*Nev*Nev)) << std::endl;
  std::cout << "check eff err for step 2 : " << sqrt((P_step2*(Nev - P_step2))/(Nev*Nev*Nev)) << std::endl;
  std::cout << "check eff err for step 3 : " << sqrt((P_step3*(Nev - P_step3))/(Nev*Nev*Nev)) << std::endl;
  std::cout << "check eff err for step 4 : " << sqrt((P_step4*(Nev - P_step4))/(Nev*Nev*Nev)) << std::endl;

  histo1D(jets, gen_jets, particles, electrons, muons, out, trees, 3);
  histo1D(jets, gen_jets, particles, electrons, muons, out, trees, 5);

  cutflow->Write();
  eff_err->Write();
  hTrack->Write();
  hPion_Track->Write();
  hPion_Track_KS->Write();
  hPion_Track_KS_3d->Write();

  tfiles->Close();
  outtr->Write();

  out->Close();

  //check cpu time (end)
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  std::cout << "CPU time used(sec): " << time_elapsed_ms/1000 << " sec\n";
  return 0;
}

const GenParticle* getLast(TClonesArray * particles, size_t iev, const GenParticle* p){
  auto mom = p;
  while(true){
    auto dau = (const GenParticle*)particles->At(mom->D1);
    if( abs(p->PID) != abs(dau->PID) ) break;
    mom = dau;
  }
  return mom;
}

Double_t getdr(const GenParticle* p1, const GenParticle* p2){
  TLorentzVector p1_tlv;
  p1_tlv.SetPtEtaPhiE(p1->PT,p1->Eta,p1->Phi,p1->E);
  TLorentzVector p2_tlv;
  p2_tlv.SetPtEtaPhiE(p2->PT,p2->Eta,p2->Phi,p2->E);
  return p1_tlv.DeltaR(p2_tlv);
}

std::vector<const GenParticle*> getMlist(TClonesArray * particles, size_t iev, const GenParticle* p){
  std::vector<const GenParticle*> mlst;
  auto idx = p->M1;
  while(true){
    auto m = (const GenParticle*)particles->At(idx);
    mlst.push_back(m);
    idx = m->M1;
    if ( idx == -1) break;
  }
  return mlst;
}

void histo1D(TClonesArray* jets,TClonesArray * gen_jets,TClonesArray * particles, TClonesArray * electrons,TClonesArray * muons, TFile * out, TTree * trees, int jetpid){

  std::string jetpid_str = std::to_string(jetpid);
  TString h0char = "number of events when abs(jetpid) = " + jetpid_str;
  TString h1char = "pT ratio Ks when abs(jetpid) = " + jetpid_str;
  TString h2char = "pT ratio lamb when abs(jetpid) = " + jetpid_str;
  TString h3char = "displaced length (rho) when abs(jetpid) = " + jetpid_str;
  TString h4char = "displaced length (r) when abs(jetpid) = " + jetpid_str;
  TString h5char = "jet id when abs(jetpid) = " + jetpid_str;
  TString h6char = "lepton isolation when abs(jetpid) = " + jetpid_str;
  TString h7char = "lepton in the jet when abs(jetpid) = " + jetpid_str;
  TString h8char = "pT ratio lepton when abs(jetpid) = " + jetpid_str;
  TString h9char = "Energy ratio lepton when abs(jetpid) = " + jetpid_str;
  TString h10char = "dipion mass when abs(jetpid) = " + jetpid_str;
  TString h11char = "dipion delta r when abs(jetpid) = " + jetpid_str;

  TH1D * h0 = new TH1D( h0char, h0char,1,0,1);
  TH1D * h1 = new TH1D( h1char, h1char,100,0,1);
  TH1D * h2 = new TH1D( h2char, h2char,100,0,1);
  TH1D * h3 = new TH1D( h3char, h3char,300,0,30);
  TH1D * h4 = new TH1D( h4char, h4char,300,0,30);
  TH1D * h5 = new TH1D( h5char, h5char,7,0,7);
  TH1D * h6 = new TH1D( h6char, h6char,200,0,2);
  TH1D * h7 = new TH1D( h7char, h7char,2,0,2);
  TH1D * h8 = new TH1D( h8char, h8char,100,0,1);
  TH1D * h9 = new TH1D( h9char, h9char,100,0,1);
  TH1D * h10 = new TH1D( h10char, h10char,1000,0,2000);
  TH1D * h11 = new TH1D( h11char, h11char,100,0,1);

  for( size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%100 ==0) std::cout << "event check : " << iev << std::endl;
    h0->Fill(0.5);
    trees->GetEntry(iev);
    auto nParticle = particles->GetEntries();
    for( unsigned j = 0; j < particles->GetEntries(); ++j) {
      auto p = (const GenParticle*) particles->At(j);
      if (p->Status > 30) continue;
      if (abs(p->PID) != 6) continue;
      auto lastTop = getLast(particles,iev,p);
      auto jet = (const GenParticle*) particles->At(lastTop->D2);
      if (abs(jet->PID) != jetpid) continue;
      h5->Fill(abs(jet->PID));
      std::vector<const GenParticle*> pions;
      for( unsigned k = 0; k < particles->GetEntries(); ++k) {
        auto q = (const GenParticle*) particles->At(k);
        if (abs(q->PID) == 211) {}
        else if (abs(q->PID) == 321) {}
        else if (abs(q->PID) == 2212) {}
        else if (abs(q->PID) == 11 ) {}
        else if (abs(q->PID) == 13) {}
        else if (q->PID == 310) {}
        else if (q->PID == 3122) {}
        else continue;
        auto deltaR = getdr(jet, q);
        if (deltaR > 0.5) continue;
        if ( abs(q->PID) == 211 || abs(q->PID) == 321 || abs(q->PID) == 2212 ) {
          pions.push_back(q);
        }
        else if ( abs(q->PID) == 310) {
          h1->Fill((q->PT/jet->PT));
        }
        else if ( abs(q->PID) == 3122) {
          h2->Fill((q->PT/jet->PT));
        }
        else if ( abs(q->PID) == 11 || abs(q->PID) == 13) {
          h6->Fill(0);
          getMlist(particles, iev, q);
          bool isJetLepton = (std::find(getMlist(particles,iev,q).begin(), getMlist(particles,iev,q).end(), jet) != getMlist(particles,iev,q).end());
          h7->Fill(isJetLepton);
          if (isJetLepton) {
            h3->Fill(sqrt(q->X*q->X + q->Y*q->Y));
            h4->Fill(sqrt(q->X*q->X + q->Y*q->Y + q->Z*q->Z));
            h8->Fill(q->PT/jet->PT);
            h9->Fill(q->E/jet->E);
          }
        }
        else continue;
      }
      if (pions.size() <= 1) continue;
      for(size_t i=0; i<pions.size(); ++i) {
        for(size_t j = i+1; j<pions.size(); ++j){
          if((pions[i]->PID * pions[j]->PID) > 0) continue;
          TLorentzVector pion1_tlv;
          pion1_tlv.SetPtEtaPhiE(pions[i]->PT, pions[i]->Eta, pions[i]->Phi, pions[i]->E);
          TLorentzVector pion2_tlv;
          pion2_tlv.SetPtEtaPhiE(pions[j]->PT, pions[j]->Eta, pions[j]->Phi, pions[j]->E);
          auto dipion_tlv = pion1_tlv + pion2_tlv;
          h10->Fill(dipion_tlv.M()*1000);
          h11->Fill(pion1_tlv.DeltaR(pion2_tlv));
        }
      }
    }
  }
  h0->Write();
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h6->Write();
  h7->Write();
  h8->Write();
  h9->Write();
  h10->Write();
  h11->Write();

  return;
}

void histo_numLepton(TClonesArray * electrons, TClonesArray * muons, TFile * out, TTree * trees){
  TH1D * histo_nLepton = new TH1D("nLepton without ev_selection", "nLeptons without event selection", 10, 0, 10);
  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%1000 == 0) std::cout << "Filling histogram for number of leptons..." << iev << std::endl;
    trees->GetEntry(iev);
    auto nLepton = muons->GetEntries() + electrons->GetEntries();
    histo_nLepton->Fill(nLepton);
  }
  histo_nLepton->Write();
  return;
}

