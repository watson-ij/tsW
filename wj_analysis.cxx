#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"

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

  auto out = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("tsW", "tsW");

// not complete yet
/*
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


  BranchI(nEvent);
  BranchI(nParticles);
*/

  histo1D(jets, gen_jets, particles, electrons, muons, out, trees, 5);
  histo1D(jets, gen_jets, particles, electrons, muons, out, trees, 3);

  tfiles->Close();
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


/*
  TH1D * h0 = new TH1D("number of events","number of events",1,0,1);
  TH1D * h1 = new TH1D("pT ratio Ks","pT ratio Ks",100,0,1);
  TH1D * h2 = new TH1D("pT ratio lamb","pT ratio lamb",100,0,1);
  TH1D * h3 = new TH1D("displaced length (rho)","displaced length (rho)",300,0,30);
  TH1D * h4 = new TH1D("displaced length (r)","displaced length (r)",300,0,30);
  TH1D * h5 = new TH1D("jet id","jet id",7,0,7);
  TH1D * h6 = new TH1D("lepton isolation","lepton isolation",200,0,2);
  TH1D * h7 = new TH1D("lepton in the jet","lepton in the jet",2,0,2);
  TH1D * h8 = new TH1D("pT ratio lepton","pT ratio lepton",100,0,1);
  TH1D * h9 = new TH1D("Energy ratio lepton","Energy ratio lepton",100,0,1);
  TH1D * h10 = new TH1D("dipion mass","dipion mass",1000,0,2000);
  TH1D * h11 = new TH1D("dipion delta r","dipion delta r",100,0,1);
*/

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

