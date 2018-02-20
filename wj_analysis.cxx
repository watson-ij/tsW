#include <stdio.h>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
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
  auto outtr = new TTree("tsW", "tsW");

  defBranchFucns(outtr);

  TString cutflow_title = "cutflow" + inf;
  TH1F * cutflow = new TH1F("cutflow", cutflow_title, 7,-1,6); // -1 : all events, 0 : events after lepton selection, 1~ : events after step
  TH1F * histo_dihadron_S = new TH1F("dihadron_mass_S", "dihadron_mass_S", 300, 0, 3000);
  TH1F * histo_dihadron_B = new TH1F("dihadron_mass_B", "dihadron_mass_B", 300, 0, 3000);

  TH1F * histo_x_KS = new TH1F("hist_x_KS", "hist_x_KS", 100, 0, 1);
  TH1F * histo_x_lamb = new TH1F("hist_x_lamb", "hist_x_lamb", 100, 0, 1);

  TH1F * histo_rho_KS = new TH1F("hist_rho_KS", "hist_rho_KS", 1000, 0, 300);
  TH1F * histo_rho_lamb = new TH1F("hist_rho_lamb", "hist_rho_lamb", 1800, 0, 18000);

  TH1F * histo_d_KS = new TH1F("hist_d_KS", "hist_d_KS", 1000, 0, 20);
  TH1F * histo_d_lamb = new TH1F("hist_d_lamb", "hist_d_lamb", 1000, 0, 100);

  //Event Loop Start!
  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%1000 == 0 ) std::cout << "event check    iev    ----> " << iev << std::endl;
    trees->GetEntry(iev);
    initValues();

    //recoParticle(cutflow);
    //genParticle(histo_dihadron_S, histo_dihadron_B);
    Finder(outtr);

    /*
    // pion track test
    for (unsigned i = 0; i < tracks->GetEntries(); ++ i){
      auto track = (Track*) tracks->At(i);
      if (abs(track->PID) == 321){
        auto genTrack = (GenParticle*) track->Particle.GetObject();
        auto genTrack_mom = (const GenParticle*) particles->At(genTrack->M1);
        std::cout << genTrack_mom->PID << std::endl;
        float pionR = sqrt(pow(track->X,2)+pow(track->Y,2)+pow(track->Z,2));
        float pionRd = sqrt(pow(track->Xd,2)+pow(track->Yd,2)+pow(track->Zd,2));
        float pionROuter = sqrt(pow(track->XOuter,2)+pow(track->YOuter,2)+pow(track->ZOuter,2));
      }
    }
    */
    if(!jetFinder)outtr->Fill();
  }

  std::cout << count << std::endl;

  tfiles->Close();

  outtr->Write();

  outtr->SetBranchAddress("x_KS", &x_KS);
  outtr->SetBranchAddress("x_lamb", &x_lamb);
  outtr->SetBranchAddress("rho_KS", &rho_KS);
  outtr->SetBranchAddress("rho_lamb", &rho_lamb);
  outtr->SetBranchAddress("d_KS", &d_KS);
  outtr->SetBranchAddress("d_lamb", &d_lamb);

  for (size_t ent = 0; ent < outtr->GetEntries(); ++ent){
    outtr->GetEntry(ent);

    histo_x_KS->Fill(x_KS);
    histo_rho_KS->Fill(rho_KS);
    histo_d_KS->Fill(d_KS);

    histo_x_lamb->Fill(x_lamb);
    histo_rho_lamb->Fill(rho_lamb);
    histo_d_lamb->Fill(d_lamb);
  }

  histo_x_KS->Write();
  histo_rho_KS->Write();
  histo_d_KS->Write();

  histo_x_lamb->Write();
  histo_rho_lamb->Write();
  histo_d_lamb->Write();

  cutflow->Write();

  histo_dihadron_S->Write();
  histo_dihadron_B->Write();

  out->Close();

  //check cpu time (end)
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  std::cout << "CPU time used(sec): " << time_elapsed_ms/1000 << " sec\n";
  return 0;
}

void initValues(){
    dilepton_mass = -99; dilepton_ch = 0; 

    x_KS = -1; x_KS_S = -1; x_KS_B = -1; rho_KS = -99; rho_KS_S = -99; rho_KS_B = -99; d_KS = -99; d_KS_S = -99; d_KS_B = -99;
    x_lamb = -1; x_lamb_S = -1; x_lamb_B = -1; rho_lamb = -999; rho_lamb_S = -999; rho_lamb_B = -999; d_lamb = -99; d_lamb_S = -99; d_lamb_B = -99;

    significance_S = -99; significance_B = -99;

    channel = 0;

    gen_step0 = false; gen_step1 = false;

    jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_energy.clear();
    jet_pid.clear();

    kshortsInjet_pt.clear(); kshortsInjet_eta.clear(); kshortsInjet_phi.clear(); kshortsInjet_energy.clear(); kshortsInjet_R.clear(); kshortsInjet_outR.clear(); kshortsInjet_rho.clear(); kshortsInjet_d.clear();
    lambdasInjet_pt.clear(); lambdasInjet_eta.clear(); lambdasInjet_phi.clear(); lambdasInjet_energy.clear(); lambdasInjet_R.clear(); lambdasInjet_outR.clear(); lambdasInjet_rho.clear(); lambdasInjet_d.clear();
    leptonsInjet_pt.clear(); leptonsInjet_eta.clear(); leptonsInjet_phi.clear(); leptonsInjet_energy.clear(); leptonsInjet_R.clear();
    nkshortsInjet.clear(); nlambdasInjet.clear(); nleptonsInjet.clear();

    jet1_diHadron_mass.clear(); jet2_diHadron_mass.clear();
}

void recoParticle(TH1F * cutflow){
    cutflow->Fill(-1);
    // object selection
    std::vector<struct Lepton> recolep;
    for (unsigned i = 0; i < muons->GetEntries(); ++i){
      auto mu = (Muon*) muons->At(i);
      if (abs(mu->Eta) > 2.4) continue;
      if (mu->PT < 20.) continue;
      recolep.push_back(toLepton(mu));
    }
    for (unsigned j = 0; j < electrons->GetEntries(); ++j){
      auto elec = (Electron*) electrons->At(j);
      if (abs(elec->Eta) > 2.4) continue;
      if (elec->PT < 20.) continue;
      recolep.push_back(toLepton(elec));
    }

    // pick highest pt lep pair
    if (recolep.size() < 2 ) return;

    sort(recolep.begin(), recolep.end(), [](struct Lepton a, struct Lepton b){return a.tlv.Pt() > b.tlv.Pt();});
    recolep.erase(recolep.begin()+2,recolep.end());

    cutflow->Fill(0);

    auto dilepton = recolep[0].tlv + recolep[1].tlv;
    dilepton_ch = abs(recolep[0].pdgid) + abs(recolep[1].pdgid); // 22 -> ee , 24 -> emu , 26 -> mumu
    dilepton_mass = dilepton.M();

    // step 1
    if(dilepton.M() < 20. || recolep[0].charge * recolep[1].charge > 0 ) return;
    cutflow->Fill(1);

    // step2
    if ( (dilepton_ch != 24) && ( dilepton.M() > 76. && dilepton.M() < 106.) ) return;
    cutflow->Fill(2);

    // step3
    auto met = ((MissingET*) missingET->At(0))->MET;
    if ( (dilepton_ch != 24) && (met < 40.) ) return;
    cutflow->Fill(3);

    // step4  
    std::vector<Jet*> selectedJets;
    for ( unsigned k = 0; k < jets->GetEntries(); ++k){
      auto jet = (Jet*) jets->At(k);
      if (jet->PT < 30.) continue;
      if (abs(jet->Eta) > 2.4) continue;
      bool hasOverlap = false;
      for (auto lep : recolep) { if ((jet->P4()).DeltaR(lep.tlv) < 0.4) hasOverlap = true; }
      if (hasOverlap) continue;
      selectedJets.push_back(jet);
    }
    if (selectedJets.size() < 2) return;
    cutflow->Fill(4);

    // step5
    std::vector<Jet*> selectedBJets;
    for (auto jet : selectedJets){
      if (jet->BTag) selectedBJets.push_back(jet);
    }
    if (selectedBJets.size() < 1) return;
    cutflow->Fill(5);

}

void genParticle(TH1F * histo_dihadron_S, TH1F * histo_dihadron_B){
    std::vector<const GenParticle*> genTops;
    std::vector<const GenParticle*> genJets;
    std::vector<struct Lepton> genLeps;
    std::vector<std::vector<float> > jet_diHadron_mass;

    for (unsigned i = 0; i < particles->GetEntries(); ++ i){
      auto p = (const GenParticle*) particles->At(i);
      if (p->Status > 30) continue;
      if (abs(p->PID) != 6) continue;
      auto top = getLast(particles,p);
      genTops.push_back(top);

      auto jet = (const GenParticle*) particles->At(top->D2);
      auto Wboson = (const GenParticle*) particles->At(top->D1);
      auto lastBoson = getLast(particles, Wboson);
      genJets.push_back(jet);

      auto lep_tmp = (const GenParticle*) particles->At(lastBoson->D1);
      auto neu = (const GenParticle*) particles->At(lastBoson->D2);

      struct Lepton lep = toLepton(lep_tmp);
      genLeps.push_back(lep);

      if (genLeps.size() == 2){
        if (abs(genLeps[0].pdgid) == 11 ||abs(genLeps[0].pdgid) == 13 ||abs(genLeps[0].pdgid) == 15){
          if (abs(genLeps[1].pdgid) == 11 ||abs(genLeps[1].pdgid) == 13 ||abs(genLeps[1].pdgid) == 15){
            channel = 2;
          }
        }
        else channel = 1;
      }


      jet_pid.push_back(jet->PID);
      auto jet_tlv = jet->P4();
      jet_pt.push_back(jet_tlv.Pt());
      jet_eta.push_back(jet_tlv.Eta());
      jet_phi.push_back(jet_tlv.Phi());
      jet_energy.push_back(jet_tlv.E());

      // collect jet constitues
      std::vector<GenParticle> jetConstitues;
      std::vector<int> pdgList = {211, 321, 2212, 310, 3122, 11, 13};
      for (unsigned i2 = 0; i2 < particles->GetEntries(); ++ i2){
        auto p2 = (const GenParticle*) particles->At(i2);

        if ( !(std::find(pdgList.begin(), pdgList.end(), abs(p2->PID)) != pdgList.end()) ) continue;
	    if ( jet_tlv.DeltaR(p2->P4()) > 0.5) continue;
        bool fromJet = false;
        std::vector<const GenParticle*> mlist = getMlist(particles, p2);
        for (auto m : mlist){ if (m == jet) fromJet = true; }
        if (!fromJet) continue;
        
        jetConstitues.push_back(*p2);
      }
      /*
      // collect jet constitues
      std::vector<GenParticle> jetConstitues;
      std::vector<int> pdgList = {211, 321, 2212, 310, 3122, 11, 13};
      for ( unsigned k = 0; k < gen_jets->GetEntries(); ++k){
        auto genjet = (Jet*) gen_jets->At(k);
        auto test = (genjet->Constituents).At(0);
        std::cout << test->PID << std::endl;
        for (c : jet.Constituents){
          if ( !(std::find(pdgList.begin(), pdgList.end(), abs(p2->PID)) != pdgList.end()) ) continue;
          jetConstitues.push_back(*p2);
        }
      }
      */

      std::vector<GenParticle> kshortsInjet;
      std::vector<GenParticle> lambdasInjet;
      std::vector<GenParticle> hadronsInjet;
      std::vector<GenParticle> leptonsInjet;
      for (auto c : jetConstitues){
        int absPid = abs(c.PID);
        if (absPid == 310)  {
/*
          auto dau1 = (GenParticle*) particles->At(c.D1);
          auto dau2 = (GenParticle*) particles->At(c.D2);

          std::cout << "[1] D1 check      : " << dau1->PID << " , ( " << dau1->Px << " , " << dau1->Py << " , " << dau1->Pz << " ) " << " , P : " << sqrt(dau1->Px*dau1->Px + dau1->Py*dau1->Py + dau1->Pz*dau1->Pz ) << " , " << dau1->P << " poistion : ( " << dau1->X << " , " << dau1->Y << " , " << dau1->Z << " ) " << std::endl;
          std::cout << "[2] D2 check      : " << dau2->PID << " , ( " << dau2->Px << " , " << dau2->Py << " , " << dau2->Pz << " ) " << " , P : " << sqrt(dau2->Px*dau2->Px + dau2->Py*dau2->Py + dau2->Pz*dau2->Pz ) << " , " << dau2->P << " poistion : ( " << dau2->X << " , " << dau2->Y << " , " << dau2->Z << " ) "  << std::endl;

	  std::cout << "[3] KS check      : " << c.PID << " , ( " << c.Px << " , " << c.Py << " , " << c.Pz << " ) " << " , P : " << sqrt(c.Px*c.Px + c.Py*c.Py + c.Pz*c.Pz ) << " , " << c.P << " poistion : ( " << c.X << " , " << c.Y << " , " << c.Z << " ) "  << std::endl;
*/
	  kshortsInjet.push_back(c);
	}
        if (absPid == 3122) {
/*
          auto dau1 = (GenParticle*) particles->At(c.D1);
          auto dau2 = (GenParticle*) particles->At(c.D2);

          std::cout << "[1] D1 check      : " << dau1->PID << " , ( " << dau1->Px << " , " << dau1->Py << " , " << dau1->Pz << " ) " << " , P : " << sqrt(dau1->Px*dau1->Px + dau1->Py*dau1->Py + dau1->Pz*dau1->Pz ) << " , " << dau1->P << " poistion : ( " << dau1->X << " , " << dau1->Y << " , " << dau1->Z << " ) " << std::endl;
          std::cout << "[2] D2 check      : " << dau2->PID << " , ( " << dau2->Px << " , " << dau2->Py << " , " << dau2->Pz << " ) " << " , P : " << sqrt(dau2->Px*dau2->Px + dau2->Py*dau2->Py + dau2->Pz*dau2->Pz ) << " , " << dau2->P << " poistion : ( " << dau2->X << " , " << dau2->Y << " , " << dau2->Z << " ) "  << std::endl;

          std::cout << "[3] lambda check      : " << c.PID << " , ( " << c.Px << " , " << c.Py << " , " << c.Pz << " ) " << " , P : " << sqrt(c.Px*c.Px + c.Py*c.Py + c.Pz*c.Pz ) << " , " << c.P << " poistion : ( " << c.X << " , " << c.Y << " , " << c.Z << " ) "  << std::endl;
*/
	  lambdasInjet.push_back(c);
	}
        if (absPid == 2212 || absPid == 321 || absPid == 211)  hadronsInjet.push_back(c);
        if (absPid == 11 || absPid == 13)  leptonsInjet.push_back(c);
      }

      // save highest pT Ks only
      nkshortsInjet.push_back(kshortsInjet.size());
      if (kshortsInjet.size() >0){
        sort(kshortsInjet.begin(), kshortsInjet.end(), [](GenParticle a, GenParticle b){return a.PT > b.PT;});

        kshortsInjet_pt.push_back(kshortsInjet[0].PT);
        kshortsInjet_eta.push_back(kshortsInjet[0].Eta);
        kshortsInjet_phi.push_back(kshortsInjet[0].Phi);
        kshortsInjet_energy.push_back(kshortsInjet[0].E);
        auto R = sqrt(pow(kshortsInjet[0].X,2)+pow(kshortsInjet[0].Y,2)+pow(kshortsInjet[0].Z,2));
        kshortsInjet_R.push_back(R);

        auto kshortDau = (const GenParticle*) particles->At(kshortsInjet[0].D1);
        auto outR = sqrt(pow(kshortDau->X,2)+pow(kshortDau->Y,2)+pow(kshortDau->Z,2));
        kshortsInjet_outR.push_back(outR);

	auto rho = sqrt(kshortsInjet[0].X*kshortsInjet[0].X + kshortsInjet[0].Y*kshortsInjet[0].Y);
	kshortsInjet_rho.push_back(rho);

	auto pion = (GenParticle*) particles->At(kshortsInjet[0].D1);
	std::vector<Double_t> PVtoKS = { pion->X, pion->Y, pion->Z };
	std::vector<Double_t> momentum_unit = { (kshortsInjet[0].Px/sqrt(kshortsInjet[0].Px*kshortsInjet[0].Px + kshortsInjet[0].Py*kshortsInjet[0].Py + kshortsInjet[0].Pz*kshortsInjet[0].Pz )), (kshortsInjet[0].Py/sqrt(kshortsInjet[0].Px*kshortsInjet[0].Px + kshortsInjet[0].Py*kshortsInjet[0].Py + kshortsInjet[0].Pz*kshortsInjet[0].Pz )), (kshortsInjet[0].Pz/sqrt(kshortsInjet[0].Px*kshortsInjet[0].Px + kshortsInjet[0].Py*kshortsInjet[0].Py + kshortsInjet[0].Pz*kshortsInjet[0].Pz )) };
	auto cross = cross3D(PVtoKS, momentum_unit);
	auto d = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
	kshortsInjet_d.push_back(d);
      }
      else {
        kshortsInjet_pt.push_back(-99);
        kshortsInjet_eta.push_back(-99);
        kshortsInjet_phi.push_back(-99);
        kshortsInjet_energy.push_back(-99);
        kshortsInjet_R.push_back(-99);
        kshortsInjet_outR.push_back(-99);
	kshortsInjet_rho.push_back(-99);
	kshortsInjet_d.push_back(-99);
      }

      if(nkshortsInjet.size() == 2){
	if(abs(jet->PID) == 3){
	  if(nkshortsInjet[1] > 0){
	    x_KS_S = (kshortsInjet_pt[1]/jet_pt[1]);
	    rho_KS_S = kshortsInjet_rho[1];
	    d_KS_S = kshortsInjet_d[1];
            if (x_KS_S < 0) x_KS_S = -1;
            if (rho_KS_S < 0) rho_KS_S = -99;
            if (d_KS_S < 0) d_KS_S = -99;
	  }
	  if(nkshortsInjet[0] > 0){
	    x_KS_B = (kshortsInjet_pt[0]/jet_pt[0]);
            rho_KS_B = kshortsInjet_rho[0];
            d_KS_B = kshortsInjet_d[0];
            if (x_KS_B < 0) x_KS_B = -1;
            if (rho_KS_B < 0) rho_KS_B = -99;
            if (d_KS_B < 0) d_KS_B = -99;
	  } 
	}
        if(abs(jet->PID) == 5){
          if(nkshortsInjet[1] > 0){
            x_KS_B = (kshortsInjet_pt[1]/jet_pt[1]);
            rho_KS_B = kshortsInjet_rho[1];
            d_KS_B = kshortsInjet_d[1];
            if (x_KS_B < 0) x_KS_B = -1;
            if (rho_KS_B < 0) rho_KS_B = -99;
            if (d_KS_B < 0) d_KS_B = -99;
          }
          if(nkshortsInjet[0] > 0){
            x_KS_S = (kshortsInjet_pt[0]/jet_pt[0]);
            rho_KS_S = kshortsInjet_rho[0];
            d_KS_S = kshortsInjet_d[0];
	    if (x_KS_S < 0) x_KS_S = -1;
            if (rho_KS_S < 0) rho_KS_S = -99;
            if (d_KS_S < 0) d_KS_S = -99;
          }
        }
	if(x_KS_S > x_KS_B) x_KS = x_KS_S;
	else x_KS = x_KS_B;
        if(rho_KS_S > rho_KS_B) rho_KS = rho_KS_S;
        else rho_KS = rho_KS_B;
	if(d_KS_S > d_KS_B) d_KS = d_KS_S;
	else d_KS = d_KS_B;
      } 

      // save highest pT lambda only
      nlambdasInjet.push_back(lambdasInjet.size());
      //std::cout << "size check : " << nlambdasInjet.size() << " , jetpid : " << jet->PID << std::endl;
      //if(nlambdasInjet.size() == 2) std::cout << "nlambdas[0] : " << nlambdasInjet[0] << " , nlambdas[1] : " << nlambdasInjet[1] << std::endl;
      if (lambdasInjet.size() >0){
        sort(lambdasInjet.begin(), lambdasInjet.end(), [](GenParticle a, GenParticle b){return a.PT > b.PT;});

        lambdasInjet_pt.push_back(lambdasInjet[0].PT);
        lambdasInjet_eta.push_back(lambdasInjet[0].Eta);
        lambdasInjet_phi.push_back(lambdasInjet[0].Phi);
        lambdasInjet_energy.push_back(lambdasInjet[0].E);

        auto R = sqrt(pow(lambdasInjet[0].X,2)+pow(lambdasInjet[0].Y,2)+pow(lambdasInjet[0].Z,2));
        lambdasInjet_R.push_back(R);

        auto lambdaDau = (const GenParticle*) particles->At(lambdasInjet[0].D1);
        auto outR = sqrt(pow(lambdaDau->X,2)+pow(lambdaDau->Y,2)+pow(lambdaDau->Z,2));
        lambdasInjet_outR.push_back(outR);

        auto rho = sqrt(lambdasInjet[0].X*lambdasInjet[0].X + lambdasInjet[0].Y*lambdasInjet[0].Y);
        lambdasInjet_rho.push_back(rho);

        auto pion = (GenParticle*) particles->At(lambdasInjet[0].D1);
        std::vector<Double_t> PVtoKS = { pion->X, pion->Y, pion->Z };
        std::vector<Double_t> momentum_unit = { (lambdasInjet[0].Px/sqrt(lambdasInjet[0].Px*lambdasInjet[0].Px + lambdasInjet[0].Py*lambdasInjet[0].Py + lambdasInjet[0].Pz*lambdasInjet[0].Pz )), (lambdasInjet[0].Py/sqrt(lambdasInjet[0].Px*lambdasInjet[0].Px + lambdasInjet[0].Py*lambdasInjet[0].Py + lambdasInjet[0].Pz*lambdasInjet[0].Pz )), (lambdasInjet[0].Pz/sqrt(lambdasInjet[0].Px*lambdasInjet[0].Px + lambdasInjet[0].Py*lambdasInjet[0].Py + lambdasInjet[0].Pz*lambdasInjet[0].Pz )) };
        auto cross = cross3D(PVtoKS, momentum_unit);
        auto d = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
        lambdasInjet_d.push_back(d);

      }
      else {
        lambdasInjet_pt.push_back(-99);
        lambdasInjet_eta.push_back(-99);
        lambdasInjet_phi.push_back(-99);
        lambdasInjet_energy.push_back(-99);
        lambdasInjet_R.push_back(-99);
        lambdasInjet_outR.push_back(-99);
        lambdasInjet_rho.push_back(-999);
	lambdasInjet_d.push_back(-99);
      }

      if(nlambdasInjet.size() == 2){
        if(abs(jet->PID) == 3){
          if(nlambdasInjet[1] > 0){
            x_lamb_S = (lambdasInjet_pt[1]/jet_pt[1]);
	    rho_lamb_S = lambdasInjet_rho[1];
            d_lamb_S = lambdasInjet_d[1];
            if (x_lamb_S < 0) x_lamb_S = -1;
            if (rho_lamb_S < 0) rho_lamb_S = -999;
            if (d_lamb_S < 0) d_lamb_S = -99;
          }
          if(nkshortsInjet[0] > 0){
            x_lamb_B = (lambdasInjet_pt[0]/jet_pt[0]);
            rho_lamb_B = lambdasInjet_rho[0];
            d_lamb_B = lambdasInjet_d[0];
            if (x_lamb_B < 0) x_lamb_B = -1;
            if (rho_lamb_B < 0) rho_lamb_B = -999;
            if (d_lamb_B < 0) d_lamb_B = -99;
          }
        }
        if(abs(jet->PID) == 5){
          if(nlambdasInjet[1] > 0){
            x_lamb_B = (lambdasInjet_pt[1]/jet_pt[1]);
            rho_lamb_B = lambdasInjet_rho[1]; 
            d_lamb_B = lambdasInjet_d[1];
            if (x_lamb_B < 0) x_lamb_B = -1;
            if (rho_lamb_B < 0) rho_lamb_B = -999;
            if (d_lamb_B < 0) d_lamb_B = -99;
          }
          if(nkshortsInjet[0] > 0){
            x_lamb_S = (lambdasInjet_pt[0]/jet_pt[0]);
            rho_lamb_S = lambdasInjet_rho[0];
            d_lamb_S = lambdasInjet_d[0];
            if (x_lamb_S < 0) x_lamb_S = -1;
            if (rho_lamb_S < 0) rho_lamb_S = -999;
            if (d_lamb_S < 0) d_lamb_S = -99;
          }
        }
        if(x_lamb_S > x_lamb_B) x_lamb = x_lamb_S;
        else x_lamb = x_lamb_B;
        if(rho_lamb_S > rho_lamb_B) rho_lamb = rho_lamb_S;
        else rho_lamb = rho_lamb_B;
	if(d_lamb_S > d_lamb_B) d_lamb = d_lamb_S;
	else d_lamb = d_lamb_B;
/*
        if(abs(jet->PID) == 3) std::cout << "pt_lamb_S check : " << lambdasInjet_pt[1] << " , " << jet_pt[1] << " , ratio : " << (lambdasInjet_pt[1]/jet_pt[1]) << std::endl;
	if(abs(jet->PID) == 5) std::cout << "pt_lamb_S check : " << lambdasInjet_pt[0] << " , " << jet_pt[0] << " , ratio : " << (lambdasInjet_pt[0]/jet_pt[0])  << std::endl;
        if(abs(jet->PID) == 3) std::cout << "pt_lamb_B check : " << lambdasInjet_pt[0] << " , " << jet_pt[0] << " , ratio : " << (lambdasInjet_pt[0]/jet_pt[0]) << std::endl;
        if(abs(jet->PID) == 5) std::cout << "pt_lamb_B check : " << lambdasInjet_pt[1] << " , " << jet_pt[1] << " , ratio : " << (lambdasInjet_pt[1]/jet_pt[1])  << std::endl;
        std::cout << "x_lamb check : " << x_lamb << " , " << x_lamb_S << " , " << x_lamb_B << std::endl;
        std::cout << "rho_lamb check : " << rho_lamb << " , " << rho_lamb_S << " , " << rho_lamb_B << std::endl;
*/
      }

      // save highest pT lepton only
      nleptonsInjet.push_back(leptonsInjet.size());
      if (leptonsInjet.size() >0){
        sort(leptonsInjet.begin(), leptonsInjet.end(), [](GenParticle a, GenParticle b){return a.PT > b.PT;});
        leptonsInjet_pt.push_back(leptonsInjet[0].PT);
        leptonsInjet_eta.push_back(leptonsInjet[0].Eta);
        leptonsInjet_phi.push_back(leptonsInjet[0].Phi);
        leptonsInjet_energy.push_back(leptonsInjet[0].E);
        auto R = sqrt(pow(leptonsInjet[0].X,2)+pow(leptonsInjet[0].Y,2)+pow(leptonsInjet[0].Z,2));
        leptonsInjet_R.push_back(R);
      }
      else {

        leptonsInjet_pt.push_back(-99);
        leptonsInjet_eta.push_back(-99);
        leptonsInjet_phi.push_back(-99);
        leptonsInjet_energy.push_back(-99);
        leptonsInjet_R.push_back(-99);
      }

      std::vector<float> diHadron_mass;
      if (hadronsInjet.size() >= 2){
        //std::vector<float> diHadron_energy;
        //std::vector<bool> diHadron_sameMother;
        for (unsigned ih1 = 0; ih1 < hadronsInjet.size(); ++ih1){
          auto hadron1 = hadronsInjet[ih1].P4();
          for (unsigned ih2 = 0; ih2 < hadronsInjet.size(); ++ih2){
            if (ih1 >= ih2) continue;
            auto hadron2 = hadronsInjet[ih2].P4();

            auto diHadron = hadron1+hadron2;
            diHadron_mass.push_back(diHadron.M());
	    if(abs(jet->PID) == 3) histo_dihadron_S->Fill(diHadron.M()*1000);
            if(abs(jet->PID) == 5) histo_dihadron_B->Fill(diHadron.M()*1000);

            //diHadron_energy.push_back(diHadron.E());
            //diHadron_sameMother.push_back(hadronsInjet[ih1].M1 == hadronsInjet[ih2].M1 );
          }
        }
      }
      jet_diHadron_mass.push_back(diHadron_mass);

    }

    jet1_diHadron_mass = jet_diHadron_mass[0];
    jet2_diHadron_mass = jet_diHadron_mass[1];

    // accepted leptons
    int selectedGenLep = 0;
    if (abs(genLeps[0].tlv.Eta()) < 2.4 && genLeps[0].tlv.Pt() > 20.) selectedGenLep++;
    if (abs(genLeps[1].tlv.Eta()) < 2.4 && genLeps[1].tlv.Pt() > 20.) selectedGenLep++;

    if (selectedGenLep == 2){
      gen_step0 = true;
      auto genDilep = genLeps[0].tlv+genLeps[1].tlv;
      if (genDilep.M() > 20. && genLeps[0].charge * genLeps[1].charge < 0 ) gen_step1 = true;
    }

}

//not completed
void Finder(TTree * outtr){
  auto dRCut = 0.5;
  bool matched = false;
  std::vector<const GenParticle*> genJets;
  for (unsigned i = 0; i < particles->GetEntries(); ++ i){
    auto p = (const GenParticle*) particles->At(i);
    if (p->Status > 30) continue;
    if (abs(p->PID) != 6) continue;
    auto top = getLast(particles,p);
    auto jet = (const GenParticle*) particles->At(top->D2);
    genJets.push_back(jet);
  }
  auto njet = jets->GetEntries();
  for(size_t i = 0; i < njet; ++i){
    auto jet = (Jet*) jets->At(i);
    for (auto& p : genJets) {
      float dR = DeltaR(jet->Eta - p->Eta, DeltaPhi(jet->Phi, p->Phi));
      if (dR < dRCut) {
	matched = true;
      }
    }
    if(!matched) continue;
    auto ndau = jet->Constituents.GetEntries();
    for(size_t j=0; j < ndau; ++j){
      auto dau = jet->Constituents.At(j);  
    }
  outtr->Fill();
  }
  jetFinder = true;  
}
std::vector<float> collectHadron(std::vector<GenParticle> hadronsInjet, bool motherCheck){
}
