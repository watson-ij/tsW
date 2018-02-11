// struct for selected leptons
struct Lepton {
   TLorentzVector tlv;
   int charge;
   int pdgid;
};

//read
TClonesArray *gen_jets = 0, *particles = 0;
TClonesArray *electrons = 0, *muons = 0, *jets = 0, *missingET = 0;
TClonesArray *tracks = 0;

//write
float dilepton_mass;
int dilepton_ch, step;

bool gen_step0, gen_step1;

std::vector<float> jet_pt, jet_eta, jet_phi, jet_energy;
std::vector<int> jet_pid;

std::vector<float> kshortsInjet_pt, kshortsInjet_eta, kshortsInjet_phi, kshortsInjet_energy, kshortsInjet_R, kshortsInjet_outR;
std::vector<float> lambdasInjet_pt, lambdasInjet_eta, lambdasInjet_phi, lambdasInjet_energy, lambdasInjet_R, lambdasInjet_outR;
std::vector<float> leptonsInjet_pt, leptonsInjet_eta, leptonsInjet_phi, leptonsInjet_energy, leptonsInjet_R;
std::vector<float> jet1_diHadron_mass, jet2_diHadron_mass;
std::vector<int> nkshortsInjet, nlambdasInjet, nleptonsInjet;


// define functions
struct Lepton toLepton(Muon* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=11*(p->Charge); return l;}
struct Lepton toLepton(Electron* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=13*(p->Charge); return l;}
struct Lepton toLepton(const GenParticle* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=p->PID; return l;}

void defBrancheFucns(TTree* outtr);
const GenParticle* getLast(TClonesArray * particles, const GenParticle* p);
std::vector<const GenParticle*> getMlist(TClonesArray * particles, const GenParticle* p);
void initValues();
void recoParticle();
void genParticle();
std::vector<float> collectHadron(std::vector<GenParticle> hadronsInjet, bool motherCheck);

void defBranchFucns(TTree* outtr){
  #define Branch_(type, name, suffix) outtr->Branch(#name, &name, #name "/" #suffix);
  #define BranchI(name) Branch_(Int_t, name, I)
  #define BranchF(name) Branch_(Float_t, name, F)
  #define BranchO(name) Branch_(Bool_t, name, O)

  #define BranchA_(type, name, size, suffix) outtr->Branch(#name, &name, #name"["#size"]/"#suffix);
  #define BranchAI(name, size) BranchA_(Int_t, name, size, I);
  #define BranchAF(name, size) BranchA_(Float_t, name, size, F);
  #define BranchAO(name, size) BranchA_(Bool_t, name, size, O);
  #define BranchVF(name) outtr->Branch(#name, "vector<float>", &name);
  #define BranchVI(name) outtr->Branch(#name, "vector<int>", &name);
  #define BranchVO(name) outtr->Branch(#name, "vector<bool>", &name);
 
  #define BranchP_(type, br, name, suffix) TBranch *br =  outtr->Branch(#name, &name, #name "/" #suffix);
  #define BranchPI(br,name) BranchP_(Int_t, br,name, I);
  #define BranchPF(br,name) BranchP_(Float_t,br, name, F);
  #define BranchPO(br,name) BranchP_(Bool_t,br, name, O);

  BranchF(dilepton_mass);
  BranchI(dilepton_ch);
  BranchI(step);

  BranchO(gen_step0);
  BranchO(gen_step1);

  BranchVF(jet_pt);
  BranchVF(jet_eta);
  BranchVF(jet_phi);
  BranchVF(jet_energy);
  BranchVI(jet_pid);

  BranchVF(kshortsInjet_pt);
  BranchVF(kshortsInjet_eta);
  BranchVF(kshortsInjet_phi);
  BranchVF(kshortsInjet_energy);
  BranchVF(kshortsInjet_R);
  BranchVF(kshortsInjet_outR);

  BranchVF(lambdasInjet_pt);
  BranchVF(lambdasInjet_eta);
  BranchVF(lambdasInjet_phi);
  BranchVF(lambdasInjet_energy);
  BranchVF(lambdasInjet_R);
  BranchVF(lambdasInjet_outR);

  BranchVF(leptonsInjet_pt);
  BranchVF(leptonsInjet_eta);
  BranchVF(leptonsInjet_phi);
  BranchVF(leptonsInjet_energy);
  BranchVF(leptonsInjet_R);
  BranchVF(leptonsInjet_outR);

  BranchVI(nkshortsInjet);
  BranchVI(nlambdasInjet);
  BranchVI(nleptonsInjet);

  BranchVF(jet1_diHadron_mass);
  BranchVF(jet2_diHadron_mass);

  return;
}

const GenParticle* getLast(TClonesArray * particles, const GenParticle* p){
  auto mom = p;
  while(true){
    auto dau = (const GenParticle*)particles->At(mom->D1);
    if( abs(p->PID) != abs(dau->PID) ) break;
    mom = dau;
  }
  return mom;
}

std::vector<const GenParticle*> getMlist(TClonesArray * particles, const GenParticle* p){
  std::vector<const GenParticle*> mlst;
  auto idx = p->M1;
  if (idx == -1) return mlst;
  while(true){
    auto m = (const GenParticle*)particles->At(idx);
    mlst.push_back(m);
    idx = m->M1;
    if ( idx == -1) break;
  }
  return mlst;
}
