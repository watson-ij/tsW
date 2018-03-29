//define
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

//struct for selected leptons
struct Lepton {
  TLorentzVector tlv;
  int charge;
  int pdgid;
};

//read data
TClonesArray *gen_jets = 0, *particles = 0;
TClonesArray *electrons = 0, *muons = 0, *jets = 0, *missingET = 0;
TClonesArray *tracks = 0;

//declare variable for branch
//  recoParticle()

//  genParticle()
float dilepton_mass;

float x_KS, x_KS_S, x_KS_B, rho_KS, rho_KS_S,rho_KS_B, d_KS, d_KS_S, d_KS_B;
float x_lamb, x_lamb_S, x_lamb_B, rho_lamb, rho_lamb_S, rho_lamb_B, d_lamb, d_lamb_S, d_lamb_B;
float significance_S, significance_B;

int dilepton_ch, channel;

bool gen_step0, gen_step1;
std::vector<float> jet_pt, jet_eta, jet_phi, jet_energy;
std::vector<int> jet_pid;

std::vector<float> kshortsInjet_pt, kshortsInjet_eta, kshortsInjet_phi, kshortsInjet_energy, kshortsInjet_R, kshortsInjet_outR, kshortsInjet_rho, kshortsInjet_d;
std::vector<float> lambdasInjet_pt, lambdasInjet_eta, lambdasInjet_phi, lambdasInjet_energy, lambdasInjet_R, lambdasInjet_outR, lambdasInjet_rho, lambdasInjet_d;
std::vector<float> leptonsInjet_pt, leptonsInjet_eta, leptonsInjet_phi, leptonsInjet_energy, leptonsInjet_R, leptonsInjet_outR;

std::vector<int> nkshortsInjet, nlambdasInjet, nleptonsInjet;

std::vector<float> jet1_diHadron_mass, jet2_diHadron_mass;

//declare functions
struct Lepton toLepton(Muon* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=11*(p->Charge); return l;}
struct Lepton toLepton(Electron* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=13*(p->Charge); return l;}
struct Lepton toLepton(const GenParticle* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=p->PID; return l;}

void defBranchGen(TTree* outtr);

const GenParticle* getLast(TClonesArray * particles, const GenParticle* p);
std::vector<const GenParticle*> getMlist(TClonesArray * particles, const GenParticle* p);
std::vector<float> collectHadron(std::vector<GenParticle> hadronsInjet, bool motherCheck);
std::vector<Double_t> cross3D(std::vector<Double_t> & a, std::vector<Double_t> & b);
Double_t DeltaPhi(Double_t phi1, Double_t phi2);
Double_t DeltaR(Double_t deta, Double_t dphi);

void initValues();
void recoParticle(TH1F*);
void genParticle(TH1F*, TH1F*);

//define functions
void defBranchGen(TTree* outtr){

  BranchF(dilepton_mass);

  BranchF(x_KS); BranchF(x_KS_S); BranchF(x_KS_B); BranchF(rho_KS); BranchF(rho_KS_S); BranchF(rho_KS_B); BranchF(d_KS); BranchF(d_KS_S); BranchF(d_KS_B);
  BranchF(x_lamb); BranchF(x_lamb_S); BranchF(x_lamb_B); BranchF(rho_lamb); BranchF(rho_lamb_S); BranchF(rho_lamb_B);  BranchF(d_lamb); BranchF(d_lamb_S); BranchF(d_lamb_B);

  BranchF(significance_S); BranchF(significance_B);

  BranchI(dilepton_ch);
  BranchI(channel);

  BranchO(gen_step0);
  BranchO(gen_step1);

  BranchVF(jet_pt); BranchVF(jet_eta); BranchVF(jet_phi); BranchVF(jet_energy); BranchVI(jet_pid);
  BranchVF(kshortsInjet_pt); BranchVF(kshortsInjet_eta); BranchVF(kshortsInjet_phi); BranchVF(kshortsInjet_energy); BranchVF(kshortsInjet_R); BranchVF(kshortsInjet_outR); BranchVF(kshortsInjet_rho); BranchVF(kshortsInjet_d); 
  BranchVF(lambdasInjet_pt); BranchVF(lambdasInjet_eta); BranchVF(lambdasInjet_phi); BranchVF(lambdasInjet_energy); BranchVF(lambdasInjet_R); BranchVF(lambdasInjet_outR); BranchVF(lambdasInjet_rho); BranchVF(lambdasInjet_d);
  BranchVF(leptonsInjet_pt); BranchVF(leptonsInjet_eta); BranchVF(leptonsInjet_phi); BranchVF(leptonsInjet_energy); BranchVF(leptonsInjet_R); BranchVF(leptonsInjet_outR);
  BranchVI(nkshortsInjet); BranchVI(nlambdasInjet); BranchVI(nleptonsInjet);

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


std::vector<Double_t> cross3D(std::vector<Double_t> & a, std::vector<Double_t> & b){
  std::vector<Double_t> c = { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
  return c;
}

Double_t DeltaPhi(Double_t phi1, Double_t phi2) {
  static const Double_t kPI = TMath::Pi();
  static const Double_t kTWOPI = 2*TMath::Pi();
  Double_t x = phi1 - phi2;
  if(TMath::IsNaN(x)){
    std::cerr << "DeltaPhi function called with NaN" << std::endl;
    return x;
  }
  while (x >= kPI) x -= kTWOPI;
  while (x < -kPI) x += kTWOPI;
  return x;
}

Double_t DeltaR(Double_t deta, Double_t dphi) {
  return TMath::Sqrt(deta*deta + dphi*dphi);
}

