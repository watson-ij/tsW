//define Branch
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

#define BranchP_(type, br, name, suffix) type name = 0; TBranch *br =  outtr->Branch(#name, &name, #name "/" #suffix);
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

//declare functions
struct Lepton toLepton(Muon* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=11*(p->Charge); return l;}
struct Lepton toLepton(Electron* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=13*(p->Charge); return l;}
struct Lepton toLepton(const GenParticle* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=p->PID; return l;}

const GenParticle* getLast(TClonesArray * particles, const GenParticle* p);
std::vector<const GenParticle*> getMlist(TClonesArray * particles, const GenParticle* p);
std::vector<float> collectHadron(std::vector<GenParticle> hadronsInjet, bool motherCheck);
std::vector<Double_t> cross3D(std::vector<Double_t> & a, std::vector<Double_t> & b);
Double_t DeltaPhi(Double_t phi1, Double_t phi2);
Double_t DeltaR(Double_t deta, Double_t dphi);

//define functions
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

