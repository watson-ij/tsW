// struct for selected leptons
struct lepton {
  TLorentzVector tlv;
  int charge;
  int pdgid;
};

// define functions
void defBranchFucns();
const GenParticle* getLast(TClonesArray * particles, size_t iev, const GenParticle* p);

void defBranchFucns(){
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
  #define BranchLV(name) TLorentzVector name; outtr->Branch(#name, "TLorentzVector", &name);
  
  #define BranchP_(type, br, name, suffix) type name = 0; TBranch *br =  outtr->Branch(#name, &name, #name "/" #suffix);
  #define BranchPI(br,name) BranchP_(Int_t, br,name, I);
  #define BranchPF(br,name) BranchP_(Float_t,br, name, F);
  #define BranchPO(br,name) BranchP_(Bool_t,br, name, O);
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

