import ROOT, math
from array import array
execfile("../loadDelphes.py")

def getLast(event,p):
    mom = p
    while True:
        dau = event.Particle[mom.D1]
        if abs(p.PID) != abs(dau.PID): break
        mom = dau 
    return mom

def getdr(p1, p2):
    p1_tlv = ROOT.TLorentzVector()
    p1_tlv.SetPtEtaPhiE(p1.PT,p1.Eta,p1.Phi,p1.E)
    p2_tlv = ROOT.TLorentzVector()
    p2_tlv.SetPtEtaPhiE(p2.PT,p2.Eta,p2.Phi,p2.E)
    return p1_tlv.DeltaR(p2_tlv)
    
def getMlist(event, p):
    mlst = []
    idx = p.M1
    while True:
        m = event.Particle[idx]
        mlst.append(m)
        #mlst.append(abs(m.PID))
        idx = m.M1
        if idx == -1: break
    return mlst

def ana(filedir, filename, jetpid):
    hlist=[]
    hlist.append(ROOT.TH1D("number of events","number of events",1,0,1)) #0
    hlist.append(ROOT.TH1D("Energy ratio Ks","Energy ratio Ks",100,0,1))
    hlist.append(ROOT.TH1D("Energy ratio lamb","Energy ratio lamb",100,0,1))
    hlist.append(ROOT.TH1D("displaced length (rho)","displaced length (rho)",300,0,30))
    hlist.append(ROOT.TH1D("displaced length (r)","displaced length (r)",300,0,30))
    hlist.append(ROOT.TH1D("jet id","jet id",7,0,7))
    hlist.append(ROOT.TH1D("lepton isolation","lepton isolation",200,0,2))
    hlist.append(ROOT.TH1D("lepton in the jet","lepton in the jet",2,0,2))
    hlist.append(ROOT.TH1D("pT ratio lepton","pT ratio lepton",100,0,1))
    hlist.append(ROOT.TH1D("Energy ratio lepton","Energy ratio lepton",100,0,1))
    hlist.append(ROOT.TH1D("dipion mass","dipion mass",1000,0,3000)) #10
    hlist.append(ROOT.TH1D("dipion delta r","dipion delta r",100,0,1))
    hlist.append(ROOT.TH1D("dipion x","dipion x",100,0,1))
    hlist.append(ROOT.TH2D("dipion x vs mass","dipion x vs mass",1000,0,3000,100,0,1))
    hlist.append(ROOT.TH2D("tripion x vs mass","tripion x vs mass",1000,0,5000,100,0,1))
    hlist.append(ROOT.TH2D("dipion x vs mass (True)","dipion x vs mass (True)",1000,0,3000,100,0,1))
    hlist.append(ROOT.TH2D("triipion x vs mass (True)","triipion x vs mass (True)",1000,0,5000,100,0,1))
    hlist.append(ROOT.TH1D("tripion mass","tripion mass",1000,0,5000))
    hlist.append(ROOT.TH1D("tripion x","tripion x",100,0,1)) #18

    hlist[13].GetXaxis().SetTitle("dipion mass [MeV]")
    hlist[13].GetYaxis().SetTitle("dipion x")
    hlist[14].GetXaxis().SetTitle("dipion mass [MeV]")
    hlist[14].GetYaxis().SetTitle("dipion x")

    tfiles = ROOT.TFile(filedir+filename+".root")
    trees = tfiles.Get("Delphes")

    for iev, event in enumerate(trees):
        hlist[0].Fill(0.5)
        for p in event.Particle:
            if p.Status > 30: continue
            if abs(p.PID) != 6 : continue
            lastTop = getLast(event,p)
            jet = event.Particle[lastTop.D2]

            if abs(jet.PID) != jetpid : continue
            hlist[5].Fill(abs(jet.PID))

            pions = []
            for q in event.Particle:
                if q.PID not in [211, -211, 321, -321, 2212, -2212, 310, 3122, 11, 13, -11, -13]: continue
                deltaR = getdr(jet, q)
                if deltaR > 0.5: continue

                if abs(q.PID) == 211 or abs(q.PID) == 321 or abs(q.PID) == 2212:
                    pions.append(q)

                elif abs(q.PID) == 310 : #Ks
                    hlist[1].Fill(q.PT/jet.PT)

                elif abs(q.PID) == 3122 : #lambda
                    hlist[2].Fill(q.PT/jet.PT)

                elif abs(q.PID) == 11 or abs(q.PID) == 13 : #lepton
                    isJetLepton = jet in getMlist(event,q)
                    hlist[6].Fill(0)# gen particle has isolation?
                    hlist[7].Fill(isJetLepton)

                    if isJetLepton:
                        hlist[3].Fill(math.sqrt(q.X**2+q.Y**2))
                        hlist[4].Fill(math.sqrt(q.X**2+q.Y**2+q.Z**2))
                        hlist[8].Fill(q.PT/jet.PT)
                        hlist[9].Fill(q.E/jet.E)

                else: continue

            if len(pions) <= 1: continue
            for ip1, pion1 in enumerate(pions):
                for ip2, pion2 in enumerate(pions):
                    if ip1 >= ip2: continue
                    if pion1.PID * pion2.PID > 0: continue
                    pion1_tlv = ROOT.TLorentzVector()
                    pion1_tlv.SetPtEtaPhiE(pion1.PT,pion1.Eta,pion1.Phi,pion1.E)
                    pion2_tlv = ROOT.TLorentzVector()
                    pion2_tlv.SetPtEtaPhiE(pion2.PT,pion2.Eta,pion2.Phi,pion2.E)
                    dipion_tlv = pion1_tlv+pion2_tlv
                    hlist[10].Fill(dipion_tlv.M()*1000)
                    hlist[11].Fill(pion1_tlv.DeltaR(pion2_tlv))
                    hlist[12].Fill(dipion_tlv.E()/jet.E)
                    hlist[13].Fill(dipion_tlv.M()*1000, dipion_tlv.E()/jet.E)
                    if pion1.M1 == pion2.M1:
                        hlist[15].Fill(dipion_tlv.M()*1000, dipion_tlv.E()/jet.E)
                    for ip3, pion3 in enumerate(pions):
                        if ip2 >= ip3: continue
                        pion3_tlv = ROOT.TLorentzVector()
                        pion3_tlv.SetPtEtaPhiE(pion3.PT,pion3.Eta,pion3.Phi,pion3.E)
                        tripion_tlv = dipion_tlv+pion3_tlv
                        hlist[17].Fill(tripion_tlv.M()*1000)
                        hlist[18].Fill(tripion_tlv.E()/jet.E)
                        hlist[14].Fill(tripion_tlv.M()*1000, tripion_tlv.E()/jet.E)
                        if pion1.M1 == pion2.M1 == pion3.M1:
                            hlist[16].Fill(tripion_tlv.M()*1000, tripion_tlv.E()/jet.E)

        if iev%1000 == 0: print iev

    tfiles.Close()
    return hlist

filedir = "/scratch/tsW/"
h_tsW = ana(filedir, "tsWbW_10k",3)
h_tbW = ana(filedir, "tsWbW_10k",5)

rt = ROOT.TFile("test.root", "RECREATE")
rt.mkdir("tsW")
rt.mkdir("tbW")

rt.cd("tsW")
for hist in h_tsW: hist.Write()
rt.cd("tbW")
for hist in h_tbW: hist.Write()

rt.Close()
