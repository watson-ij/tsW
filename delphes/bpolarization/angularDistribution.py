import ROOT, math
from array import array
execfile("../loadDelphes.py")



def ana(filedir, filename, jetpid):
    hlist=[]
    hlist.append(ROOT.TH1D("momentum sum","momentum sum",20,0,200))
    hlist.append(ROOT.TH1D("angle","angle",100,0,3.0))
    hlist.append(ROOT.TH1D("d2-d1","d2-d1",10,0,10))
    tfiles = ROOT.TFile(filedir+filename+".root")
    trees=tfiles.Get("Delphes")
    for iev, event in enumerate(trees):
        for ip, p in enumerate(event.Particle):
            if p.Status>30: continue
            if abs(p.PID) != 5122: continue # lambda b
            daughters=[]
            jpsy=0
            lamb=0
            hlist[2].Fill(p.D2-p.D1)
            if (p.D2-p.D1)!=1:continue
            for i in range(p.D1,p.D2+1):
                daughters.append(event.Particle[i])
            if daughters[0].PID!=3122 and daughters[1].PID!=3122:continue
            (jpsy, lamb) = (daughters[0], daughters[1]) if daughters[1].PID==3122 else (daughters[1], daughters[0])
            if jpsy.PID != 443: continue
            p_tlv=ROOT.TLorentzVector()
            jpsy_tlv=ROOT.TLorentzVector()
            lamb_tlv=ROOT.TLorentzVector()
            p_tlv.SetPtEtaPhiE(p.PT,p.Eta,p.Phi,p.E)
            jpsy_tlv.SetPtEtaPhiE(jpsy.PT,jpsy.Eta,jpsy.Phi,jpsy.E)
            lamb_tlv.SetPtEtaPhiE(lamb.PT,lamb.Eta,lamb.Phi,lamb.E)
            bv=p_tlv.BoostVector()
            jpsy_tlv.Boost(-bv)
            lamb_tlv.Boost(-bv)
            hlist[0].Fill((jpsy_tlv+lamb_tlv).Pt())
            #hlist[1].Fill(jpsy_tlv.Angle(lamb_tlv.BoostVector()))
            mom=event.Particle[p.M1]
            if True:
                print 'bottom!'
                b_tlv=ROOT.TLorentzVector()
                b_tlv.SetPtEtaPhiE(mom.PT,mom.Eta,mom.Phi,mom.E)
                n_vector=b_tlv.BoostVector().Cross(p_tlv.BoostVector())
                angle=n_vector.Angle(lamb_tlv.BoostVector())
                print "angle"+str( angle )
                hlist[1].Fill(angle)
    return hlist

             


#filedir = "/mnt/delphes/" #-B /home/scratch/tsW:/mnt
filedir = "/mnt/" #-B /home/iwatson/tsW/:/mnt
h_plz = ana(filedir,"l",5)
#h_plz = ana(filedir,"tsWbW_1k",5)
#h_plz = ana(filedir,"tsWbW_20k",5)
#h_plz = ana(filedir,"tsWbW_100k",5)

rt = ROOT.TFile("test_l.root","RECREATE")
#rt = ROOT.TFile("test_1k.root","RECREATE")
#rt = ROOT.TFile("test_100k.root","RECREATE")
#rt = ROOT.TFile("test_100k.root","RECREATE")
#rt = ROOT.TFile("deltaRdistribution_100k.root","RECREATE")
rt.mkdir("plz")

rt.cd("plz")
for hist in h_plz: hist.Write()
rt.Close()
