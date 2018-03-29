import ROOT, math
from array import array
execfile("../loadDelphes.py")

def getLast(event,p): #returns the last particle whose daughter is not itself.
    mom = p
    while True:
        dau = event.Particle[mom.D1]
        if abs(p.PID) != abs(dau.PID): break
        mom = dau
    return mom

def getdr(p1, p2): #returns value of Delta R between two particles, p1 and p2.
    p1_tlv = ROOT.TLorentzVector()
    p1_tlv.SetPtEtaPhiE(p1.PT,p1.Eta,p1.Phi,p1.E)
    p2_tlv = ROOT.TLorentzVector()
    p2_tlv.SetPtEtaPhiE(p2.PT,p2.Eta,p2.Phi,p2.E)
    return p1_tlv.DeltaR(p2_tlv)

def getMlist(event, p): #returns a list of mother particles of particle p.
    mlst = []
    idx = p.M1
    while True:
        m = event.Particle[idx]
        mlst.append(m)
        #mlist.append(abs(m.PID)) #returns a list of PID of mother particles of particle p.
        idx = m.M1
        if idx == -1: break
    return mlst

#def getDecayChannels(lambdab,depth):
    #if(lambdab.D1.isFinal()):


def ana(filedir, filename, jetpid):
    hlist=[]
    hlist.append(ROOT.TH1D("number of events","number of events",1,0,1)) #0
    hlist.append(ROOT.TH1D("number of lambda b","number of lambda b",1,0,1)) #1
    hlist.append(ROOT.TH1D("number of b jet","number of b jet",1,0,1)) #2
    hlist.append(ROOT.TH1D("b jet PT","b jet PT",100,0,1)) #3
    hlist.append(ROOT.TH1D("lambda b PT","lambda b PT",100,0,1)) #4 
    hlist.append(ROOT.TH1D("lambda b PT / b PT","lambda b PT / b PT",100,0,1)) #5
    hlist.append(ROOT.TH1D("lambda b (5122)PT / b PT","lambda b (5122)PT / b PT",100,0,1)) #5
    hlist.append(ROOT.TH1D("lambda b (-5122)PT / b PT","lambda b (-5122)PT / b PT",100,0,1)) #5
    hlist.append(ROOT.TH1D("deltaR","deltaR",100,0,1)) #8
    hlist.append(ROOT.TH1D("deltaR2","deltaR2",100,0,1)) #9
    #hlist.append(ROOT.TH2D("tripion x vs mass","tripion x vs mass",1000,0,5000,100,0,1))
    #hlist[].GetYaxis().SetTitle("tripion x")

    tfiles = ROOT.TFile(filedir+filename+".root")
    trees=tfiles.Get("Delphes")
    daughters=[]
    score=[0,0]
    for iev, event in enumerate(trees):
        for ip, p in enumerate(event.Particle):
            if p.Status>30: continue
            if abs(p.PID) != 5122: continue # lambda b
            mlist=getMlist(event,p)
            jet=-1
            for m in mlist:
                if m.PID == jetpid:
                    jet=m
                    break
            if jet!=-1:
                deltaR=getdr(jet,p)
                hlist[8].Fill(deltaR)
            deltaR=1
            deltaR_jet=0
            for iq, q in enumerate(event.Particle):
                if q.PID != jetpid: continue
                temp_deltaR=getdr(q,p)
                if temp_deltaR > 0.5: continue
                if temp_deltaR<deltaR : 
                    deltaR=temp_deltaR
                    deltaR_jet=q
            hlist[9].Fill(deltaR)
            if jet==deltaR_jet:
                score[0]=score[0]+1
            else: score[1]=score[1]+1
    print '# of b jet selected by minimum deltaR accurately : ' score[0]
    print '# of b jet selected by minimum deltaR wrongly : ' score[1]
    return hlist

'''
    for iev, event in enumerate(trees):
        hlist[0].Fill(0.5) # count number of events
        for ip, p in enumerate(event.Particle):
            if p.Status>30: continue
            if abs(p.PID) == 5122: # lambda b
                temp=[event.Particle[p.D1].PID,event.Particle[p.D2].PID]
                temp.sort()
                daughters.append(temp)
                hlist[1].Fill(0.5) # count number of lambda b
            if abs(p.PID) != 6 : continue #not top
            lastTop = getLast(event,p)
            jet = event.Particle[lastTop.D2]

            if abs(jet.PID) != jetpid : continue #not b jet
            hlist[2].Fill(0.5)
            hlist[3].Fill(abs(jet.PT))

            for iq, q in enumerate(event.Particle):
                if q.PID not in [-5122,5122]: continue
                deltaR = getdr(jet,q)
                if deltaR > 0.5: continue #out of b jet
                if abs(q.PID)==5122:
                    hlist[4].Fill(q.PT)
                    hlist[5].Fill(q.PT/jet.PT)
                    if q.PID==5122:
                        hlist[6].Fill(q.PT/jet.PT)
                    elif q.PID==-5122:
                        hlist[7].Fill(q.PT/jet.PT)
    daughters.sort()
    print daughters 
    return hlist
'''

#filedir = "/mnt/delphes/" #-B /home/scratch/tsW:/mnt
filedir = "/mnt/" #-B /home/iwatson/tsW/:/mnt
h_plz = ana(filedir,"l",5)
#h_plz = ana(filedir,"tsWbW_1k",5)
#h_plz = ana(filedir,"tsWbW_20k",5)
#h_plz = ana(filedir,"tsWbW_100k",5)

rt = ROOT.TFile("test_l.root","RECREATE")
#rt = ROOT.TFile("test_1k.root","RECREATE")
#rt = ROOT.TFile("test_100k.root","RECREATE")
rt.mkdir("plz")

rt.cd("plz")
for hist in h_plz: hist.Write()
rt.Close()
