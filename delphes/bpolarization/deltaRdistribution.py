import ROOT, math
from array import array
execfile("../loadDelphes.py")

def getLast(event,p): #returns the last particle whose daughter is not itself.
    mom = p
    i=0
    while True:
        dau = event.Particle[mom.D1]
        i+=1
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
def initHlist():
    hlist=[]
    hlist.append(ROOT.TH1D("deltaR","deltaR",100,0,0.5)) #0
    hlist.append(ROOT.TH1D("deltaR(b->lambdab)","deltaR(b->lambdab)",100,0,0.5)) #1
    hlist.append(ROOT.TH1D("deltaR(smallest deltaR)","deltaR(smallest deltaR)",100,0,0.5)) #2
    hlist.append(ROOT.TH2D("deltaR vs pt ratio of lambdab","deltaR vs pt ratio of lambdab",100,0,3,100,0,0.5)) #3
    hlist.append(ROOT.TH2D("deltaR vs pt ratio of lambdab (b->lambdab)","deltaR vs pt ratio of lambdab (b->lambdab)",100,0,3,100,0,0.5)) #4
    hlist.append(ROOT.TH2D("deltaR vs pt ratio of lambdab (smallest deltaR)","deltaR vs pt ratio of lambdab (smallest deltaR)",100,0,3,100,0,0.5)) #5
    hlist.append(ROOT.TH1D("pt ratio","pt ratio",100,0,4)) #6
    hlist.append(ROOT.TH1D("pt ratio(b->lambdab)","pt ratio(b->lambdab)",100,0,4)) #7
    hlist.append(ROOT.TH1D("pt ratio(smallest deltaR)","pt ratio(smallest deltaR)",100,0,4)) #8
    
    hlist.append(ROOT.TH2D("deltaR vs pt ratio of lambdab (smallest deltaR with cut)","deltaR vs pt ratio of lambdab (smallest deltaRcut )",100,0,4,100,0,0.5)) #9
    
    hlist[3].GetXaxis().SetTitle("pt(lambda_b)/pt(b)")
    hlist[3].GetYaxis().SetTitle("deltaR")
    hlist[4].GetXaxis().SetTitle("pt(lambda_b)/pt(b)")
    hlist[4].GetYaxis().SetTitle("deltaR")
    hlist[5].GetXaxis().SetTitle("pt(lambda_b)/pt(b)")
    hlist[5].GetYaxis().SetTitle("deltaR")
    hlist[6].GetXaxis().SetTitle("pt(lambda_b)/pt(b)")
    hlist[7].GetXaxis().SetTitle("pt(lambda_b)/pt(b)")
    hlist[8].GetXaxis().SetTitle("pt(lambda_b)/pt(b)")
    hlist[9].GetXaxis().SetTitle("pt(lambda_b)/pt(b)")
    hlist[9].GetYaxis().SetTitle("deltaR")
    return hlist

def ana(filedir, filename, jetpid):
    hlist=initHlist()

    tfiles = ROOT.TFile(filedir+filename+".root")
    trees=tfiles.Get("Delphes")
    daughters=[]
    score=[0,0]
    for iev, event in enumerate(trees):
        for ip, p in enumerate(event.Particle):
            # from mother
            if p.Status>30: continue
            if abs(p.PID)==5:
                getLast(event,p)
            if abs(p.PID) != 5122: continue # lambda b
            mlist=getMlist(event,p)
            jet=0
            for m in mlist:
                if m.PID == jetpid and m.Status <= 30:
                    jet=m
                    break

            if jet!=0:
                deltaR=getdr(jet,p)
                hlist[0].Fill(deltaR)
                hlist[3].Fill(p.PT/jet.PT,deltaR)
                hlist[6].Fill(p.PT/jet.PT)
                if mlist[0].PID == jetpid:
                    hlist[1].Fill(deltaR)
                    hlist[4].Fill(p.PT/jet.PT,deltaR)
                    hlist[7].Fill(p.PT/jet.PT)
            
            # from a list of b jet sorted by deltaR
            deltaRlist=[]
            for iq, q in enumerate(event.Particle):
                if q.PID != jetpid: continue
                if q.Status > 30: continue
                if q.PT==0: continue
                ptratio=p.PT/q.PT
                if ptratio<0.6 and ptratio>1.3: continue
                temp_deltaR=getdr(q,p)
                deltaRlist.append([q,temp_deltaR,ptratio])
            
            deltaRlist.sort(key=lambda x:x[0]) # sorting by deltaR
            #deltaR_jet=getLast(event,deltaRlist[1][1]) # get the last b of the smallest deltaR b
            if len(deltaRlist)>0:
                deltaR_jet=deltaRlist[0][0]
                deltaR=deltaRlist[0][1]
                ptratio=deltaRlist[0][2]
                hlist[2].Fill(deltaR)
                hlist[5].Fill(ptratio,deltaR)
                hlist[8].Fill(ptratio)
                hlist[9].Fill(ptratio,deltaR)
                if deltaR_jet in mlist:
                    score[0]=score[0]+1
                else: score[1]=score[1]+1
    print '# of b jet selected by minimum deltaR accurately : ' + str(score[0]) +" , " + str(int(score[0]*100/(score[0]+score[1]))) + " %"
    print '# of b jet selected by minimum deltaR wrongly : ' + str(score[1]) +" , " + str(int(score[1]*100/(score[0]+score[1]))) + " %"
    
    
    return hlist


#filedir = "/mnt/delphes/" #-B /home/scratch/tsW:/mnt
filedir = "/mnt/" #-B /home/iwatson/tsW/:/mnt
h_plz = ana(filedir,"l",5)
#h_plz = ana(filedir,"tsWbW_1k",5)
#h_plz = ana(filedir,"tsWbW_20k",5)
#h_plz = ana(filedir,"tsWbW_100k",5)

#rt = ROOT.TFile("test_l.root","RECREATE")
#rt = ROOT.TFile("test_1k.root","RECREATE")
#rt = ROOT.TFile("test_100k.root","RECREATE")
#rt = ROOT.TFile("deltaRdistribution_100k.root","RECREATE")
rt = ROOT.TFile("deltaRdistribution_l.root","RECREATE")
rt.mkdir("plz")
rt.cd("plz")

for hist in h_plz: hist.Write()
rt.Close()
