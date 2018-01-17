import ROOT, math
from array import array
execfile("../loadDelphes.py")

def ana(filedir, filename, jetpid):
    hlist=[]
    hlist.append(ROOT.TH1D("number of events","number of events",1,0,1))
    hlist.append(ROOT.TH1D("# of s mesons","# of s mesons",20,0,20))

    tfiles = ROOT.TFile(filedir+filename+".root")
    trees = tfiles.Get("Delphes")

    pids = []
    for iev, event in enumerate(trees):
        hlist[0].Fill(0.5)
        for p in event.Particle:
            strpdg = str(abs(p.PID))
            if abs(p.PID) < 100 : continue
            if not strpdg.startswith("3") : continue
            a = ROOT.TDatabasePDG()
            b = a.GetParticle(abs(p.PID))
            if b == None :
                print p.PID, "is None (why?)"
                continue
            pids.append(b.GetName())

        if iev%100 == 0: print iev

    ranks = {}
    for pid in set(pids):
        ranks[pid] = pids.count(pid)
    sortedpids = sorted(ranks.items(), key=lambda x:x[1], reverse=True)

    for i, pid in enumerate(sortedpids):
        if i==20: break
        hlist[1].GetXaxis().SetBinLabel(i+1,pid[0])
        hlist[1].SetBinContent(i+1,pid[1])

    tfiles.Close()
    return hlist

filedir = "./1000events/"
h_tsW = ana(filedir, "tsW_dilep",3)
h_tbW = ana(filedir, "tbW_dilep",5)

rt = ROOT.TFile("res.root", "RECREATE")
rt.mkdir("tsW")
rt.mkdir("tbW")

rt.cd("tsW")
for hist in h_tsW: hist.Write()

rt.cd("tbW")
for hist in h_tbW: hist.Write()

rt.Close()
