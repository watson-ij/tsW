import ROOT, math
from array import array

rt = ROOT.TFile("deltaRdistribution_l.root")
rtdir = rt.Get('plz')

namelist=[]
namelist.append(['deltaR' ,'deltaR(b->lambdab)','deltaR(smallest deltaR)'])
namelist.append(['deltaR vs pt ratio of lambdab','deltaR vs pt ratio of lambdab (b->lambdab)','deltaR vs pt ratio of lambdab (smallest deltaR)'])
namelist.append(['pt ratio','pt ratio(b->lambdab)','pt ratio(smallest deltaR)'])

def selectHList(rtdir, namelist):
    selectedHList=[]
    keyList=rtdir.GetListOfKeys()
    for key in keyList:
        tobj=key.ReadObj()
        if type(namelist)==str:
            if tobj.GetName().startswith(namelist):
                selectedHList.append(tobj)
            else : continue
        elif type(namelist)==list:
            if tobj.GetName() in namelist:
                selectedHList.append(tobj)
            else : continue

    if len(selectedHList)==0: print "There isn't any histogram whose name starts with " + name
    return selectedHList

def draw(rtdir, namelist):
    for name in namelist:
        hlist=selectHList(rtdir,name)
        canv = ROOT.TCanvas(name[0],name[0],600,800)
        canv.cd()
        hs = ROOT.THStack(name[0],name[0])
        color=[ROOT.kBlack,ROOT.kRed,ROOT.kViolet]
        hlist.sort(key=lambda x:-x.GetMaximum())
        for i, h in enumerate(hlist):
            h.SetMarkerStyle(21)
            h.SetFillColor(color[i])
            h.SetMarkerColor(ROOT.kBlue)
            hs.Add(h)
        hs.Draw('nostack')
        canv.BuildLegend()
        canv.Write()

tf = ROOT.TFile('result.root','RECREATE')
draw(rtdir,namelist)
tf.Close()
