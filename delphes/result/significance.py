import ROOT, sys, math
from histoHelper import *

bfiledir = './new_wj_out_tbWbW_100k.root'
sfiledir = './new_wj_out_tsWbW_100k.root'
treename = "tsW"

name = "significance"
binning = [100,0,1]

x0 = "kshortsInjet_energy[0]/jet_energy[0]"
x1 = "kshortsInjet_energy[1]/jet_energy[1]"
#x = "max\(kshortsInjet_energy[0]/jet_energy[0],kshortsInjet_energy[1]/jet_energy[1]\)"
h_b = makeTH1(bfiledir, treename, "bW", binning, x0, "%s>%s"%(x0,x1))
h_s = makeTH1(sfiledir, treename, "sW", binning, x0, "%s>%s"%(x0,x1))
h_b.Add(makeTH1(bfiledir, treename, "bW", binning, x1, "%s<%s"%(x0,x1)))
h_s.Add(makeTH1(sfiledir, treename, "sW", binning, x1, "%s<%s"%(x0,x1)))
h_b.Scale(4453997/h_b.Integral())
h_s.Scale(14529/h_s.Integral())
h_b.SetLineColor(1)
h_s.SetLineColor(2)
h_b.SetLineWidth(2)
h_s.SetLineWidth(2)

y=array.array('f',[])
x=array.array('f',[i/10. for i in range(0,100)])
"""
for i in range(0,100):
    s = h_s.Integral(0,i+1)
    b = h_b.Integral(0,i+1)
    significance = s/math.sqrt(s+b)
    #print significance
    y.append(significance)
"""
x=array.array('f',[i/10. for i in range(100,0,-1)])
for i in range(100,0,-1):
    s = h_s.Integral(i,101)
    b = h_b.Integral(i,101)
    significance = s/math.sqrt(s+b)
    #print significance
    y.append(significance)
gr = ROOT.TGraph(len(x), x, y)
gr.SetTitle("Significance (x_Ks)")
gr.SetLineWidth(2)

c = ROOT.TCanvas()
h_b.SetStats(0)
h_b.Draw("hist")
h_s.Draw("samehist")
c.SetLogy()
c.SaveAs("x_Ks.png")

c1 = ROOT.TCanvas()
gr.Draw()
c1.SaveAs("sig.png")
