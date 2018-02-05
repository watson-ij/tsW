import ROOT

f = ROOT.TFile("result/wj_out_tsWbW_100k.root","READ")
t = f.Get("tsW")

h = ROOT.TH1D("x","x",100,0,1)
t.Draw("Ks_energy[0]/jet_energy[0] >> x", "abs(jet_pid[0])==3")
h.Scale(h.Integral()/14529.)

cnv = ROOT.TCanvas()
cnv.SetLogy()
h.SetMinimum(1)
h.Draw("hist")
cnv.Print("test.png")

