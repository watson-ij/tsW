import ROOT, sys
from histoHelper import *

'''
python drawer.py nKs nkshortsInjet[0] 8,0,8
python drawer.py nLambs nlambdasInjet[0] 8,0,8
python drawer.py x_Ks kshortsInjet_energy[0]/jet_energy[0] 200,0,1
python drawer.py x_Ks_withCut kshortsInjet_energy[0]/jet_energy[0] 200,0,1
python drawer.py x_Lamb lambdasInjet_energy[0]/jet_energy[0] 200,0,1
python drawer.py displacement_ks kshortsInjet_outR[0]-kshortsInjet_R[0] 200,0,2000
python drawer.py r_lepton leptonsInjet_R[0] 200,0,2000
'''

bfiledir = './wj_out_tsWbW_100k.root'
sfiledir = './wj_out_tsWbW_100k.root'
treename = "tsW"

name = sys.argv[1] 
plotvar = sys.argv[2] 
binning = [int(l) for l in (sys.argv[3]).split(',')]

for i in range(0,100):
    cut = "1==1"
    #cut = "(abs(kshortsInjet_eta[%d])<1.5&&kshortsInjet_outR[%d]<868)||(abs(kshortsInjet_eta[%d])>1.5&&kshortsInjet_outR[%d]<2200)"%(i,i,i,i)
    #cut = "((kshortsInjet_energy[%d]/jet_energy[%d])>0.5)"%(i,i)
    print cut
    if i == 0:
        h_b = makeTH1(bfiledir, treename, "bW", binning, plotvar, "abs(jet_pid[0])==5&&%s"%cut)
        h_s = makeTH1(sfiledir, treename, "sW", binning, plotvar, "abs(jet_pid[0])==3&&%s"%cut)
        continue
    h_b.Add(makeTH1(bfiledir, treename, "bW", binning, plotvar.replace('[0]', '[%d]'%i), "abs(jet_pid[%d])==5&&%s"%(i,cut)))
    h_s.Add(makeTH1(sfiledir, treename, "sW", binning, plotvar.replace('[0]', '[%d]'%i), "abs(jet_pid[%d])==3&&%s"%(i,cut)))
    if "jet_" in plotvar or "Injet" in plotvar: break

print h_b.Integral()
print h_s.Integral()
hlist = [h_b, h_s]
drawTH1(name, hlist, name, "# of jets", False)
#drawTH1(name, hlist, name, "dN/dx", True)
