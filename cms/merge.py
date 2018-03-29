#!/usr/bin/env python

import ROOT
from ROOT import TFile, TTree
from glob import glob
import os

if __name__ == "__main__":
    if False:
        print "Running V0Sim" 
        flist = glob("result/201803*ttsWbW/step3.root")
        for f in flist:
            d = os.path.dirname(f)
            os.system("cd %s && V0Sim step3.root" % d)
            os.system("hadd -f9 %s %s %s" % (d+"/nv.root", d+"/nanoAOD.root", d+"/v0sim.root"))
        print "Merging V0Sim/NanoAOD"
        flist = glob("result/201803*ttsWbW/nv.root")
        os.system("hadd -f9 cmssw_ttswbw_v0sim.root "+" ".join(flist))

        print "Continuing with others"
        # flist = glob("result/*14_*tsW/nano*root")
        # cmd = "hadd -f9 cmssw_tsw.root "+" ".join(flist)
        # # print cmd
        # os.system(cmd)
        # ###
        # flist = glob("result/*14_*tbW/nano*root")
        # cmd = "hadd -f9 cmssw_tbw.root "+" ".join(flist)
        # # print cmd
        # os.system(cmd)
        flist = glob("result/20180222_*ttsWbW/nanoAOD.root")
        flist.extend(glob("result/20180306_*ttsWbW/nanoAOD.root"))
        print (len(flist))
        # find result -name "nanoAOD.root" -exec ls -l \{\} \; | grep "Feb 23" | grep "15:"
        # flist = []
        cmd = "hadd -f9 cmssw_ttswbw.root "+" ".join(flist)
        # print cmd
        os.system(cmd)

    # Lambda b
    flist = glob("result/*lam/nanoAOD.root")
    cmd = "hadd -f9 cmssw_tt_lambdab.root "+" ".join(flist)
    # print cmd
    os.system(cmd)
