#!/usr/bin/env python

import ROOT
from ROOT import TFile, TTree
from glob import glob
import os

if __name__ == "__main__":
    # flist = glob("result/*14_*tsW/nano*root")
    # cmd = "hadd -f9 cmssw_tsw.root "+" ".join(flist)
    # # print cmd
    # os.system(cmd)
    # ###
    # flist = glob("result/*14_*tbW/nano*root")
    # cmd = "hadd -f9 cmssw_tbw.root "+" ".join(flist)
    # # print cmd
    # os.system(cmd)
    flist = glob("result/*222_*ttsWbW/nanoAOD.root")
    # find result -name "nanoAOD.root" -exec ls -l \{\} \; | grep "Feb 23" | grep "15:"
    # flist = []
    cmd = "hadd -f9 cmssw_ttswbw.root "+" ".join(flist)
    # print cmd
    os.system(cmd)
