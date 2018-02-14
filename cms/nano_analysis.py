#!/usr/bin/env python

import ROOT
from ROOT import TFile, TCanvas, TTree, TH1F
from math import hypot, pi

f = TFile("result/20180212_135259_wn4045.sdfarm.kr_tbW/nanoAOD.root")
f = TFile("result/20180212_135216_wn4045.sdfarm.kr_tsW/nanoAOD.root")
f = TFile("result/20180214_132212_cms-t3-wn3022.sdfarm.kr_ttsWbW/nanoAOD.root")
# # f.Events.Print()

def deltaPhi(phi1,phi2):
    ## Catch if being called with two objects
    if type(phi1) != float and type(phi1) != int:
        phi1 = phi1.phi
    if type(phi2) != float and type(phi2) != int:
        phi2 = phi2.phi
    ## Otherwise
    dphi = (phi1-phi2)
    while dphi >  pi: dphi -= 2*pi
    while dphi < -pi: dphi += 2*pi
    return dphi

def deltaR(eta1,phi1,eta2=None,phi2=None):
    ## catch if called with objects
    if eta2 == None:
        return deltaR(eta1.eta,eta1.phi,phi1.eta,phi1.phi)
    ## otherwise
    return hypot(eta1-eta2, deltaPhi(phi1,phi2))

import collections
nS = 0; nSM = 0
nB = 0; nBM = 0
nSB = collections.defaultdict(lambda: 0)
for e in f.Events:
    nS = nB = 0
    q, qb = None, None
    isS = False
    for j in xrange(f.Events.nGenPart):
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == 3:  q = j; nS += 1; isS = True
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == -3: qb = j; nS += 1
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == 5:  q = j; nB += 1
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == -5: qb = j; nB += 1
    if q is None: continue
    nSB[(nS, nB)] += 1
    for j in xrange(f.Events.nKshort):
        if deltaR(f.Events.Kshort_eta[j], f.Events.Kshort_phi[j],
                  f.Events.GenPart_eta[q], f.Events.GenPart_phi[q]) < 0.5:
            if f.Events.Kshort_pt[j] / f.Events.GenPart_eta[q] < 0.2: continue
            if isS: nSM += 1
            else: nBM += 1
#            print "Match q ", f.Events.Kshort_pt[j] / f.Events.GenPart_pt[q], f.Events.Kshort_pt[j]
        if deltaR(f.Events.Kshort_eta[j], f.Events.Kshort_phi[j],
                  f.Events.GenPart_eta[qb], f.Events.GenPart_phi[qb]) < 0.5:
            if f.Events.Kshort_pt[j] / f.Events.GenPart_eta[qb] < 0.2: continue
            if isS: nSM += 1
            else: nBM += 1
#            print "Match qb", f.Events.Kshort_pt[j] / f.Events.GenPart_pt[qb], f.Events.Kshort_pt[j]

#print "S ", nS, "B ", nB
#print "SM", nSM, "BM", nBM

for k, v in nSB.items():
    print k, v
