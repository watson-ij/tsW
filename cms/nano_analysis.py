#!/usr/bin/env python

import ROOT
from ROOT import TFile, TCanvas, TTree, TH1F
from math import hypot, pi

f = TFile("nanoAOD.root")
# f.Events.Print()

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

for e in f.Events:
    q, qb = None, None
    for j in xrange(f.Events.nGenPart):
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == 3:  q = j
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == -3: qb = j
    if q is None: continue
    for j in xrange(f.Events.nKshort):
        if deltaR(f.Events.Kshort_eta[j], f.Events.Kshort_phi[j],
                  f.Events.GenPart_eta[q], f.Events.GenPart_phi[q]) < 0.5:
            print "Match q ", f.Events.Kshort_pt[j] / f.Events.GenPart_pt[q], f.Events.Kshort_pt[j]
        if deltaR(f.Events.Kshort_eta[j], f.Events.Kshort_phi[j],
                  f.Events.GenPart_eta[qb], f.Events.GenPart_phi[qb]) < 0.5:
            print "Match qb", f.Events.Kshort_pt[j] / f.Events.GenPart_pt[qb], f.Events.Kshort_pt[j]
