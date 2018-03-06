#!/usr/bin/env python2

import ROOT
from ROOT import TFile, TCanvas, TTree, TH1F
from math import hypot, pi, sqrt
ROOT.gROOT.SetBatch(True)

ext = "_dj"
f = TFile("/cms/scratch/iwatson/tsW/cms/cmssw_ttswbw.root")

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

def genMatch(producer, idx, jIdx):
    l_dr=[]
    for ig in xrange(e.nV0GenPart):
        #if e.V0GenPart_pt[ig] / e.Jet_pt[jIdx] < 0.15: continue
        #if abs(e.V0GenPart_mass[ig]-0.49) > 0.03: continue
        if e.V0GenPart_pdgId[ig] !=310: continue
        dr_gen = deltaR(e.V0GenPart_eta[ig], e.V0GenPart_phi[ig],
                        e.Jet_eta[jIdx], e.Jet_phi[jIdx])
        if dr_gen > 0.3 : continue
        if producer == "V": dr_withgen = deltaR(e.V0GenPart_eta[ig], e.V0GenPart_phi[ig],
                                               e.Kshort_eta[idx], e.Kshort_phi[idx])
        elif producer == "M": dr_withgen = deltaR(e.V0GenPart_eta[ig], e.V0GenPart_phi[ig],
                                                 e.meson_eta[idx], e.meson_phi[idx])
        l_dr.append(dr_withgen)
    if len(l_dr) == 0: return 100
    else: return min(l_dr)

hV = TH1F("V", "V0Meson collection;Mass [GeV];", 50, 0.43, 0.57); hV.SetDirectory(0)
hM = TH1F("M", "Meson collection;Mass [GeV];", 50, 0.43, 0.57); hM.SetDirectory(0)

h_genPt = TH1F("genpt", "gen Kshort pT; [GeV];", 50, 0, 50); h_genPt.SetDirectory(0)
h_genMatchedPt_V = TH1F("genpt", "gen Kshort pT; [GeV];", 50, 0, 50); h_genMatchedPt_V.SetDirectory(0)
h_genMatchedPt_M = TH1F("genpt", "gen Kshort pT; [GeV];", 50, 0, 50); h_genMatchedPt_M.SetDirectory(0)

hV_reco_dr = TH1F("dr", "dr;", 100, 0, 1); hV_reco_dr.SetDirectory(0)
hM_reco_dr = TH1F("dr", "dr;", 100, 0, 1); hM_reco_dr.SetDirectory(0)

drQuarkJet = 0.5 #maching quark, jet by dr is right...?
drKsJet = 0.3
drGenjetRecojet = 0.3

for iEv, e in enumerate(f.Events):
    if iEv%1000 == 0: print "iEv:", iEv
    KsV = []
    KsM = []
    q, qb = None, None
    for j in xrange(e.nGenPart):
        if e.GenPart_status[j] > 30 or e.GenPart_status[j] < 20: continue
        if e.GenPart_pdgId[j] == 3:  q = j
        if e.GenPart_pdgId[j] == -3: qb = j
        if e.GenPart_pdgId[j] == 5:  q = j
        if e.GenPart_pdgId[j] == -5: qb = j
    if q is None: continue

    qj = None
    qbj = None
    for j in xrange(e.nJet):
        if deltaR(e.Jet_eta[j], e.Jet_phi[j],
                  e.GenPart_eta[q], e.GenPart_phi[q]) < drQuarkJet:
            qj = j
        elif deltaR(e.Jet_eta[j], e.Jet_phi[j],
                    e.GenPart_eta[qb], e.GenPart_phi[qb]) < drQuarkJet:
            qbj = j
    if qj is None or qbj is None: continue

    gqj = None
    gqbj = None
    for j in xrange(e.nGenJet):
        if deltaR(e.GenJet_eta[j], e.GenJet_phi[j],
                  e.GenPart_eta[q], e.GenPart_phi[q]) < drQuarkJet:
            gqj = j
        elif deltaR(e.GenJet_eta[j], e.GenJet_phi[j],
                    e.GenPart_eta[qb], e.GenPart_phi[qb]) < drQuarkJet:
            gqbj = j
    if gqj is None or gqbj is None: continue

    # V0 producer
    vL2D = lambda j: sqrt(e.Kshort_x[j]*e.Kshort_x[j] + e.Kshort_y[j]*e.Kshort_y[j])
    vL3D = lambda j: sqrt(e.Kshort_x[j]*e.Kshort_x[j] + e.Kshort_y[j]*e.Kshort_y[j] + e.Kshort_z[j]*e.Kshort_z[j])
    for j in xrange(e.nKshort):
        dr = deltaR(e.Kshort_eta[j], e.Kshort_phi[j],
                    e.Jet_eta[qj], e.Jet_phi[qj])
        drb = deltaR(e.Kshort_eta[j], e.Kshort_phi[j],
                     e.Jet_eta[qbj], e.Jet_phi[qbj])

        if   dr  < drKsJet: x = e.Kshort_pt[j] / e.Jet_pt[qj]
        elif drb < drKsJet: x = e.Kshort_pt[j] / e.Jet_pt[qbj]
        else: continue

        #if x < 0.15: continue
        #if abs(e.Kshort_mass[j]-0.49) > 0.03: continue
        KsV.append(j)

        hV.Fill(e.Kshort_mass[j])
        #hV_reco_dr.Fill(mindr)

    # meson producer
    for j in xrange(e.nmeson):
        if e.meson_pdgId[j] != 310: continue
        dr = deltaR(e.meson_eta[j], e.meson_phi[j],
                    e.Jet_eta[qj], e.Jet_phi[qj])
        drb = deltaR(e.meson_eta[j], e.meson_phi[j],
                     e.Jet_eta[qbj], e.Jet_phi[qbj])

        if   dr  < drKsJet: x = e.meson_pt[j] / e.Jet_pt[qj]
        elif drb < drKsJet: x = e.meson_pt[j] / e.Jet_pt[qbj]
        else: continue

        #if x < 0.15: continue
        #if abs(e.meson_mass[j]-0.49) > 0.03: continue
        KsM.append(j)

        hM.Fill(e.meson_mass[j])
        #hM_reco_dr.Fill(mindr)

    #matchig
    for j in xrange(e.nV0GenPart):
        if e.V0GenPart_pdgId[j] !=310: continue
        dr_gen = deltaR(e.V0GenPart_eta[j], e.V0GenPart_phi[j],
                        e.GenJet_eta[gqj], e.GenJet_phi[gqj])
        drb_gen = deltaR(e.V0GenPart_eta[j], e.V0GenPart_phi[j],
                         e.GenJet_eta[gqbj], e.GenJet_phi[gqbj])

        if   dr_gen  < drKsJet: x = e.V0GenPart_pt[j] / e.GenJet_pt[gqj]
        elif drb_gen < drKsJet: x = e.V0GenPart_pt[j] / e.GenJet_pt[gqbj]
        else: continue

        #if x < 0.15: continue
        #if abs(e.V0GenPart_mass[j]-0.49) > 0.03: continue
        h_genPt.Fill(e.V0GenPart_pt[j])

        matchV = False
        for i in KsV:
            dr_withgen = deltaR(e.Kshort_eta[i], e.Kshort_phi[i],
                                e.V0GenPart_eta[j], e.V0GenPart_phi[j])
            if dr_withgen < drGenjetRecojet: matchV = True; break;
        if matchV: h_genMatchedPt_V.Fill(e.V0GenPart_pt[j])

        matchM = False
        for i in KsM:
            dr_withgen = deltaR(e.meson_eta[i], e.meson_phi[i],
                                e.V0GenPart_eta[j], e.V0GenPart_phi[j])
            if dr_withgen < drGenjetRecojet: matchM = True; break;
        if matchM: h_genMatchedPt_M.Fill(e.V0GenPart_pt[j])

c=TCanvas()

#hV_reco_dr.SetLineColor(2)
#hM_reco_dr.Scale(1/hM_reco_dr.Integral())
#hV_reco_dr.Scale(1/hV_reco_dr.Integral())
#hM_reco_dr.Draw()
#hV_reco_dr.Draw("same")
#c.Print("dr%s.png" % ext)

h_genPt.Draw()
h_genPt.SetLineColor(1)
h_genMatchedPt_V.Draw("same")
h_genMatchedPt_V.SetLineColor(2)
h_genMatchedPt_M.Draw("same")
c.Print("ksPt%s.png" % ext)
h_genMatchedPt_M.SetMinimum(0)
h_genMatchedPt_M.Divide(h_genPt)
h_genMatchedPt_M.Draw()
h_genMatchedPt_V.Divide(h_genPt)
h_genMatchedPt_V.Draw("same")
c.Print("eff%s.png" % ext)

