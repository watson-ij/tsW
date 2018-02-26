#!/usr/bin/env python2

import ROOT
from ROOT import TFile, TCanvas, TTree, TH1F
from math import hypot, pi, sqrt
ROOT.gROOT.SetBatch(True)

ext = ""
f = TFile("cmssw_ttswbw.root")

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
nVC = collections.defaultdict(lambda: 0)
hC = TH1F("C", "C Meson collection;Mass [GeV];", 50, 0.43, 0.57); hC.SetDirectory(0)
hCc = TH1F("Cc", "C Meson collection;Mass [GeV];", 50, 0.43, 0.57); hCc.SetDirectory(0)
hV = TH1F("V", "V0Meson collection;Mass [GeV];", 50, 0.43, 0.57); hV.SetDirectory(0)
hM = TH1F("M", "Meson collection;Mass [GeV];", 50, 0.43, 0.57); hM.SetDirectory(0)
hMc = TH1F("Mc", "Meson collection;Mass [GeV];", 50, 0.43, 0.57); hMc.SetDirectory(0)

hC2D = TH1F("C2D", "C Meson collection;Decay Length 2D [cm];", 50, 0., 10); hC2D.SetDirectory(0)
hV2D = TH1F("V2D", "V0Meson collection;Decay Length 2D [cm];", 50, 0., 10); hC2D.SetDirectory(0)
hM2D = TH1F("M2D", "Meson collection;Decay Length 2D [cm];", 50, 0., 10); hM2D.SetDirectory(0)

hC3D = TH1F("C3D", "C Meson collection;Decay Length 3D [cm];", 50, 0., 25); hC3D.SetDirectory(0)
hV3D = TH1F("V3D", "V0Meson collection;Decay Length 3D [cm];", 50, 0., 25); hC3D.SetDirectory(0)
hM3D = TH1F("M3D", "Meson collection;Decay Length 3D [cm];", 50, 0., 25); hM3D.SetDirectory(0)

hCx = TH1F("Cx", "C Meson collection;x;", 50, 0., 1.); hCx.SetDirectory(0)
hVx = TH1F("Vx", "V0Meson collection;x;", 50, 0., 1.); hVx.SetDirectory(0)
hMx = TH1F("Mx", "Meson collection;x;", 50, 0., 1.); hMx.SetDirectory(0)


hCj = TH1F("Cj", "C Meson collection;j;", 50, 0., 1.); hCx.SetDirectory(0)
hVj = TH1F("Vx", "V0Meson collection;j;", 50, 0., 1.); hVx.SetDirectory(0)
hMj = TH1F("Mx", "Meson collection;j;", 50, 0., 1.); hMx.SetDirectory(0)

for iEv, e in enumerate(f.Events):
    nS = nB = 0
    q, qb = None, None
    isS = False
    for j in xrange(f.Events.nGenPart):
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == 3:  q = j; nS += 1; isS = True
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == -3: qb = j; nS += 1
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == 5:  q = j; nB += 1
        if 20 < f.Events.GenPart_status[j] < 30 and f.Events.GenPart_pdgId[j] == -5: qb = j; nB += 1
    if q is None: continue
    qj = None
    qbj = None
    for j in xrange(f.Events.nJet):
        if deltaR(f.Events.Jet_eta[j], f.Events.Jet_phi[j],
                  f.Events.GenPart_eta[q], f.Events.GenPart_phi[q]) < 0.5:
            qj = j
        elif deltaR(f.Events.Jet_eta[j], f.Events.Jet_phi[j],
                    f.Events.GenPart_eta[qb], f.Events.GenPart_phi[qb]) < 0.5:
            qbj = j
    if qj is None or qbj is None: continue
    nSB[(nS, nB)] += 1
    match = []
    vL2D = lambda j: sqrt(f.Events.Kshort_x[j]*f.Events.Kshort_x[j] + f.Events.Kshort_y[j]*f.Events.Kshort_y[j])
    vL3D = lambda j: sqrt(f.Events.Kshort_x[j]*f.Events.Kshort_x[j] + f.Events.Kshort_y[j]*f.Events.Kshort_y[j] + f.Events.Kshort_z[j]*f.Events.Kshort_z[j])
    for j in xrange(f.Events.nKshort):
        dr = deltaR(f.Events.Kshort_eta[j], f.Events.Kshort_phi[j],
#                  f.Events.GenPart_eta[q], f.Events.GenPart_phi[q]) < 0.5:
                  f.Events.Jet_eta[qj], f.Events.Jet_phi[qj])
        drb = deltaR(f.Events.Kshort_eta[j], f.Events.Kshort_phi[j],
                     # f.Events.GenPart_eta[qb], f.Events.GenPart_phi[qb]) < 0.5:
                     f.Events.Jet_eta[qbj], f.Events.Jet_phi[qbj])
        if dr < 0.3:
            if f.Events.Kshort_pt[j] / f.Events.Jet_pt[qj] < 0.15: continue
            if isS:
                nSM += 1; match.append(j)
                hV.Fill(f.Events.Kshort_mass[j])
                if 0.48 < f.Events.Kshort_mass[j] < 0.52:
                    hV2D.Fill(vL2D(j))
                    hVx.Fill(f.Events.Kshort_pt[j] / f.Events.Jet_pt[qj])
                    hVj.Fill(dr)
            else:
                nBM += 1
                hV.Fill(f.Events.Kshort_mass[j])
                if 0.48 < f.Events.Kshort_mass[j] < 0.52:
                    hV2D.Fill(vL2D(j))
                    hVx.Fill(f.Events.Kshort_pt[j] / f.Events.Jet_pt[qj])
                    hVj.Fill(dr)
        elif drb < 0.3:
            if f.Events.Kshort_pt[j] / f.Events.Jet_pt[qbj] < 0.15: continue
            if isS:
                nSM += 1; match.append(j)
                hV.Fill(f.Events.Kshort_mass[j])
                if 0.48 < f.Events.Kshort_mass[j] < 0.52:
                    hV2D.Fill(vL2D(j))
                    hV3D.Fill(vL3D(j))
                    hVx.Fill(f.Events.Kshort_pt[j] / f.Events.Jet_pt[qbj])
                    hVj.Fill(drb)
            else:
                nBM += 1
                hV.Fill(f.Events.Kshort_mass[j])
                if 0.48 < f.Events.Kshort_mass[j] < 0.52:
                    hV2D.Fill(vL2D(j))
                    hV3D.Fill(vL3D(j))
                    hVx.Fill(f.Events.Kshort_pt[j] / f.Events.Jet_pt[qbj])
                    hVj.Fill(drb)


    vL2D = lambda j: sqrt(f.Events.meson_x[j]*f.Events.meson_x[j] + f.Events.meson_y[j]*f.Events.meson_y[j])
    vL3D = lambda j: sqrt(f.Events.meson_x[j]*f.Events.meson_x[j] + f.Events.meson_y[j]*f.Events.meson_y[j] + f.Events.meson_z[j]*f.Events.meson_z[j])
    for j in xrange(f.Events.nmeson):
        if f.Events.meson_pdgId[j] != 310: continue
        dr = deltaR(f.Events.meson_eta[j], f.Events.meson_phi[j],
#                  f.Events.GenPart_eta[q], f.Events.GenPart_phi[q]) < 0.5:
                  f.Events.Jet_eta[qj], f.Events.Jet_phi[qj])
        drb = deltaR(f.Events.meson_eta[j], f.Events.meson_phi[j],
                     # f.Events.GenPart_eta[qb], f.Events.GenPart_phi[qb]) < 0.5:
                     f.Events.Jet_eta[qbj], f.Events.Jet_phi[qbj])
        if dr < 0.3:
            if f.Events.meson_pt[j] / f.Events.Jet_pt[qj] < 0.15: continue
            if isS:
                nSM += 1; match.append(j)
                hM.Fill(f.Events.meson_mass[j])
                if abs(f.Events.meson_chi2[j]) < 1.0 and f.Events.meson_dca[j]<0.2 and f.Events.meson_lxy[j]>0.25 and f.Events.meson_angleXY[j] > 0:
                    hMc.Fill(f.Events.meson_mass[j])
                    if 0.48 < f.Events.meson_mass[j] < 0.52:
                        hM2D.Fill(vL2D(j))
                        hMx.Fill(f.Events.meson_pt[j] / f.Events.Jet_pt[qj])
                        hMj.Fill(dr)
            else:
                nBM += 1
                hM.Fill(f.Events.meson_mass[j])
                if abs(f.Events.meson_chi2[j]) < 1.0 and f.Events.meson_dca[j]<0.2 and f.Events.meson_lxy[j]>0.25 and f.Events.meson_angleXY[j] > 0:
                    hMc.Fill(f.Events.meson_mass[j])
                    if 0.48 < f.Events.meson_mass[j] < 0.52:
                        hM2D.Fill(vL2D(j))
                        hMx.Fill(f.Events.meson_pt[j] / f.Events.Jet_pt[qj])
                        hMj.Fill(dr)
        elif drb < 0.3:
            if f.Events.meson_pt[j] / f.Events.Jet_pt[qbj] < 0.15: continue
            if isS:
                nSM += 1; match.append(j)
                hM.Fill(f.Events.meson_mass[j])
                if abs(f.Events.meson_chi2[j]) < 1.0 and f.Events.meson_dca[j]<0.2 and f.Events.meson_lxy[j]>0.25 and f.Events.meson_angleXY[j] > 0:
                    hMc.Fill(f.Events.meson_mass[j])
                    if 0.48 < f.Events.meson_mass[j] < 0.52:
                        hM2D.Fill(vL2D(j))
                        hM3D.Fill(vL3D(j))
                        hMx.Fill(f.Events.meson_pt[j] / f.Events.Jet_pt[qbj])
                        hMj.Fill(drb)
            else:
                nBM += 1
                hM.Fill(f.Events.meson_mass[j])
                if abs(f.Events.meson_chi2[j]) < 1.0 and f.Events.meson_dca[j]<0.2 and f.Events.meson_lxy[j]>0.25 and f.Events.meson_angleXY[j] > 0:
                    hMc.Fill(f.Events.meson_mass[j])
                    if 0.48 < f.Events.meson_mass[j] < 0.52:
                        hM2D.Fill(vL2D(j))
                        hM3D.Fill(vL3D(j))
                        hMx.Fill(f.Events.meson_pt[j] / f.Events.Jet_pt[qbj])
                        hMj.Fill(drb)

                    # for j in xrange(f.Events.nV0GenPart):
    #     for m in match:
    #         if m >= 0 and deltaR(f.Events.Kshort_eta[m], f.Events.Kshort_phi[m],
    #                              f.Events.V0GenPart_eta[j], f.Events.V0GenPart_phi[j]) < 0.1 and \
    #                              abs(1.0 - f.Events.Kshort_pt[m] / f.Events.V0GenPart_pt[j]) < 0.1:
    #             pass
    cmeson = []
    cL2D = lambda j: sqrt(f.Events.cmeson_x[j]*f.Events.cmeson_x[j] + f.Events.cmeson_y[j]*f.Events.cmeson_y[j])
    cL3D = lambda j: sqrt(f.Events.cmeson_x[j]*f.Events.cmeson_x[j] + f.Events.cmeson_y[j]*f.Events.cmeson_y[j] + f.Events.cmeson_z[j]*f.Events.cmeson_z[j])
    for j in xrange(f.Events.ncmeson):
        if f.Events.cmeson_angleXY[j] < 0: continue
        if f.Events.cmeson_pdgId[j] == 310:
            dr =  deltaR(f.Events.cmeson_eta[j], f.Events.cmeson_phi[j],
                      # f.Events.GenPart_eta[q], f.Events.GenPart_phi[q]) < 0.5:
                         f.Events.Jet_eta[qj], f.Events.Jet_phi[qj])
            drb = deltaR(f.Events.cmeson_eta[j], f.Events.cmeson_phi[j],
                        # f.Events.GenPart_eta[qb], f.Events.GenPart_phi[qb]) < 0.5:
                        f.Events.Jet_eta[qbj], f.Events.Jet_phi[qbj])
            if dr < 0.3:
                if f.Events.cmeson_pt[j] / f.Events.Jet_pt[qj] < 0.15: continue
                cmeson.append(j)
                hC.Fill(f.Events.cmeson_mass[j])
                if abs(f.Events.cmeson_chi2[j]) < 1.0 and f.Events.cmeson_dca[j]<0.2 and f.Events.cmeson_lxy[j]>0.25 and f.Events.cmeson_angleXY[j] > 0:
                    hCc.Fill(f.Events.cmeson_mass[j])
                    if 0.48 < f.Events.cmeson_mass[j] < 0.52:
                        hC2D.Fill(cL2D(j))
                        hC3D.Fill(cL3D(j))
                        hCx.Fill(f.Events.cmeson_pt[j] / f.Events.Jet_pt[qj])
                        hCj.Fill(dr)
            elif drb < 0.3:
                if f.Events.cmeson_pt[j] / f.Events.Jet_pt[qbj] < 0.15: continue
                cmeson.append(j)
                hC.Fill(f.Events.cmeson_mass[j])
                if abs(f.Events.cmeson_chi2[j]) < 1.0 and f.Events.cmeson_dca[j]<0.2 and f.Events.cmeson_lxy[j]>0.25 and f.Events.cmeson_angleXY[j] > 0:
                    hCc.Fill(f.Events.cmeson_mass[j])
                    if 0.48 < f.Events.cmeson_mass[j] < 0.52:
                        hC2D.Fill(cL2D(j))
                        hC3D.Fill(cL3D(j))
                        hCx.Fill(f.Events.cmeson_pt[j] / f.Events.Jet_pt[qbj])
                        hCj.Fill(drb)
    nVC[(len(match), len(cmeson))] += 1

for i in nVC: print i, nVC[i],

c=TCanvas()

hC.SetMinimum(0)
hC.Draw()
hCc.SetLineColor(ROOT.kGreen)
hCc.Draw("same")
hV.SetLineColor(ROOT.kRed)
hV.Draw("same")
hM.SetLineColor(ROOT.kCyan+3)
hM.Draw("same")
hMc.SetLineColor(ROOT.kMagenta+1)
hMc.Draw("same")
c.Print("c%s.png" % ext)
hV.SetMinimum(0)
hV.Draw()
c.Print("v%s.png" % ext)

l = ROOT.TLegend(0.5,0.4,0.95,0.7)

hC2D.SetMinimum(0)
hC2D.SetMaximum(200)
hC2D.Draw()
hV2D.SetLineColor(ROOT.kRed)
hV2D.Draw("same")
hM2D.SetLineColor(ROOT.kMagenta+1)
hM2D.Draw("same")
l.AddEntry(hC2D, "CMesonProducer Collection")
l.AddEntry(hV2D, "V0Producer Collection")
l.AddEntry(hM2D, "MesonProducer Collection")
l.Draw()
c.Print("c2D%s.png" % ext)
hV2D.SetMinimum(0)
hV2D.Draw()
c.Print("v2D%s.png" % ext)

hC3D.SetMinimum(0)
hC3D.SetMaximum(200)
hC3D.Draw()
hV3D.SetLineColor(ROOT.kRed)
hV3D.Draw("same")
hM3D.SetLineColor(ROOT.kMagenta+1)
hM3D.Draw("same")
c.Print("c3D%s.png" % ext)
hV3D.SetMinimum(0)
hV3D.Draw()
c.Print("v3D%s.png" % ext)

hC2D.SetMinimum(0.01)
hC2D.Draw()
c.SetLogy(True)
c.Print("c2D_l%s.png" % ext)
hV2D.SetMinimum(0.01)
hV2D.Draw()
c.Print("v2D_l%s.png" % ext)

hCx.SetMinimum(0.01)
hCx.Draw()
hVx.SetLineColor(ROOT.kRed)
hVx.Draw("same")
hMx.SetLineColor(ROOT.kMagenta+1)
hMx.Draw("same")
#l.Draw()
c.Print("cx%s.png" % ext)

c.SetLogy(False)

hCj.SetMinimum(0.0)
hCj.Draw()
hVj.SetLineColor(ROOT.kRed)
hVj.Draw("same")
hMj.SetLineColor(ROOT.kMagenta+1)
hMj.Draw("same")
# l.Draw()
c.Print("cj%s.png" % ext)
