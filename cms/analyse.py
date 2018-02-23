import ROOT
from ROOT import TFile, TLorentzVector
from glob import glob

# Must have loaded the CMSSW environment to read the files properly
ROOT.gSystem.Load("libPhysicsToolsFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable()

fN = TFile("result/20180219_171825_cms-t3-wn3022.sdfarm.kr_ttsWbW/nanoAOD.root")
f = TFile("result/20180219_171825_cms-t3-wn3022.sdfarm.kr_ttsWbW/step3_inMINIAODSIM.root")

isS = None
for i in range(f.Events.GetEntries()):
    _ = f.Events.GetEntry(i)
    _ = fN.Events.GetEntry(i)
    s4 = None
    ## GenParticle, search for our s-quark
    # EDM
    # hepmcData = f.Events.edmHepMCProduct_generatorSmeared__SIM.product().getHepMCData()
    # genParticles = f.Events.recoGenParticles_genParticles__HLT.product()
    # MiniAOD
    genParticles = f.Events.recoGenParticles_prunedGenParticles__RECO.product()
    for g in genParticles:
        if 20 < g.status() < 30 and g.pdgId() == 3:
            s4 = g.p4()
            isS = True
        elif 20 < g.status() < 30 and g.pdgId() == 5:
            s4 = g.p4()
            isS = False
    if s4 is None:
        print "BAD EVENT"
        continue
    TLs4 = TLorentzVector(); TLs4.SetPtEtaPhiE(s4.pt(), s4.eta(), s4.phi(), s4.e())
    ### Reconstructed KS/Lambda collections
    # EDM
    # KS = f.Events.recoVertexCompositeCandidates_generalV0Candidates_Kshort_RECO.product()
    # MiniAOD
    nKSL = 0
    KS = f.Events.recoVertexCompositePtrCandidates_slimmedKshortVertices__RECO.product()
    TLk = TLorentzVector()
    mK = []
    for k in KS:
        kl = k.p4()
        TLk.SetPtEtaPhiM(kl.pt(), kl.eta(), kl.phi(), kl.mass())
        if TLk.DeltaR(TLs4) < 0.5:
            print k.numberOfDaughters(), k.daughter(0), k.daughter(0).pt(), k.daughter(0).eta(), k.daughter(0).phi(), k.daughter(1), k.daughter(1).pt(), k.daughter(1).eta(), k.daughter(1).phi(), "V", k.vx(), k.vy(), k.vz() #, dir(k)
            print "Match K", k.pt() / s4.pt(), k.mass(), isS
            nKSL += 1
            mK.append(TLk.Clone())
    # EDM
    # Lam = f.Events.recoVertexCompositeCandidates_generalV0Candidates_Lambda_RECO.product()\
    # MiniAOD
    Lam = f.Events.recoVertexCompositePtrCandidates_slimmedLambdaVertices__RECO.product()
    for k in Lam:
        k = k.p4()
        TLk.SetPtEtaPhiE(k.pt(), k.eta(), k.phi(), k.e())
        if TLk.DeltaR(TLs4) < 0.5:
            print "Match L", k.pt() / s4.pt(), k.mass(), isS
            nKSL += 1
    if nKSL > 0: print "Total KS/Lambda: ", nKSL
    if len(mK) == 0: continue
    jtl = TLorentzVector()
    Jets = f.Events.patJets_slimmedJets__RECO.product()
    for j in Jets:
        jtl.SetPtEtaPhiM(j.pt(), j.eta(), j.phi(), j.mass())
        if any(jtl.DeltaR(m) < 0.5 for m in mK): print "J", jtl.Pt(), jtl.Eta(), jtl.Phi(), "Q", TLs4.Pt(), TLs4.Eta(), TLs4.Phi()
        else: continue
        # print dir(j)
        l = None
        for m in mK:
            print " KS", m.Pt(), m.Eta(), m.Phi(), m.M()
        ptrs = []
        ntrs = []
        ptl = TLorentzVector()
        for p in j.getJetConstituents():
            # print " C", p, p.charge(), p.pt(), p.eta(), p.phi()
            #
            # print " C", p.pdgId(), p.charge(), p.isIsolatedChargedHadron(), p.hasTrackDetails(), p.mass()
            ptl.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass())
            if p.charge() == 1: ptrs.append(ptl.Clone())
            elif p.charge() == -1: ntrs.append(ptl.Clone())
            # print dir(p.get())
        for p in ptrs:
            for n in ntrs:
                cm = p+n
                if 0.43 < cm.M() < 0.57:
                    print " CC", cm.Pt(), cm.Eta(), cm.Phi(), cm.M(), ",", p.Pt(), n.Pt()
    for n in xrange(fN.Events.ncmeson):
        if fN.Events.cmeson_pdgId[n] == 310:
            print " N", fN.Events.cmeson_pt[n], fN.Events.cmeson_eta[n], fN.Events.cmeson_phi[n], fN.Events.cmeson_mass[n]
        # exit()
        # for p in j.getPFConstituents():
        #     print " P", p.pdgId()
        # exit()
