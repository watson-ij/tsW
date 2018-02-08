import ROOT
from ROOT import TFile, TLorentzVector
from glob import glob

# Must have loaded the CMSSW environment to read the files properly
ROOT.gSystem.Load("libPhysicsToolsFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable()

f = TFile("result/tmp.djMzv2e7Ml/step3_inMINIAODSIM.root")

isS = None
for i in range(f.Events.GetEntries()):
    _ = f.Events.GetEntry(i)
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
    KS = f.Events.recoVertexCompositePtrCandidates_slimmedKshortVertices__RECO.product()
    TLk = TLorentzVector()
    for k in KS:
        k = k.p4()
        TLk.SetPtEtaPhiE(k.pt(), k.eta(), k.phi(), k.e())
        if TLk.DeltaR(TLs4) < 0.5:
            print "Match K", k.pt() / s4.pt(), k.mass()
    # EDM
    # Lam = f.Events.recoVertexCompositeCandidates_generalV0Candidates_Lambda_RECO.product()
    # MiniAOD
    Lam = f.Events.recoVertexCompositePtrCandidates_slimmedLambdaVertices__RECO.product()
    for k in Lam:
        k = k.p4()
        TLk.SetPtEtaPhiE(k.pt(), k.eta(), k.phi(), k.e())
        if TLk.DeltaR(TLs4) < 0.5:
            print "Match L", k.pt() / s4.pt(), k.mass()
    
