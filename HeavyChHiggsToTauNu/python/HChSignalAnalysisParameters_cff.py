import FWCore.ParameterSet.Config as cms

trigger = cms.untracked.PSet(
    src = cms.untracked.InputTag("patTriggerEvent"),
    trigger = cms.untracked.string("HLT_SingleLooseIsoTau20") # in 36X/35X MC and Run2010A data
#    trigger = cms.untracked.string("HLT_SingleIsoTau20_Trk5_MET20") # in 38X MC and Run2010B data
)

tauSelectionBase = cms.untracked.PSet(
    src = cms.untracked.InputTag("selectedPatTausShrinkingConePFTau"),
    selection = cms.untracked.string(""),
    ptCut = cms.untracked.double(30),
    etaCut = cms.untracked.double(2.4), #no change
    leadingTrackPtCut = cms.untracked.double(20),
    rtauCut = cms.untracked.double(0.8), #no change
    invMassCut = cms.untracked.double(1.5) #no change
)

tauSelectionCaloTauCutBased = tauSelectionBase.clone()
tauSelectionCaloTauCutBased.src = cms.untracked.InputTag("selectedPatTausCaloRecoTau")
tauSelectionCaloTauCutBased.selection = cms.untracked.string("CaloTauCutBased")

tauSelectionShrinkingConeCutBased = tauSelectionBase.clone()
tauSelectionShrinkingConeCutBased.src = cms.untracked.InputTag("selectedPatTausShrinkingConePFTau")
tauSelectionShrinkingConeCutBased.selection = cms.untracked.string("ShrinkingConePFTauCutBased")

tauSelectionShrinkingConeTaNCBased = tauSelectionBase.clone()
tauSelectionShrinkingConeTaNCBased.src = cms.untracked.InputTag("selectedPatTausShrinkingConePFTau")
tauSelectionShrinkingConeTaNCBased.selection = cms.untracked.string("ShrinkingConePFTauTaNCBased")

tauSelectionHPSTauBased = tauSelectionBase.clone()
tauSelectionHPSTauBased.src = cms.untracked.InputTag("selectedPatTausHpsPFTau")
tauSelectionHPSTauBased.selection = cms.untracked.string("HPSTauBased")

tauSelection = tauSelectionShrinkingConeCutBased

jetSelection = cms.untracked.PSet(
    src = cms.untracked.InputTag("selectedPatJets"),
    #src = cms.untracked.InputTag("selectedPatJetsAK5JPT"),
    cleanTauDR = cms.untracked.double(0.5), #no change
    ptCut = cms.untracked.double(30),
    etaCut = cms.untracked.double(2.4), #no change
    minNumber = cms.untracked.uint32(3)
)

MET = cms.untracked.PSet(
    #src = cms.untracked.InputTag("patMETs"), # calo MET
    src = cms.untracked.InputTag("patMETsPF"), # PF MET
    #src = cms.untracked.InputTag("patMETsTC"), # tc MET
    METCut = cms.untracked.double(50.0)
)

bTagging = cms.untracked.PSet(
    discriminator = cms.untracked.string("trackCountingHighEffBJetTags"),
    discriminatorCut = cms.untracked.double(1.5), #no change
    ptCut = cms.untracked.double(30), #no change
    etaCut = cms.untracked.double(1.5), #change to 2.4 (was 1.5)
    minNumber = cms.untracked.uint32(1)
)

EvtTopology = cms.untracked.PSet(
    #discriminator = cms.untracked.string("test"),
    #discriminatorCut = cms.untracked.double(0.0),
    #alphaT = cms.untracked.double(-5.00)
    alphaT = cms.untracked.double(-5.0)
)
