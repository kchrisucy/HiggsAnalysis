DataCardName    = 'myDummyTestName'
#Path            = '/mnt/flustre/slehti/HplusDataForLands'
Path            = '/home/wendland/data/v25b/freeze_forPREAPP'
#Path            = '/home/wendland/data/v25b/test3loose'
#MassPoints      = [80,90,100,120,140,150,155,160]
MassPoints      = [80,90,100,120,140,150,155,160]
#MassPoints      = [80,120,160]
#MassPoints      = [120]

BlindAnalysis   = True

# Specify name of EDFilter or EDAnalyser process that produced the root files
SignalAnalysis  = "signalAnalysis"
#EmbeddingAnalysis     = "signalAnalysisCaloMet60TEff"
EmbeddingAnalysis     = "signalAnalysisCaloMet60"
QCDFactorisedAnalysis = "QCDMeasurement"
QCDInvertedAnalysis = None # FIXME

RootFileName    = "histograms.root"

# Rate counter definitions
SignalRateCounter = "Selected events"
FakeRateCounter = "EWKfaketaus:SelectedEvents"

# Options
OptionMassShape = "TransverseMass"
#OptionMassShape = "FullMass"
OptionReplaceEmbeddingByMC = False
OptionIncludeSystematics = True # Set to true if the JES and PU uncertainties were produced
OptionPurgeReservedLines = True # Makes limit running faster, but cannot combine leptonic datacards
OptionDoControlPlots = True

# Options for reports and article
OptionBr = 0.03  # Br(t->bH+)

# Shape histogram definitions
SignalShapeHisto = ""
FakeShapeHisto = ""
ShapeHistogramsDimensions = {}

if OptionMassShape == "TransverseMass":
    SignalShapeHisto = "transverseMass"
    FakeShapeHisto = "EWKFakeTausTransverseMass"
    ShapeHistogramsDimensions = { "bins": 20,
                                  "rangeMin": 0.0,
                                  "rangeMax": 400.0,
                                  "variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                                  #"variableBinSizeLowEdges": [0,20,40,60,80,100,120,140,160,200], # if an empty list is given, then uniform bin width is used
                                  "xtitle": "Transverse mass / GeV",
                                  "ytitle": "Events" }
elif OptionMassShape == "FullMass":
    SignalShapeHisto = "FullHiggsMass/HiggsMass"
    FakeShapeHisto = "FIXME" #### FIXME
    ShapeHistogramsDimensions = { "bins": 25,
                                  "rangeMin": 0.0,
                                  "rangeMax": 500.0,
                                  "variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                                  #"variableBinSizeLowEdges": [0,20,40,60,80,100,120,140,160,200], # if an empty list is given, then uniform bin width is used
                                  "xtitle": "Full mass / GeV",
                                  "ytitle": "Events" }

DataCardName += "_"+OptionMassShape

QCDFactorisedMETShapeHistogramsDimensions = {  "bins": 10,
                                               "rangeMin": 0.0,
                                               "rangeMax": 100.0,
                                               #"variableBinSizeLowEdges": [0,20,40,60,80,100,120,140,160,180,200,250,300], # if an empty list is given, then uniform bin width is used
                                               "variableBinSizeLowEdges": [0,10,20,30,40,50,60,70,80,90], # if an empty list is given, then uniform bin width is used
                                               "xtitle": "E_{T}^{miss}, GeV/c^{2}",
                                               "ytitle": "Events"}

##############################################################################
# Observation definition (how to retrieve number of observed events)
#
from HiggsAnalysis.HeavyChHiggsToTauNu.datacardtools.InputClasses import ObservationInput
Observation = ObservationInput(dirPrefix=SignalAnalysis,
                               rateCounter=SignalRateCounter,
                               datasetDefinitions=["Tau_"],
                               shapeHisto=SignalShapeHisto)
#Observation.setPaths(signalPath,signalDataPaths)

##############################################################################
# DataGroup (i.e. columns in datacard) definitions
#
from HiggsAnalysis.HeavyChHiggsToTauNu.datacardtools.InputClasses import DataGroup
DataGroups = []
EmbeddingIdList = []
EWKFakeIdList = []

signalTemplate = DataGroup(datasetType="Signal",
                           shapeHisto=SignalShapeHisto,
                           dirPrefix=SignalAnalysis,
                           rateCounter=SignalRateCounter)

for mass in MassPoints:
    myMassList = [mass]
    hhx = signalTemplate.clone()
    hhx.setLabel("HH"+str(mass)+"_a")
    hhx.setLandSProcess(-1)
    hhx.setValidMassPoints(myMassList)
    hhx.setNuisances(["01","02","03","45","46","47","09","10","17","28","33","34"])
    hhx.setDatasetDefinitions(["TTToHplusBHminusB_M"+str(mass)]),
    DataGroups.append(hhx)

    hwx = signalTemplate.clone()
    hwx.setLabel("HW"+str(mass)+"_a")
    hwx.setLandSProcess(0)
    hwx.setValidMassPoints(myMassList)
    hwx.setNuisances(["01","02","03","45","46","47","09","10","18","28","33","34"])
    hwx.setDatasetDefinitions(["TTToHplusBWB_M"+str(mass)]),
    DataGroups.append(hwx)


if OptionMassShape == "TransverseMass":
    DataGroups.append(DataGroup(
        label        = "QCDfact",
        landsProcess = 3,
        validMassPoints = MassPoints,
        dirPrefix   = QCDFactorisedAnalysis,
        datasetType  = "QCD factorised",
        datasetDefinitions = ["Tau_"],
        MCEWKDatasetDefinitions = ["TTJets","WJets","DY","WW","WZ","ZZ","T_","Tbar_"],
        #MCEWKDatasetDefinitions = ["TTJets","W2Jets","DY","WW","WZ","ZZ","T_","Tbar_"],
        nuisances    = ["12","13","40b"],
        QCDfactorisedInfo = { "afterBigboxSource": "factorisation/AfterJetSelection",
                              "afterMETLegSource": "factorisation/Leg1AfterTopSelection",
                              "afterTauLegSource": "factorisation/Leg2AfterTauID",
                              "validationMETShapeSource": ["shape_CtrlLeg1METAfterStandardSelections/CtrlLeg1METAfterStandardSelections",
                                                           "shape_CtrlLeg1METAfterStandardSelectionsMET20/CtrlLeg1METAfterStandardSelectionsMET20",
                                                           "shape_CtrlLeg1METAfterStandardSelectionsMET30/CtrlLeg1METAfterStandardSelectionsMET30",
                                                            "shape_CtrlLeg1METAfterTauIDNoRtau/CtrlLeg1METAfterTauIDNoRtau",
                                                            "shape_CtrlLeg1METAfterTauIDNoRtauMET20/CtrlLeg1METAfterTauIDNoRtauMET20",
                                                            "shape_CtrlLeg1METAfterTauIDNoRtauMET30/CtrlLeg1METAfterTauIDNoRtauMET30",
                                                            "shape_CtrlLeg1METAfterFullTauID/CtrlLeg1METAfterFullTauID",
                                                            "shape_CtrlLeg1METAfterFullTauIDMET20/CtrlLeg1METAfterFullTauIDMET20",
                                                            "shape_CtrlLeg1METAfterFullTauIDMET30/CtrlLeg1METAfterFullTauIDMET30",
                                                            ],
                              "validationMETShapeDetails": QCDFactorisedMETShapeHistogramsDimensions,
                              "basicMtHisto": "shape_MtShapesAfterFullMETLeg/MtShapesAfterFullMETLeg", # prefix for shape histograms in MET leg (will be weighted by tau leg efficiency)
                              "validationMtShapeSource": ["shape_MtShapesAfterStandardSelection/MtShapesAfterStandardSelection",
                                                          "shape_MtShapesAfterStandardSelectionMET20/MtShapesAfterStandardSelectionMET20",
                                                          "shape_MtShapesAfterStandardSelectionMET30/MtShapesAfterStandardSelectionMET30",
                                                           "shape_MtShapesAfterTauIDNoRtau/MtShapesAfterTauIDNoRtau",
                                                           "shape_MtShapesAfterTauIDNoRtauMET20/MtShapesAfterTauIDNoRtauMET20",
                                                           "shape_MtShapesAfterTauIDNoRtauMET30/MtShapesAfterTauIDNoRtauMET30",
                                                           "shape_MtShapesAfterTauID/MtShapesAfterTauID",
                                                           "shape_MtShapesAfterTauIDMET20/MtShapesAfterTauIDMET20",
                                                           "shape_MtShapesAfterTauIDMET30/MtShapesAfterTauIDMET30",
                                                           ],
                              "assumedMCEWKSystUncertainty": 0.20, # has no effect anymore ...
                              "factorisationMapAxisLabels": ["#tau p_{T}, GeV", "#tau #eta", "N_{vertices}"],
        }
    ))
elif OptionMassShape == "FullMass":
    DataGroups.append(DataGroup(
        label        = "QCDfact",
        landsProcess = 3,
        validMassPoints = MassPoints,
        dirPrefix   = QCDFactorisedAnalysis,
        datasetType  = "QCD factorised",
        datasetDefinitions = ["Tau_"],
        MCEWKDatasetDefinitions = ["TTJets","WJets","DY","WW","WZ","ZZ","T_","Tbar_"],
        #MCEWKDatasetDefinitions = ["TTJets","W2Jets","W3Jets","W4Jets","DY","WW","WZ","ZZ","T_","Tbar_"],
        nuisances    = ["12","13","40b"],
        QCDfactorisedInfo = { "afterBigboxSource": "factorisation/AfterJetSelection",
                              "afterMETLegSource": "factorisation/Leg1AfterTopSelection",
                              "afterTauLegSource": "factorisation/Leg2AfterTauID",
                              "validationMETShapeSource": ["shape_CtrlLeg1METAfterStandardSelections/CtrlLeg1METAfterStandardSelections",
                                                            "shape_CtrlLeg1METAfterTauIDNoRtau/CtrlLeg1METAfterTauIDNoRtau",
                                                            "shape_CtrlLeg1METAfterFullTauID/CtrlLeg1METAfterFullTauID"],
                              "basicMtHisto": "shape_FullMassShapesAfterFullMETLeg/FullMassShapesAfterFullMETLeg", # prefix for shape histograms in MET leg (will be weighted by tau leg efficiency)
                              "validationMtShapeSource": [],
                              "assumedMCEWKSystUncertainty": 0.20,
                              "factorisationMapAxisLabels": ["#tau p_{T}, GeV", "#tau #eta", "N_{vertices}"],
        }
    ))

DataGroups.append(DataGroup(
    label        = "QCDinv",
    landsProcess = 3,
    validMassPoints = MassPoints,
    datasetType  = "QCD inverted",
    datasetDefinitions   = "^Tau\S+",
    shapeHisto   = "mtSum",
    dirPrefix   = QCDInvertedAnalysis,
    rateCounter  = "deltaPhiTauMET160 limit",
    additionalNormalisation= 0.0066,
    nuisances    = ["41","42","43","44"] # FIXME: add shape stat, i.e. 40x
))

if not OptionReplaceEmbeddingByMC:
    # EWK + ttbar with genuine taus
    EmbeddingIdList = [4]
    DataGroups.append(DataGroup(
        label        = "EWK_Tau",
        landsProcess = 4,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinitions   = ["SingleMu"],
        dirPrefix   = EmbeddingAnalysis,
        rateCounter  = SignalRateCounter,
        validMassPoints = MassPoints,
        nuisances    = ["01b","03","14","15","16","19","40","48"]
        #nuisances    = ["01b","03","45","14","15","16","19","40"]
    ))
    
    #DataGroups.append(DataGroup(
        #label        = "res_DY.",
        #landsProcess = 5,
        #datasetType  = "None",
        #validMassPoints = MassPoints,
    #))
    #DataGroups.append(DataGroup(
        #label        = "res_WW.",
        #landsProcess = 6,
        #datasetType  = "None",
        #validMassPoints = MassPoints,
    #))
    #DataGroups.append(DataGroup(
        #label        = "EWK_DY",
        #landsProcess = 5,
        #shapeHisto   = SignalShapeHisto,
        #datasetType  = "Embedding",
        #datasetDefinitions   = ["DYJetsToLL"],
        #dirPrefix   = SignalAnalysis,
        #rateCounter  = SignalRateCounter,
        #validMassPoints = MassPoints,
        #nuisances    = ["01c","03","45","46","47","09","11b","15b","16b","31","33","34","24"]
    #))
    #DataGroups.append(DataGroup(
        #label        = "EWK_VV",
        #landsProcess = 6,
        #shapeHisto   = SignalShapeHisto,
        #datasetType  = "Embedding",
        #datasetDefinitions   = ["WW"], #,"WZ","ZZ"],
        #dirPrefix   = SignalAnalysis,
        #rateCounter  = SignalRateCounter,
        #validMassPoints = MassPoints,
        #nuisances    = ["01c","03","45","46","47","09","11b","15b","16b","32","33","34","27"]
    #))

    # EWK + ttbar with fake taus
    EWKFakeIdList = [1,5,6]
    DataGroups.append(DataGroup(
        label        = "EWK_tt_faketau",
        landsProcess = 1,
        shapeHisto   = FakeShapeHisto,
        datasetType  = "Signal",
        datasetDefinitions   = ["TTJets_"],
        dirPrefix   = SignalAnalysis,
        rateCounter  = FakeRateCounter,
        validMassPoints = MassPoints,
        nuisances    = ["01","02","04","45","46","47","09","10","28","33","34b","35"]
    ))
    DataGroups.append(DataGroup(
        label        = "EWK_W_faketau",
        landsProcess = 5,
        shapeHisto   = FakeShapeHisto,
        datasetType  = "Signal",
        datasetDefinitions   = ["WJets"],
        dirPrefix   = SignalAnalysis,
        rateCounter  = FakeRateCounter,
        validMassPoints = MassPoints,
        nuisances    = ["01","02","04","45","46","47","09","11","29","33","34b","37"]
    ))
    DataGroups.append(DataGroup(
        label        = "EWK_t_faketau",
        landsProcess = 6,
        shapeHisto   = FakeShapeHisto,
        datasetType  = "Signal",
        datasetDefinitions   = ["T_", "Tbar_"],
        dirPrefix   = SignalAnalysis,
        rateCounter  = FakeRateCounter,
        validMassPoints = MassPoints,
        nuisances    = ["01","02","04","45","46","47","09","10","30","33","34b","38"]
    ))
else:
    # Mimic embedding with MC analysis (introduces double counting of EWK fakes, but that should be small effect)
    EmbeddingIdList = [4,5,6]
    DataGroups.append(DataGroup(
        label        = "MC_EWK_Tau",
        landsProcess = 4,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Signal",
        #datasetDefinitions = ["TTJets", "W2Jets","W3Jets","W4Jets", "Tbar_", "T_"],
        datasetDefinitions = ["TTJets", "WJets", "Tbar_", "T_"],
        dirPrefix   = SignalAnalysis,
        rateCounter  = SignalRateCounter,
        validMassPoints = MassPoints,
        nuisances    = ["01b","03","45","46","47","14","15","16","19","40"]
    ))
    DataGroups.append(DataGroup(
        label        = "MC_EWK_DY",
        landsProcess = 5,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Signal",
        datasetDefinitions   = ["DYJetsToLL"],
        dirPrefix   = SignalAnalysis,
        rateCounter  = SignalRateCounter,
        validMassPoints = MassPoints,
        nuisances    = ["01","03","45","46","47","09","11","15b","16b","31","33","34","24"]
    ))
    DataGroups.append(DataGroup(
        label        = "MC_EWK_VV",
        landsProcess = 6,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Signal",
        datasetDefinitions   = ["WW","WZ","ZZ"],
        dirPrefix   = SignalAnalysis,
        rateCounter  = SignalRateCounter,
        validMassPoints = MassPoints,
        nuisances    = ["01","03","45","46","47","09","11","15b","16b","32","33","34","27"]
    ))
    DataGroups.append(DataGroup(
        label        = "empty",
        landsProcess = 1,
        datasetType  = "None",
        validMassPoints = MassPoints,
    ))


# Reserve column 2
DataGroups.append(DataGroup(
    label        = "res.",
    landsProcess = 2,
    datasetType  = "None",
    validMassPoints = MassPoints,
))


##############################################################################
# Definition of nuisance parameters
#
# Note: Remember to include 'stat.' into the label of nuistances of statistical nature
#
from HiggsAnalysis.HeavyChHiggsToTauNu.datacardtools.InputClasses import Nuisance
ReservedNuisances = []
ReservedNuisances.append(["05", "reserved for leptonic"])
ReservedNuisances.append(["06", "reserved for leptonic"])
ReservedNuisances.append(["08", "reserved for leptonic"])
ReservedNuisances.append(["20", "reserved for leptonic"])
ReservedNuisances.append(["21", "reserved for leptonic"])
ReservedNuisances.append(["23", "reserved for leptonic"])
if OptionPurgeReservedLines:
    ReservedNuisances = []

Nuisances = []

Nuisances.append(Nuisance(
    id            = "01",
    label         = "tau+MET trg tau part",
    distr         = "lnN",
    function      = "ScaleFactor",
    histoDir      = ["ScaleFactorUncertainties"],
    histograms    = ["TriggerScaleFactorAbsUncert_AfterSelection"],
    normalisation = ["TriggerScaleFactorAbsUncertCounts_AfterSelection"],
))

Nuisances.append(Nuisance(
    id            = "01b",
    label         = "tau+MET trg MET part for EWKtau (temp)",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.066
))

Nuisances.append(Nuisance(
    id            = "01c",
    label         = "tau+MET trg scale factor",
    distr         = "lnN",
    function      = "ScaleFactor",
    histoDir      = ["ScaleFactorUncertainties",
                     "ScaleFactorUncertainties"],
    histograms    = ["TriggerScaleFactorAbsUncert_AfterDeltaPhi160",
                     "TriggerScaleFactorAbsUncert_AfterDeltaPhi160"],
    normalisation = ["TriggerScaleFactorAbsUncertCounts_AfterSelection",
                     "TriggerScaleFactorAbsUncertCounts_AfterSelection"],
    #addUncertaintyInQuadrature = 0.10 # MET leg uncertainty
))

Nuisances.append(Nuisance(
    id            = "02",
    label         = "tau+MET trg MET part",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.10
))

Nuisances.append(Nuisance(
    id            = "03",
    label         = "tau-jet ID (no Rtau)",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.06
))

Nuisances.append(Nuisance(
    id            = "04", 
    label         = "tau-jet mis ID (no Rtau)",  
    distr         = "lnN",
    function      = "ScaleFactor",
    histoDir      = ["ScaleFactorUncertainties"],
    histograms    = ["FakeTauAbsUncert_AfterSelection"],
    normalisation = ["FakeTauAbsUncertCounts_AfterSelection"],
    #function      = "Constant",
    #value         = 0.15
))

if OptionIncludeSystematics:
    Nuisances.append(Nuisance(
        id            = "45",
        label         = "TES effect on shape",
        distr         = "shapeQ",
        function      = "Shape",
        counter       = SignalRateCounter,
        histoDir      = ["TESPlus",
                        "TESMinus"],
        histograms    = [SignalShapeHisto,
                        SignalShapeHisto]
    ))
    #Nuisances.append(Nuisance(
        #id            = "45b",
        #label         = "TES effect on shape",
        #distr         = "shapeQ",
        #function      = "Shape",
        #counter       = SignalRateCounter,
        #histoDir      = ["TESPlus",
                        #"TESMinus"],
        #histograms    = [SignalShapeHisto,
                        #SignalShapeHisto]
    #))
    Nuisances.append(Nuisance(
        id            = "46",
        label         = "JES effect on shape",
        distr         = "shapeQ",
        function      = "Shape",
        counter       = SignalRateCounter,
        histoDir      = ["JESPlus",
                        "JESMinus"],
        histograms    = [SignalShapeHisto,
                        SignalShapeHisto]
    ))

    Nuisances.append(Nuisance(
        id            = "47",
        label         = "MET unclustered scale effect on shape",
        distr         = "shapeQ",
        function      = "Shape",
        counter       = SignalRateCounter,
        histoDir      = ["METPlus",
                        "METMinus"],
        histograms    = [SignalShapeHisto,
                        SignalShapeHisto]
    ))
    Nuisances.append(Nuisance(
        id            = "48",
        label         = "Temporary TES for EWKtau",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.066,
    ))
else:
    Nuisances.append(Nuisance(
        id            = "45",
        label         = "NON-EXACT VALUE for JES/JER/MET/Rtau effect on mT shape",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.03,
    ))

    Nuisances.append(Nuisance(
        id            = "46",
        label         = "NON-EXACT VALUE for JES/JER/MET/Rtau effect on mT shape",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.03,
    ))

    Nuisances.append(Nuisance(
        id            = "47",
        label         = "NON-EXACT VALUE for JES/JER/MET/Rtau effect on mT shape",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.01,
    ))

Nuisances.append(Nuisance(
    id            = "09",
    label         = "lepton veto",
    distr         = "lnN",
    function      = "Ratio",
    numerator     = "muon veto",
    denominator   = "trigger scale factor",
    scaling       = 0.02
))

Nuisances.append(Nuisance(
    id            = "10",
    label         = "btagging",
    distr         = "lnN",
    function      = "ScaleFactor",
    histoDir      = ["ScaleFactorUncertainties"],
    histograms    = ["BtagScaleFactorAbsUncert_AfterSelection"],
    normalisation = ["BtagScaleFactorAbsUncertCounts_AfterSelection"]
))

Nuisances.append(Nuisance(
    id            = "11",
    label         = "mistagging",
    distr         = "lnN",
    function      = "ScaleFactor",
    histoDir      = ["ScaleFactorUncertainties"],
    histograms    = ["BtagScaleFactorAbsUncert_AfterSelection"],
    normalisation = ["BtagScaleFactorAbsUncertCounts_AfterSelection"]
))

Nuisances.append(Nuisance(
    id            = "11b",
    label         = "mistagging",
    distr         = "lnN",
    function      = "ScaleFactor",
    histoDir      = ["ScaleFactorUncertainties/"],
    histograms    = ["BtagScaleFactorAbsUncert_AfterSelection"],
    normalisation = ["BtagScaleFactorAbsUncertCounts_AfterSelection"]
))

Nuisances.append(Nuisance(
    id            = "12",
    label         = "QCD stat.",
    distr         = "lnN",
    function      = "QCDFactorised",
    QCDmode       = "statistics",
))

Nuisances.append(Nuisance(
    id            = "13",
    label         = "QCD syst.",
    distr         = "lnN",
    function      = "QCDFactorised",
    QCDmode       = "systematics",
))

Nuisances.append(Nuisance(
    id            = "14",
    label         = "EWK with taus QCD contamination",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.012
))

Nuisances.append(Nuisance(
    id            = "15",
    label         = "EWK with taus W->tau->mu",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.007
))

Nuisances.append(Nuisance(
    id            = "15b",
    label         = "EWK with taus W->tau->mu",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.001
))

Nuisances.append(Nuisance(
    id            = "16",
    label         = "EWK with taus muon selection+ditau+mu trg",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.031
))

Nuisances.append(Nuisance(
    id            = "16b",
    label         = "EWK with taus muon selection",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.001
))

Nuisances.append(Nuisance(
    id            = "17",
    label         = "MC signal stat., HH",
    distr         = "lnN",
    function      = "Counter",
    counter       = SignalRateCounter,
))

Nuisances.append(Nuisance(
    id            = "18",
    label         = "MC signal stat., HW",
    distr         = "lnN",
    function      = "Counter",
    counter       = SignalRateCounter,
))

Nuisances.append(Nuisance(
    id            = "19",
    label         = "EWK with taus stat.",
    distr         = "lnN",
    function      = "Counter",
    counter       = SignalRateCounter,
))

Nuisances.append(Nuisance(
    id            = "22",
    label         = "tt->jjtau MC stat.",
    distr         = "lnN",
    function      = "Constant",
    value         = 0
))

Nuisances.append(Nuisance(
    id            = "24",
    label         = "Z->tautau MC stat.",
    distr         = "lnN",
    function      = "Counter",
    counter       = SignalRateCounter,
))

Nuisances.append(Nuisance(
    id            = "25",
    label         = "W+jets MC stat.",
    distr         = "lnN",
    function      = "Constant",
    value         = 0
))

Nuisances.append(Nuisance(
    id            = "26",
    label         = "Single top MC stat.",
    distr         = "lnN",
    function      = "Constant",
    value         = 0
))

Nuisances.append(Nuisance(
    id            = "27",
    label         = "diboson MC stat.",
    distr         = "lnN",
    function      = "Counter",
    counter       = FakeRateCounter
))

Nuisances.append(Nuisance(
    id            = "28",
    label         = "ttbar cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.096,
    upperValue    = 0.070,
))

Nuisances.append(Nuisance(
    id            = "29",
    label         = "W+jets cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.05
))

Nuisances.append(Nuisance(
    id            = "30",
    label         = "single top cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.08
))

Nuisances.append(Nuisance(
    id            = "31",
    label         = "Z->ll cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.04
))

Nuisances.append(Nuisance(
    id            = "32",
    label         = "diboson cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.04
))

Nuisances.append(Nuisance(
    id            = "33",
    label         = "luminosity",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.022
))

if OptionIncludeSystematics:
    Nuisances.append(Nuisance(
        id            = "34",
        label         = "pileup",
        distr         = "lnN",
        function      = "pileupUncertainty",
        histoDir      = ["", # nominal
                        "PUWeightMinus", # minus
                        "PUWeightPlus"], # up
        counter       = SignalRateCounter
    ))

    Nuisances.append(Nuisance(
        id            = "34b",
        label         = "pileup",
        distr         = "lnN",
        function      = "pileupUncertainty",
        histoDir      = ["", # nominal
                        "PUWeightMinus", # minus
                        "PUWeightPlus"], # up
        counter       = FakeRateCounter
    ))
else:
    Nuisances.append(Nuisance(
        id            = "34",
        label         = "FAKE pileup",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.05
    ))

    Nuisances.append(Nuisance(
        id            = "34b",
        label         = "FAKE pileup",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.05
    ))

Nuisances.append(Nuisance(
    id            = "35",
    label         = "ttbar fake tau MC stat.",
    distr         = "lnN",
    function      = "Counter",
    counter       = FakeRateCounter
))

Nuisances.append(Nuisance(
    id            = "36", 
    label         = "Z->tautau fake tau MC stat.",
    distr         = "lnN",
    function      = "Counter",
    counter       = FakeRateCounter
))

Nuisances.append(Nuisance(
    id            = "37",  
    label         = "W+jets fake tau MC stat.",
    distr         = "lnN",
    function      = "Counter",
    counter       = FakeRateCounter
))

Nuisances.append(Nuisance(
    id            = "38",
    label         = "single top fake tau MC stat.",
    distr         = "lnN",
    function      = "Counter",
    counter       = FakeRateCounter
))

Nuisances.append(Nuisance(
    id            = "39", 
    label         = "diboson fake tau MC stat.",
    distr         = "lnN",
    function      = "Counter",
    counter       = FakeRateCounter
))

Nuisances.append(Nuisance(
    id            = "40",
    label         = "Stat. uncertainty on the shape",
    distr         = "shapeStat",
    function      = "Shape",
    counter       = SignalRateCounter,
    histoDir      = [""],
    histograms    = [SignalShapeHisto],
))

Nuisances.append(Nuisance(
    id            = "40b",
    label         = "Stat. uncertainty on the shape",
    distr         = "shapeStat",
    function      = "QCDFactorised",
    QCDmode       = "shapestat"
))

Nuisances.append(Nuisance(
    id            = "40c",
    label         = "Stat. uncertainty on the shape",
    distr         = "shapeStat",
    function      = "Shape",
    counter       = FakeRateCounter,
    histoDir      = [""],
    histograms    = [FakeShapeHisto],
))

Nuisances.append(Nuisance(
    id            = "41",
    label         = "QCDInv: stat.",
    distr         = "lnN",
    function      = "QCDInverted",
    counter       = "deltaPhiTauMET160 limit"
))

Nuisances.append(Nuisance(
    id            = "42",
    label         = "QCDInv: JES/JER/MET/Rtau effect on normalisation",
    distr         = "lnN",
    function      = "Constant",
    value         = 0.057
))

Nuisances.append(Nuisance(
    id            = "43",
    label         = "QCDInv: MET shape", 
    distr         = "lnN",
    function      = "Constant",
    value         = 0.055
))

Nuisances.append(Nuisance(
    id            = "44",
    label         = "QCDInv: fit", 
    distr         = "lnN",
    function      = "Constant", 
    value         = 0.0043
))

MergeNuisances = []
MergeNuisances.append(["01","01b","01c"])
#MergeNuisances.append(["07","07b","07c"])
MergeNuisances.append(["11","11b"])
MergeNuisances.append(["15","15b"])
MergeNuisances.append(["16","16b"])
MergeNuisances.append(["34","34b"])
MergeNuisances.append(["40","40b","40c"])

# Control plots
from HiggsAnalysis.HeavyChHiggsToTauNu.datacardtools.InputClasses import ControlPlotInput
ControlPlots = []
ControlPlots.append(ControlPlotInput(
    title            = "Njets",
    signalHHid       = [-1],
    signalHWid       = [0],
    QCDid            = [3],
    embeddingId      = EmbeddingIdList,
    EWKfakeId        = EWKFakeIdList,
    signalHistoPath  = "ControlPlots",
    signalHistoName  = "Njets_AfterStandardSelections",
    EWKfakeHistoPath  = "ControlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "Njets_AfterStandardSelections",
    QCDFactNormalisation = "factorisation/AfterJetSelection",
    QCDFactHistoPath = "shape_CtrlLeg1AfterNjets",
    QCDFactHistoName = "CtrlLeg1AfterNjets",
    details          = { "bins": 5,
                         "rangeMin": 3.0,
                         "rangeMax": 8.0,
                         "variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         "xtitle": "Number of selected jets",
                         "ytitle": "Events",
                         "unit": "",
                         "logy": True,
                         "DeltaRatio": 0.5,
                         "ymin": 0.9,
                         "ymax": -1},
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "MET",
    signalHHid       = [-1],
    signalHWid       = [0],
    QCDid            = [3],
    embeddingId      = EmbeddingIdList,
    EWKfakeId        = EWKFakeIdList,
    signalHistoPath  = "ControlPlots",
    signalHistoName  = "MET",
    EWKfakeHistoPath  = "ControlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "MET",
    QCDFactNormalisation = "factorisation/AfterJetSelection",
    QCDFactHistoPath = "shape_CtrlLeg1AfterMET",
    QCDFactHistoName = "CtrlLeg1AfterMET",
    details          = { "bins": 13,
                         "rangeMin": 0.0,
                         "rangeMax": 500.0,
                         #"variableBinSizeLowEdges": [0,20,40,60,80,100,120,140,160,180,200,250,300], # if an empty list is given, then uniform bin width is used
                         "variableBinSizeLowEdges": [0,20,40,60,80,100,120,140,160,180,200,250,300], # if an empty list is given, then uniform bin width is used
                         "xtitle": "E_{T}^{miss}",
                         "ytitle": "Events",
                         "unit": "GeV",
                         "logy": True,
                         "DeltaRatio": 0.5,
                         "ymin": 0.9,
                         "ymax": -1},
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "#tau_{h}+#geq3j", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "BJetSelection",
    signalHHid       = [-1],
    signalHWid       = [0],
    QCDid            = [3],
    embeddingId      = EmbeddingIdList,
    EWKfakeId        = EWKFakeIdList,
    signalHistoPath  = "ControlPlots",
    signalHistoName  = "NBjets",
    EWKfakeHistoPath  = "ControlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "NBjets",
    QCDFactNormalisation = "factorisation/Leg1AfterMET",
    QCDFactHistoPath = "shape_CtrlLeg1AfterNbjets",
    QCDFactHistoName = "CtrlLeg1AfterNbjets",
    details          = { "bins": 5,
                         "rangeMin": 0.0,
                         "rangeMax": 5.0,
                         "variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         "xtitle": "Number of selected b jets",
                         "ytitle": "Events",
                         "unit": "",
                         "logy": True,
                         "DeltaRatio": 0.5,
                         "ymin": 0.9,
                         "ymax": -1},
    blindedRange     = [1.5,10], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "E_{T}^{miss}", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "DeltaPhi",
    signalHHid       = [-1],
    signalHWid       = [0],
    QCDid            = [3],
    embeddingId      = EmbeddingIdList,
    EWKfakeId        = EWKFakeIdList,
    signalHistoPath  = "",
    signalHistoName  = "deltaPhi",
    EWKfakeHistoPath  = "",
    EWKfakeHistoName  = "EWKFakeTausDeltaPhi",
    QCDFactNormalisation = "factorisation/Leg1AfterBTagging",
    QCDFactHistoPath = "shape_CtrlLeg1AfterDeltaPhiTauMET",
    QCDFactHistoName = "CtrlLeg1AfterDeltaPhiTauMET",
    details          = { "bins": 18,
                         "rangeMin": 0.0,
                         "rangeMax": 180.0,
                         "variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         "xtitle": "#Delta#phi(#tau_{h},E_{T}^{miss})",
                         "ytitle": "Events",
                         "unit": "^{o}",
                         "logy": True,
                         "DeltaRatio": 0.5,
                         "ymin": 0.9,
                         "ymax": -1},
    blindedRange     = [-1, 300], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "N_{b jets}", # Leave blank if you don't want to include the item to the selection flow plot
))

#ControlPlots.append(ControlPlotInput(
    #title            = "MaxDeltaPhi",
    #signalHHid       = [-1],
    #signalHWid       = [0],
    #QCDid            = [3],
    #embeddingId      = EmbeddingIdList,
    #EWKfakeId        = EWKFakeIdList,
    #signalHistoPath  = "",
    #signalHistoName  = "maxDeltaPhiJetMet",
    #QCDFactNormalisation = "factorisation/Leg1AfterDeltaPhiTauMET",
    #QCDFactHistoPath = "shape_CtrlLeg1AfterMaxDeltaPhiJetMET",
    #QCDFactHistoName = "CtrlLeg1AfterMaxDeltaPhiJetMET",
    #details          = { "bins": 18,
                         #"rangeMin": 0.0,
                         #"rangeMax": 180.0,
                         #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         #"xtitle": "max(#Delta#phi(jet,E_{T}^{miss})",
                         #"ytitle": "Events",
                         #"unit": "^{o}",
                         #"logy": True,
                         #"DeltaRatio": 0.5,
                         #"ymin": 0.9,
                         #"ymax": -1},
    #blindedRange     = [-1, 300], # specify range min,max if blinding applies to this control plot
    #evaluationRange  = [], # specify range to be evaluated and saved into a file
    #flowPlotCaption  = "#Delta#phi(#tau_{h},E_{T}^{miss})", # Leave blank if you don't want to include the item to the selection flow plot
#))

#ControlPlots.append(ControlPlotInput(
    #title            = "WMass",
    #signalHHid       = [-1],
    #signalHWid       = [0],
    #QCDid            = [3],
    #embeddingId      = EmbeddingIdList,
    #EWKfakeId        = EWKFakeIdList,
    #signalHistoPath  = "TopChiSelection",
    #signalHistoName  = "WMass",
    #QCDFactNormalisation = "factorisation/Leg1AfterDeltaPhiTauMET",
    #QCDFactHistoPath = "shape_CtrlLeg1AfterTopMass",
    #QCDFactHistoName = "CtrlLeg1AfterTopMass",
    #details          = { "bins": 20,
                         #"rangeMin": 0.0,
                         #"rangeMax": 200.0,
                         #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         #"xtitle": "m_{jj}",
                         #"ytitle": "Events",
                         #"unit": "GeV/c^{2}",
                         #"logy": True,
                         #"DeltaRatio": 0.5,
                         #"ymin": 0.9,
                         #"ymax": -1},
    #blindedRange     = [-1, 400], # specify range min,max if blinding applies to this control plot
    #evaluationRange  = [], # specify range to be evaluated and saved into a file
    #flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
#))

#ControlPlots.append(ControlPlotInput(
    #title            = "TopMass",
    #signalHHid       = [-1],
    #signalHWid       = [0],
    #QCDid            = [3],
    #embeddingId      = EmbeddingIdList,
    #EWKfakeId        = EWKFakeIdList,
    #signalHistoPath  = "TopChiSelection",
    #signalHistoName  = "TopMass",
    #QCDFactNormalisation = "factorisation/Leg1AfterDeltaPhiTauMET",
    #QCDFactHistoPath = "shape_CtrlLeg1AfterTopMass",
    #QCDFactHistoName = "CtrlLeg1AfterTopMass",
    #details          = { "bins": 20,
                         #"rangeMin": 0.0,
                         #"rangeMax": 400.0,
                         #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         #"xtitle": "m_{bjj}",
                         #"ytitle": "Events",
                         #"unit": "GeV/c^{2}",
                         #"logy": True,
                         #"DeltaRatio": 0.5,
                         #"ymin": 0.9,
                         #"ymax": -1},
    #blindedRange     = [-1, 400], # specify range min,max if blinding applies to this control plot
    #evaluationRange  = [], # specify range to be evaluated and saved into a file
    #flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
#))

ControlPlots.append(ControlPlotInput(
    title            = "TransverseMass",
    signalHHid       = [-1],
    signalHWid       = [0],
    QCDid            = [3],
    embeddingId      = EmbeddingIdList,
    EWKfakeId        = EWKFakeIdList,
    signalHistoPath  = "",
    signalHistoName  = "transverseMass",
    EWKfakeHistoPath  = "",
    EWKfakeHistoName  = "EWKFakeTausTransverseMass",
    QCDFactNormalisation = "factorisation/Leg1AfterTopSelection",
    QCDFactHistoPath = "shape_MtShapesAfterFullMETLeg",
    QCDFactHistoName = "MtShapesAfterFullMETLeg",
    details          = { "bins": 13,
                         "rangeMin": 0.0,
                         "rangeMax": 400.0,
                         "variableBinSizeLowEdges": [0,20,40,60,80,100,120,140,160,180,200,250,300], # if an empty list is given, then uniform bin width is used
                         "xtitle": "mT(#tau_{h},E_{T}^{miss})",
                         "ytitle": "Events",
                         "unit": "GeV/c^{2}",
                         "logy": False,
                         "DeltaRatio": 0.5,
                         "ymin": 0.9,
                         "ymax": -1},
    blindedRange     = [-1, 1000], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [60, 180], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "#Delta#phi(#tau_{h},E_{T}^{miss})", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "FullMass",
    signalHHid       = [-1],
    signalHWid       = [0],
    QCDid            = [3],
    embeddingId      = EmbeddingIdList,
    EWKfakeId        = EWKFakeIdList,
    signalHistoPath  = "FullHiggsMass",
    signalHistoName  = "HiggsMass",
    EWKfakeHistoPath  = "",
    EWKfakeHistoName  = "EWKFakeTausFullMass",
    QCDFactNormalisation = "factorisation/Leg1AfterTopSelection",
    QCDFactHistoPath = "shape_FullMassShapesAfterFullMETLeg",
    QCDFactHistoName = "FullMassShapesAfterFullMETLeg",
    details          = { "bins": 13,
                         "rangeMin": 0.0,
                         "rangeMax": 500.0,
                         "variableBinSizeLowEdges": [0,20,40,60,80,100,120,140,160,180,200,250,300], # if an empty list is given, then uniform bin width is used
                         "xtitle": "m(#tau_{h},E^{miss})",
                         "ytitle": "Events",
                         "unit": "GeV/c^{2}",
                         "logy": False,
                         "DeltaRatio": 0.5,
                         "ymin": 0.9,
                         "ymax": -1},
    blindedRange     = [-1, 1000], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [80, 180], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))


