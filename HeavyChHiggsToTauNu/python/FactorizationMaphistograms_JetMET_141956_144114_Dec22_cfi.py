import FWCore.ParameterSet.Config as cms

tauIDFactorizationCoefficients = cms.untracked.PSet(
factorizationSourceName = cms.untracked.string('histograms_JetMET_141956_144114_Dec22'),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.00666667, 0.00380228, 0.003861, 1e-99, 0.00956938, 0.00520833, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.00666667, 0.00380228, 0.003861, 0.00380228, 0.00676657, 0.00520833, 0.00666667, 0.00485437, 0.00549451, 0.00595238
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0114943, 1e-99, 1e-99, 0.0103093, 1e-99, 1e-99, 1e-99, 0.00704225, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.00854701, 0.0149254, 1e-99, 1e-99, 0.027027, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.0217391, 0.0208333, 0.0114943, 0.0140845, 0.0149254, 0.0103093, 0.0126582, 0.0243902, 0.015873, 0.00704225, 0.016129, 0.00884956, 0.00847458, 0.0060241, 0.0102041, 0.0185185, 0.00694444, 0.0227273, 0.0136986, 0.0178571, 0.00854701, 0.0149254, 0.0294118, 0.00847458, 0.027027, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.0625, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 0.0714286, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0526316, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.05, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.166667, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.166667, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.5, 0.2, 0.142857, 0.2, 0.2, 0.25, 0.25, 0.125, 0.333333, 0.333333,
  1e-99, 1e-99, 0.5, 0.0909091, 0.1, 0.2, 0.166667, 0.2, 1, 0.25, 0.5, 0.5,
  1e-99, 1e-99, 0.2, 0.1, 0.0833333, 0.0909091, 0.0625, 0.111111, 0.1, 0.111111, 0.25, 1,
  1e-99, 1e-99, 1, 0.111111, 0.111111, 0.1, 0.166667, 0.166667, 0.142857, 0.125, 0.142857, 0.125,
  1e-99, 1e-99, 0.142857, 0.0769231, 0.5, 0.0833333, 0.2, 0.2, 0.142857, 0.125, 0.333333, 0.2,
  1e-99, 1e-99, 0.1, 0.0714286, 0.0714286, 0.0833333, 0.111111, 0.1, 0.25, 0.166667, 0.125, 0.1,
  1e-99, 1e-99, 0.1, 0.111111, 0.1, 0.0714286, 0.0909091, 0.142857, 0.25, 0.142857, 0.25, 0.333333,
  1e-99, 1e-99, 0.2, 0.142857, 0.166667, 0.142857, 1, 0.5, 0.333333, 0.333333, 0.25, 0.333333,
  1e-99, 1e-99, 0.25, 0.0714286, 0.0833333, 0.166667, 0.142857, 0.2, 0.333333, 0.333333, 0.125, 1,
  1e-99, 1e-99, 0.142857, 0.0526316, 0.0769231, 0.0714286, 0.0625, 0.0714286, 0.0909091, 0.0625, 0.0666667, 0.0588235,
  1e-99, 1e-99, 1, 0.2, 0.1, 0.0909091, 0.125, 0.333333, 0.333333, 0.142857, 0.125, 0.166667,
  1e-99, 1e-99, 0.111111, 0.0769231, 0.125, 0.1, 0.0909091, 0.0769231, 0.0714286, 0.0714286, 0.0909091, 0.1,
  1e-99, 1e-99, 0.142857, 0.0833333, 0.0714286, 0.0555556, 0.0909091, 0.1, 0.166667, 0.0666667, 0.0625, 0.111111,
  1e-99, 1e-99, 1, 0.0666667, 0.0526316, 0.0588235, 0.0625, 0.0555556, 0.05, 0.0384615, 0.0666667, 0.0526316,
  1e-99, 1e-99, 0.125, 0.0588235, 0.0769231, 0.142857, 0.166667, 0.111111, 0.166667, 0.1, 0.1, 0.0833333,
  1e-99, 1e-99, 0.25, 0.166667, 1, 0.142857, 0.166667, 0.111111, 0.333333, 0.2, 0.25, 0.111111,
  1e-99, 1e-99, 0.2, 0.0588235, 0.0625, 0.0526316, 0.0714286, 0.111111, 0.111111, 0.0666667, 0.05, 0.05,
  1e-99, 1e-99, 0.333333, 0.25, 0.142857, 0.166667, 0.2, 0.5, 0.2, 0.333333, 0.5, 0.142857,
  1e-99, 1e-99, 0.2, 0.1, 0.0833333, 0.125, 0.25, 0.111111, 0.142857, 0.166667, 0.142857, 0.2,
  1e-99, 1e-99, 0.142857, 0.2, 0.0769231, 0.125, 0.25, 0.333333, 0.25, 0.166667, 0.25, 0.5,
  1e-99, 1e-99, 0.05, 0.0526316, 0.0769231, 0.0555556, 0.0714286, 0.111111, 0.2, 0.142857, 0.166667, 0.166667,
  1e-99, 1e-99, 0.125, 0.125, 0.0769231, 0.125, 0.166667, 0.166667, 0.2, 0.2, 0.2, 0.333333,
  1e-99, 1e-99, 0.2, 0.142857, 0.333333, 0.2, 0.333333, 0.5, 1, 0.25, 0.25, 1e-99,
  1e-99, 1e-99, 0.0769231, 0.0833333, 0.0625, 0.0526316, 0.0769231, 0.0526316, 0.2, 0.125, 0.142857, 0.166667,
  1e-99, 1e-99, 1, 0.5, 0.166667, 0.166667, 0.166667, 0.25, 0.333333, 0.333333, 0.2, 1,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0133333, 0.00760456, 1e-99, 1e-99, 0.0143541, 0.00520833, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.00942809, 0.00537724, 0.003861, 0.00380228, 0.00828732, 0.00520833, 0.00666667, 0.00485437, 0.00549451, 0.00595238
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0114943, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.00704225, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.00854701, 0.0447761, 1e-99, 0.00847458, 0.027027, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.0217391, 0.0208333, 0.0114943, 0.0140845, 0.0149254, 0.0103093, 0.0126582, 0.0243902, 0.015873, 0.00704225, 0.016129, 0.00884956, 0.00847458, 0.0060241, 0.0102041, 0.0185185, 0.00694444, 0.0227273, 0.0136986, 0.0178571, 0.00854701, 0.0258515, 0.0294118, 0.00847458, 0.027027, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.0625, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.0526316, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.05, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.125, 1e-99, 1e-99, 0.166667, 0.166667, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0769231, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.166667, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.5, 0.2, 0.142857, 0.2, 0.2, 0.25, 0.25, 0.125, 0.333333, 0.333333,
  1e-99, 1e-99, 0.5, 0.0909091, 0.1, 0.2, 0.166667, 0.2, 1, 0.25, 0.5, 0.5,
  1e-99, 1e-99, 0.2, 0.1, 0.0833333, 0.0909091, 0.0625, 0.111111, 0.1, 0.111111, 0.25, 1,
  1e-99, 1e-99, 1, 0.111111, 0.111111, 0.1, 0.166667, 0.166667, 0.142857, 0.125, 0.142857, 0.125,
  1e-99, 1e-99, 0.142857, 0.0769231, 0.5, 0.0833333, 0.2, 0.2, 0.142857, 0.125, 0.333333, 0.2,
  1e-99, 1e-99, 0.1, 0.0714286, 0.0714286, 0.0833333, 0.111111, 0.1, 0.25, 0.166667, 0.125, 0.1,
  1e-99, 1e-99, 0.1, 0.111111, 0.1, 0.0714286, 0.0909091, 0.142857, 0.25, 0.142857, 0.25, 0.333333,
  1e-99, 1e-99, 0.2, 0.142857, 0.166667, 0.142857, 1, 0.5, 0.333333, 0.333333, 0.25, 0.333333,
  1e-99, 1e-99, 0.25, 0.0714286, 0.0833333, 0.166667, 0.142857, 0.2, 0.333333, 0.333333, 0.125, 1,
  1e-99, 1e-99, 0.142857, 0.0526316, 0.0769231, 0.0714286, 0.0625, 0.0714286, 0.0909091, 0.0625, 0.0666667, 0.0588235,
  1e-99, 1e-99, 1, 0.2, 0.1, 0.0909091, 0.125, 0.333333, 0.333333, 0.142857, 0.125, 0.166667,
  1e-99, 1e-99, 0.111111, 0.0769231, 0.125, 0.1, 0.0909091, 0.0769231, 0.0714286, 0.0714286, 0.0909091, 0.1,
  1e-99, 1e-99, 0.142857, 0.0833333, 0.0714286, 0.0555556, 0.0909091, 0.1, 0.166667, 0.0666667, 0.0625, 0.111111,
  1e-99, 1e-99, 1, 0.0666667, 0.0526316, 0.0588235, 0.0625, 0.0555556, 0.05, 0.0384615, 0.0666667, 0.0526316,
  1e-99, 1e-99, 0.125, 0.0588235, 0.0769231, 0.142857, 0.166667, 0.111111, 0.166667, 0.1, 0.1, 0.0833333,
  1e-99, 1e-99, 0.25, 0.166667, 1, 0.142857, 0.166667, 0.111111, 0.333333, 0.2, 0.25, 0.111111,
  1e-99, 1e-99, 0.2, 0.0588235, 0.0625, 0.0526316, 0.0714286, 0.111111, 0.111111, 0.0666667, 0.05, 0.05,
  1e-99, 1e-99, 0.333333, 0.25, 0.142857, 0.166667, 0.2, 0.5, 0.2, 0.333333, 0.5, 0.142857,
  1e-99, 1e-99, 0.2, 0.1, 0.0833333, 0.125, 0.25, 0.111111, 0.142857, 0.166667, 0.142857, 0.2,
  1e-99, 1e-99, 0.142857, 0.2, 0.0769231, 0.125, 0.25, 0.333333, 0.25, 0.166667, 0.25, 0.5,
  1e-99, 1e-99, 0.05, 0.0526316, 0.0769231, 0.0555556, 0.0714286, 0.111111, 0.2, 0.142857, 0.166667, 0.166667,
  1e-99, 1e-99, 0.125, 0.125, 0.0769231, 0.125, 0.166667, 0.166667, 0.2, 0.2, 0.2, 0.333333,
  1e-99, 1e-99, 0.2, 0.142857, 0.333333, 0.2, 0.333333, 0.5, 1, 0.25, 0.25, 1e-99,
  1e-99, 1e-99, 0.0769231, 0.0833333, 0.0625, 0.0526316, 0.0769231, 0.0526316, 0.2, 0.125, 0.142857, 0.166667,
  1e-99, 1e-99, 1, 0.5, 0.166667, 0.166667, 0.166667, 0.25, 0.333333, 0.333333, 0.2, 1,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.00249377, 0.003125, 1e-99, 1e-99, 0.00970874, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0025974, 0.00249377, 0.003125, 0.00458716, 0.00666667, 0.00970874, 0.0128205, 0.00952381, 0.0178571, 0.025
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.0116279, 1e-99, 1e-99, 1e-99, 0.00757576, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.0169492, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.025, 0.0222222, 0.0126582, 0.015873, 0.0175439, 0.0116279, 0.015873, 0.0243902, 0.0185185, 0.00757576, 0.0196078, 0.00900901, 0.00952381, 0.00574713, 0.0123457, 0.0204082, 0.00775194, 0.0243902, 0.0135135, 0.02, 0.00990099, 0.0169492, 0.0357143, 0.00925926, 0.0285714, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 0.05, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.03125, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.166667, 0.142857, 0.111111, 0.2, 0.333333, 0.2, 0.5, 0.5, 1e-99, 1,
  1e-99, 1e-99, 0.0909091, 0.0833333, 0.0833333, 0.5, 0.2, 0.5, 1e-99, 1e-99, 1e-99, 1,
  1e-99, 1e-99, 0.0909091, 0.0714286, 0.0714286, 0.0714286, 0.142857, 0.111111, 0.2, 0.25, 1, 1e-99,
  1e-99, 1e-99, 0.125, 0.0769231, 0.0909091, 0.111111, 0.2, 0.333333, 1, 0.166667, 0.25, 0.333333,
  1e-99, 1e-99, 0.0833333, 0.0769231, 0.1, 0.25, 0.2, 0.25, 1e-99, 0.25, 0.25, 1,
  1e-99, 1e-99, 0.0588235, 0.0555556, 0.05, 0.1, 0.166667, 0.2, 0.5, 0.333333, 0.25, 1,
  1e-99, 1e-99, 0.05, 0.0526316, 0.125, 0.2, 0.25, 1, 0.25, 0.5, 1e-99, 1e-99,
  1e-99, 1e-99, 0.142857, 0.0909091, 0.125, 0.2, 1e-99, 0.5, 0.25, 1, 1e-99, 0.333333,
  1e-99, 1e-99, 0.0588235, 0.0769231, 0.0833333, 0.25, 0.25, 1e-99, 1, 1, 0.5, 1e-99,
  1e-99, 1e-99, 0.0434783, 0.03125, 0.0454545, 0.0909091, 0.0833333, 0.1, 0.2, 0.125, 0.166667, 0.333333,
  1e-99, 1e-99, 0.142857, 0.0909091, 0.0714286, 0.25, 1, 0.25, 0.5, 0.2, 1e-99, 0.333333,
  1e-99, 1e-99, 0.0434783, 0.047619, 0.0714286, 0.0555556, 0.0833333, 0.166667, 0.25, 0.166667, 0.25, 0.333333,
  1e-99, 1e-99, 0.047619, 0.0555556, 0.047619, 0.0909091, 0.0714286, 0.5, 0.142857, 0.111111, 0.5, 1e-99,
  1e-99, 1e-99, 0.0357143, 0.0277778, 0.0357143, 0.0344828, 0.0833333, 0.125, 0.0833333, 0.125, 0.142857, 0.166667,
  1e-99, 1e-99, 0.0526316, 0.0526316, 0.0833333, 0.125, 0.125, 0.25, 1e-99, 0.125, 0.333333, 1e-99,
  1e-99, 1e-99, 0.1, 0.125, 0.125, 0.142857, 0.2, 0.5, 0.25, 0.333333, 1, 1,
  1e-99, 1e-99, 0.04, 0.03125, 0.0625, 0.1, 0.0714286, 0.125, 0.166667, 0.142857, 0.2, 0.166667,
  1e-99, 1e-99, 0.142857, 0.111111, 0.333333, 0.166667, 0.25, 0.333333, 0.333333, 0.25, 1e-99, 0.5,
  1e-99, 1e-99, 0.05, 0.0833333, 0.0625, 0.142857, 0.25, 0.25, 0.5, 0.25, 0.333333, 0.5,
  1e-99, 1e-99, 0.0833333, 0.0666667, 0.142857, 0.2, 0.25, 0.333333, 0.5, 1e-99, 1, 1,
  1e-99, 1e-99, 0.0333333, 0.047619, 0.0555556, 0.0909091, 0.5, 0.25, 0.5, 0.111111, 0.25, 1e-99,
  1e-99, 1e-99, 0.0769231, 0.0909091, 0.1, 0.125, 0.2, 0.333333, 0.5, 0.25, 0.5, 1,
  1e-99, 1e-99, 0.166667, 0.25, 0.111111, 0.25, 1, 1, 1e-99, 0.333333, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0384615, 0.0434783, 0.0714286, 0.0555556, 0.125, 0.125, 0.25, 0.5, 0.333333, 0.5,
  1e-99, 1e-99, 0.166667, 0.111111, 0.25, 0.333333, 0.2, 0.5, 0.25, 0.5, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0149254, 0.0070922, 0.00609756, 1e-99, 0.00574713, 0.0133333, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0149254, 0.0070922, 0.00609756, 0.00617284, 0.00574713, 0.00942809, 0.00699301, 0.00425532, 0.00396825, 0.0034965
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0126582, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.00826446, 1e-99, 1e-99, 1e-99, 0.00621118, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.016129, 1e-99, 0.00934579, 0.0357143, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.0217391, 0.0232558, 0.0126582, 0.0163934, 0.015625, 0.0144928, 0.0144928, 0.0243902, 0.0172414, 0.00826446, 0.0263158, 0.00917431, 0.00970874, 0.00621118, 0.0136986, 0.0208333, 0.00746269, 0.0227273, 0.0175439, 0.0163934, 0.0149254, 0.016129, 0.0322581, 0.00934579, 0.0357143, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.0909091, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.111111, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 0.0769231, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.111111, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.111111, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.166667, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.5, 0.333333, 0.166667, 0.333333, 0.333333, 0.25, 0.25, 0.111111, 0.166667, 0.166667,
  1e-99, 1e-99, 0.5, 0.25, 0.111111, 0.333333, 1e-99, 0.25, 0.25, 0.142857, 0.5, 0.125,
  1e-99, 1e-99, 1, 0.125, 0.142857, 0.25, 0.0909091, 0.1, 0.142857, 0.111111, 0.0833333, 0.1,
  1e-99, 1e-99, 1, 0.166667, 0.25, 0.5, 0.142857, 0.5, 0.125, 0.125, 0.0769231, 0.1,
  1e-99, 1e-99, 0.5, 0.166667, 0.333333, 0.25, 0.0909091, 0.2, 0.142857, 0.0666667, 0.2, 0.166667,
  1e-99, 1e-99, 0.333333, 0.111111, 0.25, 0.111111, 0.1, 0.25, 0.25, 0.2, 0.0833333, 0.111111,
  1e-99, 1e-99, 0.25, 0.142857, 0.125, 0.125, 0.142857, 0.142857, 0.25, 0.0833333, 0.2, 0.142857,
  1e-99, 1e-99, 1e-99, 0.2, 0.142857, 0.25, 0.25, 0.25, 1e-99, 0.333333, 0.2, 0.111111,
  1e-99, 1e-99, 1e-99, 0.1, 0.166667, 0.166667, 0.125, 0.166667, 0.25, 0.25, 0.142857, 0.142857,
  1e-99, 1e-99, 0.2, 0.111111, 0.2, 0.111111, 0.0666667, 0.125, 0.0833333, 0.0714286, 0.047619, 0.0434783,
  1e-99, 1e-99, 1e-99, 1e-99, 0.333333, 1e-99, 0.333333, 0.2, 1e-99, 0.111111, 0.0909091, 0.142857,
  1e-99, 1e-99, 0.333333, 0.125, 0.111111, 0.1, 0.125, 0.111111, 0.0909091, 0.0769231, 0.05, 0.0555556,
  1e-99, 1e-99, 0.25, 0.142857, 0.166667, 0.0666667, 0.142857, 0.142857, 0.1, 0.0909091, 0.05, 0.0625,
  1e-99, 1e-99, 0.333333, 0.166667, 0.0769231, 0.0769231, 0.0833333, 0.0588235, 0.1, 0.0434783, 0.0333333, 0.0294118,
  1e-99, 1e-99, 0.333333, 0.333333, 0.0714286, 0.125, 0.333333, 0.25, 0.166667, 0.0909091, 0.111111, 0.0833333,
  1e-99, 1e-99, 0.2, 0.5, 1, 1, 0.5, 0.142857, 0.166667, 0.25, 0.125, 0.0833333,
  1e-99, 1e-99, 0.5, 0.0833333, 0.111111, 0.0769231, 0.0769231, 0.166667, 0.0909091, 0.0526316, 0.0526316, 0.0333333,
  1e-99, 1e-99, 1e-99, 0.25, 0.333333, 0.2, 0.333333, 0.333333, 0.333333, 0.166667, 0.25, 0.0769231,
  1e-99, 1e-99, 1e-99, 0.333333, 0.125, 0.166667, 0.142857, 0.142857, 0.166667, 0.2, 0.125, 0.142857,
  1e-99, 1e-99, 0.2, 0.166667, 0.142857, 0.125, 0.142857, 0.2, 0.25, 0.2, 0.142857, 0.142857,
  1e-99, 1e-99, 0.166667, 0.25, 0.0833333, 0.125, 0.111111, 0.333333, 0.25, 0.142857, 0.111111, 0.2,
  1e-99, 1e-99, 0.2, 0.166667, 0.166667, 0.1, 0.25, 0.111111, 0.25, 0.1, 0.333333, 0.2,
  1e-99, 1e-99, 1, 0.25, 1, 0.333333, 0.166667, 1, 0.333333, 0.333333, 0.333333, 0.166667,
  1e-99, 1e-99, 0.111111, 0.125, 0.111111, 0.111111, 0.0909091, 0.142857, 0.142857, 0.0454545, 0.1, 0.0666667,
  1e-99, 1e-99, 1, 1, 0.25, 1, 0.333333, 0.166667, 0.25, 1, 0.333333, 0.25,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.00243902, 0.00337838, 1e-99, 1e-99, 0.0119048, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.00254453, 0.00243902, 0.00337838, 0.00456621, 0.00746269, 0.0119048, 0.015873, 0.0116279, 0.025, 0.0322581
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.0117647, 1e-99, 1e-99, 1e-99, 0.00757576, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.0175439, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.0263158, 0.0238095, 0.0128205, 0.015873, 0.0196078, 0.0117647, 0.0149254, 0.0277778, 0.02, 0.00757576, 0.0232558, 0.0104167, 0.00943396, 0.00641026, 0.0126582, 0.0208333, 0.00854701, 0.0243902, 0.016129, 0.0208333, 0.0105263, 0.0175439, 0.0344828, 0.00925926, 0.0344828, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 0.0555556, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.037037, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.142857, 0.111111, 0.166667, 0.2, 0.5, 0.25, 1, 0.333333, 1e-99, 1,
  1e-99, 1e-99, 0.142857, 0.0714286, 0.0833333, 0.5, 0.5, 0.5, 1e-99, 1, 1, 1,
  1e-99, 1e-99, 0.0714286, 0.111111, 0.0526316, 0.0714286, 0.166667, 0.2, 0.2, 0.25, 0.5, 1e-99,
  1e-99, 1e-99, 0.1, 0.0909091, 0.0625, 0.1, 0.2, 1, 1, 0.25, 0.333333, 0.5,
  1e-99, 1e-99, 0.125, 0.0714286, 0.1, 0.25, 0.25, 0.25, 1e-99, 0.333333, 0.5, 0.5,
  1e-99, 1e-99, 0.047619, 0.0454545, 0.0555556, 0.111111, 0.25, 0.333333, 1, 0.25, 0.5, 1,
  1e-99, 1e-99, 0.0526316, 0.0555556, 0.0909091, 0.2, 0.166667, 0.5, 0.2, 1, 1e-99, 1e-99,
  1e-99, 1e-99, 0.125, 0.1, 0.25, 0.142857, 1e-99, 0.5, 0.333333, 1, 1e-99, 1,
  1e-99, 1e-99, 0.0555556, 0.0833333, 0.1, 0.333333, 0.333333, 1, 1, 1e-99, 0.5, 1e-99,
  1e-99, 1e-99, 0.030303, 0.037037, 0.0526316, 0.047619, 0.1, 0.0833333, 0.333333, 0.25, 1, 0.5,
  1e-99, 1e-99, 0.2, 0.0833333, 0.111111, 0.2, 0.333333, 1, 1e-99, 0.2, 1, 0.5,
  1e-99, 1e-99, 0.0416667, 0.05, 0.0625, 0.0909091, 0.111111, 0.25, 0.5, 0.2, 0.333333, 0.5,
  1e-99, 1e-99, 0.0434783, 0.05, 0.0555556, 0.0909091, 0.0769231, 0.333333, 0.166667, 0.1, 0.5, 1e-99,
  1e-99, 1e-99, 0.0322581, 0.0263158, 0.0384615, 0.0526316, 0.0833333, 0.25, 0.1, 0.142857, 0.333333, 0.166667,
  1e-99, 1e-99, 0.05, 0.047619, 0.0833333, 0.142857, 0.166667, 0.333333, 1, 0.25, 0.2, 1e-99,
  1e-99, 1e-99, 0.0909091, 0.0833333, 0.2, 0.2, 0.166667, 0.5, 0.5, 1, 0.333333, 1,
  1e-99, 1e-99, 0.0384615, 0.027027, 0.0769231, 0.0909091, 0.125, 0.25, 0.333333, 0.166667, 0.25, 0.2,
  1e-99, 1e-99, 0.125, 0.142857, 0.142857, 0.25, 0.2, 0.2, 1, 1, 1e-99, 0.333333,
  1e-99, 1e-99, 0.0833333, 0.0769231, 0.0909091, 0.0909091, 0.166667, 0.5, 0.333333, 0.25, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0555556, 0.0909091, 0.142857, 0.2, 0.333333, 0.5, 1, 1, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0526316, 0.04, 0.0666667, 0.0666667, 0.333333, 0.25, 0.25, 0.125, 0.5, 1e-99,
  1e-99, 1e-99, 0.0909091, 0.0833333, 0.111111, 0.0909091, 0.25, 0.333333, 1, 0.25, 0.5, 1e-99,
  1e-99, 1e-99, 0.142857, 0.166667, 0.2, 0.2, 1, 0.5, 1, 0.5, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0357143, 0.04, 0.0714286, 0.0625, 0.111111, 0.142857, 0.25, 1, 0.5, 0.5,
  1e-99, 1e-99, 0.2, 0.2, 0.25, 0.333333, 0.25, 0.5, 0.25, 0.5, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

)
