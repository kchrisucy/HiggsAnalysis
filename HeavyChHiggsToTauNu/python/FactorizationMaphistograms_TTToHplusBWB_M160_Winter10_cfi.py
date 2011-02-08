import FWCore.ParameterSet.Config as cms

tauIDFactorizationCoefficients = cms.untracked.PSet(
factorizationSourceName = cms.untracked.string('histograms_TTToHplusBWB_M160_Winter10'),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.326141, 0.238667, 0.200883, 0.20618, 0.179048, 0.17138, 0.149721, 0.145679, 0.134843, 0.0917823
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.01418, 0.0102992, 0.00941754, 0.0104807, 0.0106621, 0.0121184, 0.0129339, 0.0109499, 0.0115204, 0.00989714
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  0.285714, 0.12, 0.24, 0.280992, 0.243243, 0.229111, 0.184466, 0.172566, 0.1622, 0.200346, 0.194553, 0.197748, 0.200676, 0.217027, 0.215076, 0.172897, 0.172414, 0.175325, 0.188119, 0.15625, 0.17931, 0.170068, 0.290909, 0.22449, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.142857, 0.069282, 0.0565685, 0.0481897, 0.0362606, 0.0248505, 0.0244331, 0.0276327, 0.0151253, 0.0107552, 0.0158851, 0.0117967, 0.0116444, 0.00978433, 0.0121958, 0.0164107, 0.0109044, 0.0238587, 0.0215787, 0.0247053, 0.0248659, 0.0340136, 0.0727273, 0.0478614, 0.0588235, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.4, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99, 1, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.142857, 0.2, 1e-99, 1e-99, 1e-99, 1, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.25, 0.4, 0.0909091, 0.222222, 0.142857, 0.5, 1, 0.166667, 1e-99, 1e-99,
  1e-99, 1e-99, 0.380952, 0.346154, 0.210526, 0.333333, 0.384615, 1e-99, 0.25, 0.125, 0.142857, 1e-99,
  1e-99, 1e-99, 0.37037, 0.325, 0.181818, 0.16129, 0.375, 0.105263, 0.125, 0.25, 0.2, 1e-99,
  1e-99, 1e-99, 0.466667, 0.346154, 0.296875, 0.204545, 0.0810811, 0.166667, 0.107143, 0.09375, 0.125, 0.111111,
  1e-99, 1e-99, 0.342857, 0.264151, 0.166667, 0.210526, 0.0740741, 0.148148, 0.105263, 0.173913, 0.181818, 1e-99,
  1e-99, 1e-99, 0.409091, 0.296296, 0.107143, 0.1875, 0.142857, 0.0714286, 0.294118, 0.111111, 1e-99, 0.0526316,
  1e-99, 1e-99, 0.329114, 0.268519, 0.144231, 0.2, 0.123288, 0.0816327, 0.116279, 0.0181818, 0.135593, 0.0338983,
  1e-99, 1e-99, 0.324324, 0.213333, 0.230469, 0.209091, 0.177778, 0.177632, 0.145631, 0.175758, 0.15748, 0.092437,
  1e-99, 1e-99, 0.360465, 0.238938, 0.168, 0.193878, 0.147059, 0.148148, 0.195652, 0.103896, 0.148148, 0.18,
  1e-99, 1e-99, 0.280303, 0.175879, 0.225225, 0.226519, 0.187845, 0.220339, 0.161616, 0.166667, 0.148936, 0.114943,
  1e-99, 1e-99, 0.305195, 0.280335, 0.165829, 0.176166, 0.24183, 0.142857, 0.101852, 0.198529, 0.168317, 0.0897436,
  1e-99, 1e-99, 0.3125, 0.241096, 0.214085, 0.223729, 0.192029, 0.203209, 0.245763, 0.157303, 0.16129, 0.142857,
  1e-99, 1e-99, 0.432927, 0.191388, 0.219298, 0.240642, 0.19186, 0.196721, 0.147727, 0.168224, 0.115385, 0.0879121,
  1e-99, 1e-99, 0.288136, 0.194444, 0.19, 0.220779, 0.177419, 0.230769, 0.075, 0.0727273, 0.109091, 0.0851064,
  1e-99, 1e-99, 0.297468, 0.243119, 0.190909, 0.188235, 0.125828, 0.125, 0.11236, 0.125984, 0.0934579, 0.0714286,
  1e-99, 1e-99, 0.21875, 0.230769, 0.192308, 0.255814, 0.0909091, 0.0869565, 0.117647, 0.130435, 0.166667, 0.136364,
  1e-99, 1e-99, 0.285714, 0.262295, 0.137931, 0.175, 0.261905, 0.275862, 0.115385, 0.166667, 0.0714286, 0.0833333,
  1e-99, 1e-99, 0.16, 0.285714, 0.205882, 0.075, 0.206897, 0.130435, 0.111111, 0.1, 0.125, 0.0588235,
  1e-99, 1e-99, 0.340426, 0.2, 0.148148, 0.176471, 0.214286, 0.190476, 1e-99, 0.0588235, 0.0769231, 0.0588235,
  1e-99, 1e-99, 0.285714, 0.0555556, 0.28125, 0.214286, 0.125, 1e-99, 0.142857, 1e-99, 1e-99, 0.111111,
  1e-99, 1e-99, 0.666667, 0.428571, 0.2, 0.333333, 1e-99, 1e-99, 0.5, 1e-99, 0.2, 1e-99,
  1e-99, 1e-99, 0.294118, 0.333333, 0.384615, 0.133333, 1e-99, 0.6, 0.142857, 0.285714, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.282843, 0.333333, 1, 1e-99, 1e-99, 1, 1, 1, 0.5,
  1e-99, 1e-99, 0.2, 0.142857, 0.2, 0.5, 1, 0.5, 1, 1e-99, 0.5, 1e-99,
  1e-99, 1e-99, 0.125, 0.163299, 0.0909091, 0.157135, 0.142857, 0.353553, 1, 0.166667, 0.333333, 0.333333,
  1e-99, 1e-99, 0.134687, 0.115385, 0.105263, 0.149071, 0.172005, 0.2, 0.25, 0.125, 0.142857, 0.333333,
  1e-99, 1e-99, 0.117121, 0.0901388, 0.0909091, 0.0721312, 0.153093, 0.0744323, 0.125, 0.144338, 0.2, 0.2,
  1e-99, 1e-99, 0.101835, 0.0815892, 0.0681078, 0.0681818, 0.0468122, 0.096225, 0.061859, 0.0541266, 0.0721688, 0.06415,
  1e-99, 1e-99, 0.0989743, 0.0705973, 0.0629941, 0.0744323, 0.0523783, 0.0740741, 0.0744323, 0.0869565, 0.0909091, 0.0434783,
  1e-99, 1e-99, 0.136364, 0.104757, 0.061859, 0.0765466, 0.0824786, 0.0505076, 0.131533, 0.0785674, 0.0714286, 0.0526316,
  1e-99, 1e-99, 0.0645446, 0.0498626, 0.0372402, 0.05, 0.0410959, 0.0408163, 0.0520016, 0.0181818, 0.0479394, 0.0239697,
  1e-99, 1e-99, 0.0418701, 0.030792, 0.0300045, 0.0308288, 0.031427, 0.0341852, 0.0376018, 0.0326374, 0.0352137, 0.0278708,
  1e-99, 1e-99, 0.0647414, 0.0459836, 0.0366606, 0.0444786, 0.0465041, 0.0523783, 0.0652174, 0.0367328, 0.0523783, 0.06,
  1e-99, 1e-99, 0.0460815, 0.029729, 0.0318517, 0.0353764, 0.0322152, 0.043212, 0.040404, 0.0392837, 0.0398049, 0.036348,
  1e-99, 1e-99, 0.0445172, 0.0342483, 0.0288671, 0.0302122, 0.0397566, 0.0346479, 0.0307095, 0.038207, 0.0408228, 0.0339199,
  1e-99, 1e-99, 0.0360844, 0.0257009, 0.0245572, 0.0275391, 0.0263772, 0.0329648, 0.045637, 0.0297275, 0.0322581, 0.0381802,
  1e-99, 1e-99, 0.051379, 0.030261, 0.0310135, 0.0358727, 0.0333986, 0.0401556, 0.0409722, 0.0396508, 0.0384615, 0.0310816,
  1e-99, 1e-99, 0.0698831, 0.0424313, 0.043589, 0.0535468, 0.0534939, 0.0769231, 0.0433013, 0.0363636, 0.0445362, 0.0425532,
  1e-99, 1e-99, 0.0433902, 0.033395, 0.0294579, 0.0332756, 0.0288669, 0.0334077, 0.0355312, 0.0314961, 0.029554, 0.0269975,
  1e-99, 1e-99, 0.0826797, 0.0769231, 0.060813, 0.0771308, 0.0524864, 0.0614875, 0.083189, 0.0753066, 0.0833333, 0.0787296,
  1e-99, 1e-99, 0.0824786, 0.0655738, 0.048766, 0.0661438, 0.0789673, 0.097532, 0.0666173, 0.0745356, 0.0505076, 0.0416667,
  1e-99, 1e-99, 0.08, 0.0903508, 0.0778162, 0.0433013, 0.0844652, 0.0753066, 0.111111, 0.0707107, 0.0721688, 0.0588235,
  1e-99, 1e-99, 0.0851064, 0.057735, 0.0523783, 0.0720438, 0.123718, 0.0952381, 0.0769231, 0.0588235, 0.0769231, 0.0588235,
  1e-99, 1e-99, 0.116642, 0.0555556, 0.09375, 0.0874818, 0.125, 0.142857, 0.142857, 0.111111, 0.125, 0.111111,
  1e-99, 1e-99, 0.333333, 0.174964, 0.11547, 0.333333, 0.25, 0.5, 0.5, 1, 0.2, 0.333333,
  1e-99, 1e-99, 0.131533, 0.166667, 0.172005, 0.0942809, 0.166667, 0.34641, 0.142857, 0.202031, 0.2, 0.0909091,
  1e-99, 1e-99, 0.25, 0.5, 0.25, 1, 1, 0.5, 1, 1e-99, 1, 1,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.411837, 0.312444, 0.269316, 0.26212, 0.23873, 0.233933, 0.202235, 0.2, 0.184055, 0.113127
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0159345, 0.0117841, 0.0109043, 0.0118173, 0.0123116, 0.0141583, 0.015032, 0.01283, 0.0134594, 0.0109879
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  0.285714, 0.2, 0.213333, 0.297521, 0.308108, 0.253369, 0.229773, 0.238938, 0.229901, 0.253464, 0.2607, 0.284307, 0.277027, 0.286723, 0.271093, 0.239875, 0.235172, 0.204545, 0.230198, 0.207031, 0.258621, 0.197279, 0.309091, 0.183673, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.142857, 0.0894427, 0.0533333, 0.0495868, 0.0408099, 0.026133, 0.0272691, 0.0325154, 0.0180073, 0.0120972, 0.0183884, 0.0141448, 0.0136814, 0.0112462, 0.0136922, 0.0193297, 0.0127353, 0.0257703, 0.0238704, 0.0284379, 0.0298629, 0.0366338, 0.0749656, 0.0432923, 0.0588235, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.2, 0.333333, 1e-99, 1e-99, 1e-99, 1, 1, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.142857, 0.2, 0.5, 1e-99, 0.5, 1, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.1875, 0.4, 0.181818, 0.222222, 1e-99, 0.5, 1e-99, 0.166667, 1e-99, 1e-99,
  1e-99, 1e-99, 0.428571, 0.307692, 0.210526, 0.4, 0.461538, 1e-99, 0.25, 0.125, 0.142857, 1e-99,
  1e-99, 1e-99, 0.481481, 0.45, 0.318182, 0.16129, 0.3125, 0.105263, 0.25, 0.333333, 0.2, 1e-99,
  1e-99, 1e-99, 0.488889, 0.365385, 0.25, 0.272727, 0.189189, 0.166667, 0.178571, 0.09375, 0.125, 0.148148,
  1e-99, 1e-99, 0.457143, 0.358491, 0.238095, 0.184211, 0.0740741, 0.185185, 0.157895, 0.217391, 0.181818, 1e-99,
  1e-99, 1e-99, 0.454545, 0.37037, 0.214286, 0.25, 0.238095, 0.178571, 0.294118, 0.166667, 1e-99, 0.105263,
  1e-99, 1e-99, 0.43038, 0.342593, 0.278846, 0.225, 0.164384, 0.142857, 0.186047, 0.0363636, 0.186441, 0.0847458,
  1e-99, 1e-99, 0.378378, 0.306667, 0.292969, 0.281818, 0.227778, 0.236842, 0.126214, 0.212121, 0.181102, 0.12605,
  1e-99, 1e-99, 0.418605, 0.300885, 0.24, 0.285714, 0.176471, 0.277778, 0.23913, 0.220779, 0.185185, 0.16,
  1e-99, 1e-99, 0.454545, 0.291457, 0.274775, 0.298343, 0.270718, 0.305085, 0.252525, 0.259259, 0.202128, 0.16092,
  1e-99, 1e-99, 0.402597, 0.355649, 0.226131, 0.26943, 0.346405, 0.226891, 0.166667, 0.272059, 0.217822, 0.115385,
  1e-99, 1e-99, 0.408333, 0.317808, 0.309859, 0.267797, 0.26087, 0.278075, 0.262712, 0.213483, 0.225806, 0.193878,
  1e-99, 1e-99, 0.506098, 0.267943, 0.298246, 0.294118, 0.261628, 0.188525, 0.181818, 0.196262, 0.217949, 0.0879121,
  1e-99, 1e-99, 0.423729, 0.259259, 0.3, 0.25974, 0.225806, 0.307692, 0.1, 0.127273, 0.181818, 0.0851064,
  1e-99, 1e-99, 0.398734, 0.288991, 0.259091, 0.247059, 0.165563, 0.196429, 0.247191, 0.181102, 0.149533, 0.0816327,
  1e-99, 1e-99, 0.21875, 0.25641, 0.25, 0.325581, 0.0909091, 0.0869565, 0.176471, 0.173913, 0.166667, 0.136364,
  1e-99, 1e-99, 0.333333, 0.393443, 0.172414, 0.1, 0.261905, 0.37931, 0.269231, 0.133333, 0.142857, 0.0833333,
  1e-99, 1e-99, 0.16, 0.285714, 0.323529, 0.225, 0.241379, 0.130435, 0.111111, 0.2, 0.125, 0.0588235,
  1e-99, 1e-99, 0.510638, 0.316667, 0.166667, 0.176471, 0.285714, 0.285714, 1e-99, 0.176471, 0.230769, 0.0588235,
  1e-99, 1e-99, 0.333333, 0.111111, 0.28125, 0.214286, 0.25, 1e-99, 0.285714, 1e-99, 1e-99, 0.111111,
  1e-99, 1e-99, 0.666667, 0.5, 0.2, 1e-99, 0.25, 1e-99, 0.5, 1e-99, 0.2, 1e-99,
  1e-99, 1e-99, 0.235294, 0.25, 0.230769, 0.133333, 1e-99, 0.6, 0.142857, 0.285714, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.2, 0.333333, 1, 1e-99, 1e-99, 1, 1, 1, 0.5,
  1e-99, 1e-99, 0.2, 0.142857, 0.2, 0.5, 1, 0.5, 1, 1e-99, 0.5, 1e-99,
  1e-99, 1e-99, 0.108253, 0.163299, 0.128565, 0.157135, 0.142857, 0.353553, 1, 0.166667, 0.333333, 0.333333,
  1e-99, 1e-99, 0.142857, 0.108786, 0.105263, 0.163299, 0.188422, 0.2, 0.25, 0.125, 0.142857, 0.333333,
  1e-99, 1e-99, 0.133539, 0.106066, 0.120261, 0.0721312, 0.139754, 0.0744323, 0.176777, 0.166667, 0.2, 0.2,
  1e-99, 1e-99, 0.104231, 0.083825, 0.0625, 0.0787296, 0.0715068, 0.096225, 0.0798596, 0.0541266, 0.0721688, 0.0740741,
  1e-99, 1e-99, 0.114286, 0.0822434, 0.0752923, 0.069625, 0.0523783, 0.0828173, 0.0911606, 0.0972203, 0.0909091, 0.0434783,
  1e-99, 1e-99, 0.14374, 0.117121, 0.0874818, 0.0883883, 0.106479, 0.0798596, 0.131533, 0.096225, 0.0714286, 0.0744323,
  1e-99, 1e-99, 0.0738095, 0.0563219, 0.0517804, 0.053033, 0.0474534, 0.0539949, 0.0657774, 0.025713, 0.056214, 0.0378995,
  1e-99, 1e-99, 0.0452249, 0.0369183, 0.0338291, 0.0357909, 0.0355729, 0.0394737, 0.0350054, 0.035855, 0.0377625, 0.0325461,
  1e-99, 1e-99, 0.0697674, 0.0516013, 0.0438178, 0.0539949, 0.0509427, 0.0717219, 0.0721005, 0.0535468, 0.0585607, 0.0565685,
  1e-99, 1e-99, 0.0586816, 0.0382702, 0.0351813, 0.0405993, 0.038674, 0.0508475, 0.0505051, 0.0489954, 0.0463713, 0.0430076,
  1e-99, 1e-99, 0.0511299, 0.0385755, 0.0337096, 0.0373632, 0.0475824, 0.0436651, 0.0392837, 0.0447262, 0.0464398, 0.0384615,
  1e-99, 1e-99, 0.0412479, 0.0295078, 0.0295439, 0.0301295, 0.0307438, 0.038562, 0.0471844, 0.0346315, 0.0381683, 0.0444786,
  1e-99, 1e-99, 0.0555514, 0.0358053, 0.0361676, 0.0396588, 0.0390012, 0.0393101, 0.0454545, 0.0428278, 0.0528603, 0.0310816,
  1e-99, 1e-99, 0.0847458, 0.0489954, 0.0547723, 0.0580797, 0.0603493, 0.0888231, 0.05, 0.0481046, 0.057496, 0.0425532,
  1e-99, 1e-99, 0.0502358, 0.0364094, 0.0343174, 0.038122, 0.0331126, 0.0418787, 0.0527013, 0.0377625, 0.0373832, 0.0288615,
  1e-99, 1e-99, 0.0826797, 0.081084, 0.0693375, 0.0870153, 0.0524864, 0.0614875, 0.101885, 0.0869565, 0.0833333, 0.0787296,
  1e-99, 1e-99, 0.0890871, 0.0803111, 0.054522, 0.05, 0.0789673, 0.114366, 0.10176, 0.0666667, 0.0714286, 0.0416667,
  1e-99, 1e-99, 0.08, 0.0903508, 0.0975478, 0.075, 0.0912328, 0.0753066, 0.111111, 0.1, 0.0721688, 0.0588235,
  1e-99, 1e-99, 0.104234, 0.0726483, 0.0555556, 0.0720438, 0.142857, 0.116642, 0.0769231, 0.101885, 0.133235, 0.0588235,
  1e-99, 1e-99, 0.125988, 0.0785674, 0.09375, 0.0874818, 0.176777, 0.142857, 0.202031, 0.111111, 0.125, 0.111111,
  1e-99, 1e-99, 0.333333, 0.188982, 0.11547, 0.333333, 0.25, 0.5, 0.5, 1, 0.2, 0.333333,
  1e-99, 1e-99, 0.117647, 0.144338, 0.133235, 0.0942809, 0.166667, 0.34641, 0.142857, 0.202031, 0.2, 0.0909091,
  1e-99, 1e-99, 0.25, 0.5, 0.25, 1, 1, 0.5, 1, 1e-99, 1, 1,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.127367, 0.108857, 0.102737, 0.115861, 0.105749, 0.090081, 0.102426, 0.0912951, 0.121076, 0.0699588
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.00740303, 0.00651709, 0.00678903, 0.00807238, 0.00866328, 0.00954856, 0.0117491, 0.00984461, 0.0134529, 0.0119978
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  0.25, 0.125, 0.19403, 0.240741, 0.168605, 0.119658, 0.0866667, 0.0801887, 0.077266, 0.102238, 0.0936639, 0.105691, 0.120338, 0.129136, 0.114754, 0.0860927, 0.08967, 0.0996564, 0.0625, 0.0871369, 0.109489, 0.122137, 0.166667, 0.164835, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.144338, 0.0721688, 0.0538142, 0.0472131, 0.0313091, 0.0184636, 0.0169967, 0.0194486, 0.0107149, 0.00786449, 0.0113584, 0.00883833, 0.00920246, 0.00770361, 0.00904389, 0.0119389, 0.00802033, 0.0185057, 0.0127578, 0.0190148, 0.0199899, 0.0305344, 0.0555556, 0.0425603, 0.0666667, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.5, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.166667, 1, 1e-99, 1e-99, 1e-99, 0.5, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.125, 0.222222, 0.125, 0.111111, 0.2, 1, 0.5, 0.2, 1e-99, 1e-99,
  1e-99, 1e-99, 0.291667, 0.24, 0.285714, 0.333333, 0.1875, 0.166667, 0.5, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0857143, 0.26087, 0.15, 0.125, 0.235294, 0.0833333, 1e-99, 0.25, 0.25, 1e-99,
  1e-99, 1e-99, 0.15873, 0.184615, 0.148148, 0.106383, 0.0294118, 0.0869565, 1e-99, 1e-99, 0.3, 0.0714286,
  1e-99, 1e-99, 0.132075, 0.0517241, 0.0952381, 0.0857143, 0.030303, 1e-99, 0.105263, 0.1875, 0.230769, 1e-99,
  1e-99, 1e-99, 0.0909091, 0.148148, 0.0666667, 0.178571, 0.05, 1e-99, 0.0714286, 1e-99, 1e-99, 0.111111,
  1e-99, 1e-99, 0.127273, 0.116279, 0.0729167, 0.0786517, 0.031746, 0.027027, 0.0277778, 1e-99, 0.0869565, 0.0416667,
  1e-99, 1e-99, 0.116541, 0.0795455, 0.10219, 0.108597, 0.113924, 0.0991736, 0.123457, 0.106557, 0.126582, 0.0298507,
  1e-99, 1e-99, 0.109244, 0.0634921, 0.0859375, 0.0813953, 0.0806452, 0.115385, 0.151515, 0.111111, 0.108108, 0.103448,
  1e-99, 1e-99, 0.0841584, 0.0758929, 0.12037, 0.14094, 0.0987654, 0.102804, 0.136842, 0.0804598, 0.193548, 0.0612245,
  1e-99, 1e-99, 0.141463, 0.159091, 0.0909091, 0.106599, 0.125874, 0.0909091, 0.0769231, 0.11, 0.142857, 0.102041,
  1e-99, 1e-99, 0.135057, 0.11271, 0.115493, 0.145522, 0.137931, 0.128378, 0.185567, 0.118421, 0.148515, 0.0862069,
  1e-99, 1e-99, 0.169421, 0.0991736, 0.1, 0.12069, 0.10274, 0.0982143, 0.0533333, 0.123457, 0.114754, 0.12,
  1e-99, 1e-99, 0.0978261, 0.0894309, 0.0795455, 0.152778, 0.0983607, 0.0238095, 0.0357143, 0.097561, 0.0285714, 0.0454545,
  1e-99, 1e-99, 0.138393, 0.0852713, 0.0986547, 0.0837989, 0.0813008, 0.0609756, 0.0561798, 0.0736842, 0.0675676, 0.0638298,
  1e-99, 1e-99, 0.12766, 0.0769231, 0.0681818, 0.162791, 0.0571429, 0.0454545, 0.1, 0.125, 0.181818, 0.0909091,
  1e-99, 1e-99, 0.078125, 0.0714286, 0.0392157, 0.0638298, 0.147059, 0.0344828, 0.0909091, 1e-99, 0.047619, 1e-99,
  1e-99, 1e-99, 0.027027, 0.205128, 0.097561, 1e-99, 0.172414, 0.0909091, 0.0909091, 1e-99, 0.111111, 1e-99,
  1e-99, 1e-99, 0.16, 0.107143, 0.04, 0.0714286, 0.214286, 0.133333, 1e-99, 0.0769231, 0.0909091, 0.2,
  1e-99, 1e-99, 0.08, 1e-99, 0.1875, 0.227273, 0.1, 1e-99, 0.125, 1e-99, 1e-99, 0.2,
  1e-99, 1e-99, 0.2, 0.1875, 0.166667, 0.5, 1e-99, 1e-99, 0.5, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.210526, 0.142857, 0.25, 0.1, 1e-99, 0.5, 0.142857, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.5, 0.353553, 0.5, 1e-99, 1e-99, 1e-99, 1, 1, 1, 1,
  1e-99, 1e-99, 0.125, 0.166667, 1, 0.5, 0.5, 0.333333, 0.5, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0883883, 0.111111, 0.125, 0.111111, 0.2, 0.707107, 0.5, 0.2, 1e-99, 0.5,
  1e-99, 1e-99, 0.11024, 0.0979796, 0.142857, 0.166667, 0.108253, 0.166667, 0.5, 0.25, 0.25, 1,
  1e-99, 1e-99, 0.0494872, 0.0753066, 0.0866025, 0.0721688, 0.117647, 0.0833333, 0.25, 0.176777, 0.25, 0.5,
  1e-99, 1e-99, 0.0501949, 0.0532939, 0.0523783, 0.0475759, 0.0294118, 0.0614875, 0.0555556, 0.0434783, 0.173205, 0.0714286,
  1e-99, 1e-99, 0.0499198, 0.0298629, 0.047619, 0.0494872, 0.030303, 0.0454545, 0.0744323, 0.108253, 0.133235, 0.111111,
  1e-99, 1e-99, 0.0524864, 0.0740741, 0.0471405, 0.0798596, 0.05, 0.037037, 0.0714286, 0.047619, 0.333333, 0.111111,
  1e-99, 1e-99, 0.0340151, 0.0300231, 0.0275599, 0.0297275, 0.0224478, 0.027027, 0.0277778, 0.0232558, 0.0434783, 0.0416667,
  1e-99, 1e-99, 0.0209314, 0.0173582, 0.0193121, 0.0221673, 0.0268522, 0.0286289, 0.0390405, 0.0295537, 0.0400288, 0.0211077,
  1e-99, 1e-99, 0.0302988, 0.0224478, 0.0259111, 0.0307646, 0.0360656, 0.0471056, 0.0677596, 0.0453609, 0.0540541, 0.0597259,
  1e-99, 1e-99, 0.0204114, 0.0184067, 0.0236066, 0.0307555, 0.0246914, 0.0309965, 0.0379532, 0.0304109, 0.0558726, 0.035348,
  1e-99, 1e-99, 0.0262691, 0.0245483, 0.020856, 0.0232618, 0.0296688, 0.030303, 0.0314037, 0.0331662, 0.043073, 0.045634,
  1e-99, 1e-99, 0.0197002, 0.0164404, 0.018037, 0.0233022, 0.024383, 0.029452, 0.0437386, 0.0279121, 0.0383464, 0.0385529,
  1e-99, 1e-99, 0.0264592, 0.0202437, 0.0213201, 0.0263366, 0.0265273, 0.0296127, 0.0266667, 0.0390405, 0.043373, 0.0489898,
  1e-99, 1e-99, 0.0326087, 0.0269644, 0.0300654, 0.0460642, 0.0401556, 0.0238095, 0.0357143, 0.0487805, 0.0285714, 0.0454545,
  1e-99, 1e-99, 0.0248561, 0.0181799, 0.0210333, 0.0216368, 0.0257096, 0.0272691, 0.0251244, 0.02785, 0.0302171, 0.0368521,
  1e-99, 1e-99, 0.0521168, 0.0384615, 0.0393648, 0.0615291, 0.0404061, 0.0454545, 0.1, 0.0883883, 0.128565, 0.0909091,
  1e-99, 1e-99, 0.0349386, 0.0319438, 0.0277297, 0.0368521, 0.0657667, 0.0344828, 0.0642824, 0.0333333, 0.047619, 0.0625,
  1e-99, 1e-99, 0.027027, 0.0725238, 0.0487805, 0.0294118, 0.0771058, 0.0909091, 0.0909091, 0.0588235, 0.111111, 0.0769231,
  1e-99, 1e-99, 0.046188, 0.0437409, 0.0282843, 0.0505076, 0.123718, 0.0942809, 0.142857, 0.0769231, 0.0909091, 0.2,
  1e-99, 1e-99, 0.0565685, 0.0625, 0.0765466, 0.101639, 0.1, 0.2, 0.125, 0.333333, 0.2, 0.2,
  1e-99, 1e-99, 0.141421, 0.108253, 0.117851, 0.5, 0.25, 0.25, 0.5, 1, 0.5, 1,
  1e-99, 1e-99, 0.105263, 0.101015, 0.125, 0.1, 0.166667, 0.288675, 0.142857, 0.125, 0.333333, 0.5,
  1e-99, 1e-99, 0.2, 0.25, 0.333333, 1e-99, 1e-99, 1, 1, 1, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.463742, 0.333765, 0.289639, 0.261979, 0.264574, 0.237915, 0.190211, 0.186699, 0.177451, 0.104414
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0203757, 0.0146932, 0.012876, 0.0127681, 0.0140619, 0.015017, 0.0145458, 0.012231, 0.0131898, 0.00941472
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  0.272727, 0.25, 0.271429, 0.330579, 0.242424, 0.239544, 0.254032, 0.226244, 0.22449, 0.251326, 0.229907, 0.284306, 0.264798, 0.287396, 0.286856, 0.245315, 0.218123, 0.214521, 0.243323, 0.240664, 0.182353, 0.226277, 0.296296, 0.241758, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.157459, 0.102062, 0.06227, 0.0522691, 0.0383306, 0.0301797, 0.0320051, 0.0319958, 0.0187728, 0.0129098, 0.02073, 0.0146815, 0.0143607, 0.0118722, 0.0159611, 0.0204429, 0.0132257, 0.0266081, 0.0268706, 0.0316007, 0.0327516, 0.0406406, 0.0740741, 0.051543, 0.0588235, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1, 1e-99, 0.666667, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.5, 1e-99, 0.166667, 1e-99, 1e-99, 0.333333, 1e-99, 0.5, 1e-99, 1e-99,
  1e-99, 1e-99, 0.2, 0.526316, 0.375, 0.0909091, 0.333333, 1e-99, 1e-99, 0.2, 1e-99, 1e-99,
  1e-99, 1e-99, 0.578947, 0.375, 0.277778, 0.357143, 0.526316, 0.125, 1e-99, 0.142857, 0.0909091, 1e-99,
  1e-99, 1e-99, 0.428571, 0.421053, 0.296296, 0.0769231, 0.157895, 0.388889, 0.1, 0.176471, 0.142857, 0.125,
  1e-99, 1e-99, 0.521739, 0.357143, 0.305556, 0.297297, 0.214286, 0.142857, 0.238095, 0.103448, 0.1, 0.0740741,
  1e-99, 1e-99, 0.535714, 0.482759, 0.25, 0.307692, 0.0740741, 0.25, 0.133333, 0.26087, 0.136364, 0.0588235,
  1e-99, 1e-99, 0.571429, 0.272727, 0.26087, 0.266667, 0.4, 0.235294, 0.181818, 0.125, 1e-99, 0.0952381,
  1e-99, 1e-99, 0.5, 0.316456, 0.304878, 0.219512, 0.234375, 0.157895, 0.190476, 0.1, 0.114754, 0.0909091,
  1e-99, 1e-99, 0.41129, 0.357143, 0.335329, 0.276042, 0.25, 0.18797, 0.145631, 0.179641, 0.209302, 0.116129,
  1e-99, 1e-99, 0.395833, 0.19697, 0.207317, 0.283582, 0.297297, 0.25, 0.233333, 0.160714, 0.137255, 0.185185,
  1e-99, 1e-99, 0.504673, 0.361111, 0.3125, 0.271676, 0.231707, 0.275229, 0.263636, 0.28, 0.222222, 0.123967,
  1e-99, 1e-99, 0.438095, 0.378049, 0.237179, 0.243902, 0.308219, 0.242991, 0.12766, 0.236111, 0.220183, 0.147368,
  1e-99, 1e-99, 0.47929, 0.302491, 0.314286, 0.312741, 0.281938, 0.280899, 0.235294, 0.198953, 0.219512, 0.168067,
  1e-99, 1e-99, 0.544715, 0.295082, 0.335443, 0.28481, 0.350427, 0.257426, 0.134146, 0.217822, 0.150685, 0.120879,
  1e-99, 1e-99, 0.461538, 0.37037, 0.303371, 0.238095, 0.230769, 0.288462, 0.216216, 0.103448, 0.152542, 0.0701754,
  1e-99, 1e-99, 0.481818, 0.315385, 0.25, 0.166667, 0.2, 0.237113, 0.186275, 0.165354, 0.128713, 0.0666667,
  1e-99, 1e-99, 0.28, 0.272727, 0.243243, 0.418605, 0.172414, 0.0416667, 0.208333, 0.142857, 0.185185, 0.075,
  1e-99, 1e-99, 0.32, 0.487179, 0.22449, 0.285714, 0.291667, 0.25, 0.333333, 0.0645161, 0.21875, 0.0434783,
  1e-99, 1e-99, 0.4, 0.333333, 0.277778, 0.136364, 0.459459, 0.26087, 0.117647, 0.15, 0.157895, 0.0588235,
  1e-99, 1e-99, 0.347826, 0.24, 0.21875, 0.181818, 0.0555556, 0.133333, 1e-99, 0.153846, 0.333333, 1e-99,
  1e-99, 1e-99, 0.4375, 0.15, 0.272727, 0.25, 0.2, 0.285714, 0.285714, 1e-99, 0.166667, 0.0909091,
  1e-99, 1e-99, 0.6, 0.363636, 0.625, 0.142857, 0.2, 1e-99, 0.142857, 1e-99, 1, 1e-99,
  1e-99, 1e-99, 0.428571, 0.428571, 0.142857, 0.222222, 0.181818, 0.25, 0.333333, 0.333333, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1, 0.5, 0.471405, 1, 1, 1e-99, 1, 1e-99, 1, 1,
  1e-99, 1e-99, 0.288675, 1, 0.166667, 1, 0.5, 0.333333, 1, 0.5, 1e-99, 0.5,
  1e-99, 1e-99, 0.141421, 0.166436, 0.216506, 0.0909091, 0.235702, 0.333333, 1, 0.2, 0.25, 0.333333,
  1e-99, 1e-99, 0.174559, 0.153093, 0.124226, 0.159719, 0.166436, 0.125, 0.25, 0.142857, 0.0909091, 0.2,
  1e-99, 1e-99, 0.174964, 0.148865, 0.104757, 0.0543928, 0.0911606, 0.146986, 0.1, 0.101885, 0.142857, 0.125,
  1e-99, 1e-99, 0.150613, 0.112938, 0.0921285, 0.0896385, 0.123718, 0.0714286, 0.106479, 0.0597259, 0.0707107, 0.0523783,
  1e-99, 1e-99, 0.138321, 0.129023, 0.0883883, 0.108786, 0.0523783, 0.144338, 0.0942809, 0.1065, 0.0787296, 0.0415945,
  1e-99, 1e-99, 0.202031, 0.11134, 0.1065, 0.0942809, 0.141421, 0.117647, 0.0909091, 0.0625, 0.05, 0.0673435,
  1e-99, 1e-99, 0.0980581, 0.0632911, 0.0609756, 0.0517395, 0.0605154, 0.0644603, 0.0673435, 0.0408248, 0.043373, 0.0343604,
  1e-99, 1e-99, 0.0575922, 0.0442981, 0.0448103, 0.0379172, 0.040032, 0.037594, 0.0376018, 0.0327978, 0.0402803, 0.0273719,
  1e-99, 1e-99, 0.0908104, 0.0546296, 0.0502818, 0.0650582, 0.0896385, 0.0753778, 0.0881917, 0.0535714, 0.0518775, 0.0585607,
  1e-99, 1e-99, 0.0686773, 0.0500771, 0.0421375, 0.0396281, 0.0375879, 0.0502498, 0.048956, 0.0473286, 0.0496904, 0.0320081,
  1e-99, 1e-99, 0.0645936, 0.0480122, 0.0389921, 0.0385644, 0.0459466, 0.0476544, 0.0368521, 0.0404927, 0.0449448, 0.0393859,
  1e-99, 1e-99, 0.0532544, 0.0328098, 0.0315869, 0.034749, 0.0352423, 0.0397251, 0.0415945, 0.0322744, 0.0365854, 0.037581,
  1e-99, 1e-99, 0.0665476, 0.0491803, 0.0460766, 0.042457, 0.0547276, 0.0504853, 0.0404466, 0.0464398, 0.0454332, 0.0364464,
  1e-99, 1e-99, 0.108786, 0.0676201, 0.0583837, 0.0614759, 0.0666173, 0.0744804, 0.076444, 0.0422326, 0.0508475, 0.0350877,
  1e-99, 1e-99, 0.0661828, 0.0492548, 0.0376889, 0.0340207, 0.04, 0.0494416, 0.0427343, 0.0360833, 0.0356985, 0.0222222,
  1e-99, 1e-99, 0.10583, 0.0909091, 0.0810811, 0.0986661, 0.0771058, 0.0416667, 0.0931695, 0.0824786, 0.0828173, 0.0433013,
  1e-99, 1e-99, 0.113137, 0.111767, 0.0676862, 0.0824786, 0.11024, 0.0944911, 0.125988, 0.0456198, 0.0826797, 0.0307438,
  1e-99, 1e-99, 0.163299, 0.136083, 0.087841, 0.0787296, 0.111435, 0.1065, 0.083189, 0.0866025, 0.0911606, 0.0415945,
  1e-99, 1e-99, 0.122975, 0.0979796, 0.0826797, 0.128565, 0.0555556, 0.0942809, 0.125, 0.108786, 0.19245, 0.0625,
  1e-99, 1e-99, 0.165359, 0.0866025, 0.11134, 0.102062, 0.11547, 0.202031, 0.202031, 0.111111, 0.166667, 0.0909091,
  1e-99, 1e-99, 0.34641, 0.181818, 0.279508, 0.142857, 0.2, 0.333333, 0.142857, 1, 1, 0.166667,
  1e-99, 1e-99, 0.174964, 0.174964, 0.142857, 0.157135, 0.128565, 0.25, 0.333333, 0.19245, 0.5, 0.0555556,
  1e-99, 1e-99, 0.5, 1, 0.5, 1, 0.333333, 0.333333, 1, 1e-99, 0.5, 0.5,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.117908, 0.0952009, 0.0915179, 0.0956772, 0.0860778, 0.0708245, 0.0701262, 0.0692841, 0.0713101, 0.0269608
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.00722939, 0.00609462, 0.00639188, 0.00742599, 0.0080268, 0.00865259, 0.00991735, 0.00894453, 0.0108747, 0.00812898
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  0.142857, 0.142857, 0.136364, 0.207547, 0.119048, 0.0919881, 0.0752688, 0.0738916, 0.0614693, 0.0864899, 0.0773639, 0.0898792, 0.0969131, 0.111374, 0.0981685, 0.070568, 0.0777358, 0.075, 0.0513514, 0.0553191, 0.0730769, 0.103175, 0.176471, 0.126437, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.101015, 0.0824786, 0.0454545, 0.0442492, 0.0266199, 0.0165216, 0.016425, 0.0190787, 0.00959989, 0.00738933, 0.0105279, 0.00823921, 0.00834095, 0.00726527, 0.00848047, 0.0110209, 0.00765954, 0.0163663, 0.0117808, 0.0153428, 0.016765, 0.0286155, 0.0588235, 0.0381221, 0.0625, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.5, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 0.25, 0.5, 1e-99, 1e-99, 1e-99, 1, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.142857, 0.2, 0.142857, 0.111111, 1e-99, 0.5, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.272727, 0.222222, 0.214286, 0.25, 0.153846, 1e-99, 0.333333, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0625, 0.209302, 0.125, 0.0434783, 0.266667, 0.0833333, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.152542, 0.145161, 0.135593, 0.0217391, 0.03125, 0.0909091, 1e-99, 1e-99, 0.111111, 1e-99,
  1e-99, 1e-99, 0.130435, 0.0508475, 0.075, 0.0967742, 0.0285714, 1e-99, 0.111111, 0.125, 0.0909091, 1e-99,
  1e-99, 1e-99, 0.0967742, 0.166667, 0.0645161, 0.133333, 0.05, 1e-99, 1e-99, 1e-99, 1e-99, 0.125,
  1e-99, 1e-99, 0.10084, 0.113821, 0.0480769, 0.0526316, 0.0163934, 0.0285714, 1e-99, 1e-99, 0.075, 1e-99,
  1e-99, 1e-99, 0.118367, 0.0643939, 0.0966543, 0.092233, 0.0745342, 0.0869565, 0.1125, 0.0877193, 0.0547945, 0.0175439,
  1e-99, 1e-99, 0.114035, 0.0454545, 0.0743802, 0.0714286, 0.0892857, 0.0740741, 0.0666667, 0.0980392, 0.0714286, 0.0714286,
  1e-99, 1e-99, 0.0769231, 0.0688073, 0.114679, 0.121019, 0.083871, 0.0792079, 0.09375, 0.075, 0.112903, 0.047619,
  1e-99, 1e-99, 0.12381, 0.134831, 0.0809524, 0.0918919, 0.102941, 0.0652174, 0.0779221, 0.0625, 0.0821918, 0.0212766,
  1e-99, 1e-99, 0.141538, 0.106796, 0.101408, 0.124528, 0.127273, 0.0965517, 0.131313, 0.0833333, 0.0824742, 0.0208333,
  1e-99, 1e-99, 0.144, 0.0962343, 0.0952381, 0.106509, 0.0769231, 0.0825688, 0.0142857, 0.113924, 0.0769231, 0.0681818,
  1e-99, 1e-99, 0.0888889, 0.0655738, 0.0777778, 0.132353, 0.0740741, 0.0232558, 1e-99, 0.0769231, 0.03125, 1e-99,
  1e-99, 1e-99, 0.12, 0.0733591, 0.1, 0.0764706, 0.0720721, 0.04, 0.0365854, 0.060241, 0.0461538, 1e-99,
  1e-99, 1e-99, 0.113636, 0.0535714, 0.06, 0.0769231, 0.0769231, 0.05, 0.0909091, 0.0769231, 0.222222, 1e-99,
  1e-99, 1e-99, 0.0833333, 0.0684932, 0.0181818, 0.0625, 0.0882353, 0.0384615, 0.05, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0333333, 0.136364, 0.0416667, 1e-99, 0.0833333, 0.1, 1e-99, 1e-99, 0.125, 1e-99,
  1e-99, 1e-99, 0.105263, 0.0877193, 0.0212766, 0.0434783, 0.125, 0.0769231, 1e-99, 0.0909091, 1e-99, 1e-99,
  1e-99, 1e-99, 0.08, 1e-99, 0.1875, 0.181818, 0.142857, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.181818, 0.2, 0.166667, 0.5, 1e-99, 1e-99, 1, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.157895, 0.133333, 0.117647, 0.111111, 1e-99, 0.428571, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.5, 0.353553, 0.5, 1e-99, 1e-99, 0.5, 1, 0.5, 1e-99, 1,
  1e-99, 1e-99, 0.142857, 0.25, 0.5, 0.5, 0.5, 0.333333, 1, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.101015, 0.1, 0.142857, 0.111111, 0.2, 0.5, 0.5, 0.2, 1e-99, 0.5,
  1e-99, 1e-99, 0.11134, 0.0907218, 0.123718, 0.125, 0.108786, 0.25, 0.333333, 0.333333, 0.25, 1e-99,
  1e-99, 1e-99, 0.0441942, 0.0697674, 0.0721688, 0.0434783, 0.133333, 0.0833333, 0.2, 0.111111, 0.5, 0.333333,
  1e-99, 1e-99, 0.0508475, 0.0483871, 0.0479394, 0.0217391, 0.03125, 0.0642824, 0.0526316, 0.0555556, 0.111111, 0.0909091,
  1e-99, 1e-99, 0.0532498, 0.0293568, 0.0433013, 0.0558726, 0.0285714, 0.05, 0.0785674, 0.0883883, 0.0909091, 0.333333,
  1e-99, 1e-99, 0.0558726, 0.0833333, 0.0456198, 0.0666667, 0.05, 0.037037, 0.0833333, 0.0666667, 0.2, 0.125,
  1e-99, 1e-99, 0.0291101, 0.03042, 0.0215007, 0.0235376, 0.0163934, 0.0285714, 0.030303, 0.027027, 0.0433013, 0.05,
  1e-99, 1e-99, 0.0219803, 0.0156178, 0.0189555, 0.0211597, 0.0215162, 0.0274981, 0.0375, 0.0277393, 0.0273973, 0.0175439,
  1e-99, 1e-99, 0.0316276, 0.0185567, 0.0247934, 0.0291606, 0.0399298, 0.037037, 0.0471405, 0.0438445, 0.0505076, 0.0505076,
  1e-99, 1e-99, 0.0198615, 0.017766, 0.0229358, 0.0277637, 0.0232616, 0.0280042, 0.03125, 0.0306186, 0.0426734, 0.0336718,
  1e-99, 1e-99, 0.024281, 0.0224719, 0.0196338, 0.0222871, 0.0275122, 0.0266249, 0.0318116, 0.0255155, 0.0335547, 0.0212766,
  1e-99, 1e-99, 0.0208687, 0.0161001, 0.0169014, 0.0216776, 0.0240523, 0.0258045, 0.0364197, 0.0240563, 0.029159, 0.0208333,
  1e-99, 1e-99, 0.024, 0.0200662, 0.0212959, 0.0251044, 0.0231932, 0.0275229, 0.0142857, 0.0379747, 0.0384615, 0.0393648,
  1e-99, 1e-99, 0.031427, 0.0231838, 0.0293972, 0.0441176, 0.037037, 0.0232558, 0.0416667, 0.0444116, 0.03125, 0.0526316,
  1e-99, 1e-99, 0.023094, 0.0168297, 0.0213201, 0.0212091, 0.0254813, 0.023094, 0.0211226, 0.0269406, 0.0266469, 0.0285714,
  1e-99, 1e-99, 0.0508197, 0.0309295, 0.034641, 0.0444116, 0.0543928, 0.05, 0.0909091, 0.0769231, 0.157135, 0.0833333,
  1e-99, 1e-99, 0.0372678, 0.0306311, 0.0181818, 0.0360844, 0.0509427, 0.0384615, 0.05, 0.037037, 0.0666667, 0.0833333,
  1e-99, 1e-99, 0.0333333, 0.0556702, 0.0294628, 0.0277778, 0.0589256, 0.1, 0.0909091, 0.0588235, 0.125, 0.142857,
  1e-99, 1e-99, 0.0372161, 0.0392293, 0.0212766, 0.0434783, 0.0883883, 0.0769231, 0.2, 0.0909091, 0.1, 0.5,
  1e-99, 1e-99, 0.0565685, 0.0526316, 0.0765466, 0.0909091, 0.142857, 0.2, 0.142857, 0.5, 0.25, 0.333333,
  1e-99, 1e-99, 0.128565, 0.11547, 0.117851, 0.5, 0.25, 0.333333, 1, 1e-99, 0.5, 1,
  1e-99, 1e-99, 0.0911606, 0.0942809, 0.083189, 0.111111, 0.2, 0.247436, 0.2, 0.2, 0.5, 0.333333,
  1e-99, 1e-99, 0.2, 0.2, 0.333333, 1e-99, 1, 1, 1, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

)
