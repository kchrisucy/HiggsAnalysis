import FWCore.ParameterSet.Config as cms

tauIDFactorizationCoefficients = cms.untracked.PSet(
factorizationSourceName = cms.untracked.string('histograms_TTToHplusBWB_M90_Winter10'),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.317169, 0.240249, 0.166875, 0.168, 0.119886, 0.109961, 0.101786, 0.0944099, 0.0808709, 0.0440587
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0144168, 0.0116538, 0.0102126, 0.0115931, 0.0106803, 0.011927, 0.0134818, 0.0108296, 0.0112148, 0.00766964
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 0.25, 0.196078, 0.175, 0.178571, 0.131783, 0.134387, 0.132948, 0.154902, 0.168558, 0.171875, 0.161914, 0.20019, 0.183134, 0.179796, 0.178926, 0.163769, 0.114286, 0.158076, 0.123656, 0.144144, 0.1875, 0.146341, 0.191176, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.2, 0.176777, 0.0620054, 0.0467707, 0.0399298, 0.0226006, 0.0230472, 0.0277216, 0.0174278, 0.0116874, 0.0183219, 0.0122047, 0.0137816, 0.0107355, 0.0129086, 0.0188605, 0.0122407, 0.0233285, 0.023307, 0.025784, 0.0254813, 0.0441942, 0.0597437, 0.0530228, 0.0909091, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.333333, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.333333, 0.0909091, 0.166667, 0.4, 0.25, 1e-99, 1e-99, 1e-99, 0.333333, 1e-99,
  1e-99, 1e-99, 0.235294, 0.352941, 0.0769231, 0.285714, 0.125, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.25, 0.363636, 0.111111, 0.333333, 0.166667, 1e-99, 1e-99, 1e-99, 0.166667, 1e-99,
  1e-99, 1e-99, 0.314286, 0.325581, 0.0638298, 0.0952381, 0.0555556, 1e-99, 0.0769231, 0.0588235, 0.0526316, 1e-99,
  1e-99, 1e-99, 0.307692, 0.166667, 0.0909091, 0.208333, 0.115385, 0.111111, 0.0588235, 0.0588235, 1e-99, 0.04,
  1e-99, 1e-99, 0.363636, 0.266667, 0.0714286, 0.1, 1e-99, 0.125, 0.0909091, 0.0625, 1e-99, 0.047619,
  1e-99, 1e-99, 0.289855, 0.265957, 0.166667, 0.146341, 0.111111, 0.0857143, 0.0869565, 0.0465116, 0.027027, 0.0444444,
  1e-99, 1e-99, 0.322981, 0.253968, 0.169312, 0.165414, 0.130769, 0.1, 0.125, 0.0784314, 0.0843373, 0.0537634,
  1e-99, 1e-99, 0.296296, 0.23913, 0.177215, 0.183099, 0.133333, 0.030303, 0.0740741, 0.0857143, 0.16, 0.025641,
  1e-99, 1e-99, 0.292517, 0.196629, 0.115152, 0.157895, 0.169492, 0.131868, 0.0571429, 0.119403, 0.129032, 0.107143,
  1e-99, 1e-99, 0.325843, 0.264706, 0.258741, 0.155556, 0.157407, 0.118421, 0.0980392, 0.0689655, 0.134615, 0.111111,
  1e-99, 1e-99, 0.312217, 0.236641, 0.192, 0.188119, 0.0965909, 0.14876, 0.127907, 0.135593, 0.0789474, 0.0779221,
  1e-99, 1e-99, 0.310127, 0.259067, 0.132075, 0.157534, 0.125, 0.136986, 0.111111, 0.136986, 0.163934, 0.0454545,
  1e-99, 1e-99, 0.381579, 0.267442, 0.163934, 0.164384, 0.113636, 0.103448, 0.0833333, 0.113636, 0.0333333, 1e-99,
  1e-99, 1e-99, 0.340741, 0.210843, 0.184049, 0.188034, 0.118182, 0.102273, 0.114754, 0.108696, 0.0595238, 0.025974,
  1e-99, 1e-99, 0.269231, 0.151515, 0.166667, 0.0769231, 0.0909091, 0.0666667, 0.0625, 0.117647, 1e-99, 1e-99,
  1e-99, 1e-99, 0.482759, 0.25, 0.2, 0.129032, 0.107143, 0.173913, 0.105263, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.291667, 0.15625, 0.172414, 0.105263, 0.111111, 0.125, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.285714, 0.162162, 0.135135, 0.16, 1e-99, 1e-99, 0.5, 0.111111, 1e-99, 1e-99,
  1e-99, 1e-99, 0.285714, 0.4, 0.227273, 0.25, 1e-99, 0.142857, 1e-99, 0.111111, 1e-99, 1e-99,
  1e-99, 1e-99, 0.4, 0.1, 0.285714, 1e-99, 0.2, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.444444, 0.357143, 1e-99, 0.25, 0.0833333, 1e-99, 1e-99, 0.2, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.5, 1, 1, 1e-99, 1e-99, 1e-99, 1,
  1e-99, 1e-99, 0.333333, 0.333333, 1e-99, 1e-99, 1, 1e-99, 1e-99, 1e-99, 1e-99, 1,
  1e-99, 1e-99, 0.19245, 0.0909091, 0.117851, 0.282843, 0.25, 1, 1, 0.5, 0.333333, 0.333333,
  1e-99, 1e-99, 0.117647, 0.144088, 0.0769231, 0.202031, 0.125, 0.333333, 0.5, 0.25, 0.166667, 0.333333,
  1e-99, 1e-99, 0.111803, 0.128565, 0.0785674, 0.19245, 0.166667, 0.1, 0.2, 0.166667, 0.166667, 0.1,
  1e-99, 1e-99, 0.0947607, 0.0870153, 0.0368521, 0.0673435, 0.0555556, 0.0666667, 0.0769231, 0.0588235, 0.0526316, 0.0333333,
  1e-99, 1e-99, 0.0888231, 0.0680414, 0.0524864, 0.0931695, 0.0666173, 0.0785674, 0.0588235, 0.0588235, 0.0555556, 0.04,
  1e-99, 1e-99, 0.128565, 0.0942809, 0.0714286, 0.0707107, 0.0526316, 0.125, 0.0909091, 0.0625, 0.0833333, 0.047619,
  1e-99, 1e-99, 0.0648136, 0.0531915, 0.046225, 0.0597437, 0.0496904, 0.0494872, 0.0614875, 0.0328887, 0.027027, 0.031427,
  1e-99, 1e-99, 0.0447895, 0.0366572, 0.0299304, 0.0352663, 0.0317162, 0.0333333, 0.0441942, 0.0277297, 0.0318765, 0.0240437,
  1e-99, 1e-99, 0.0604812, 0.0509828, 0.0473628, 0.0507824, 0.0666667, 0.030303, 0.0523783, 0.0494872, 0.08, 0.025641,
  1e-99, 1e-99, 0.0446084, 0.0332364, 0.0264176, 0.0344555, 0.0378995, 0.0380671, 0.0285714, 0.0422153, 0.0456198, 0.0437409,
  1e-99, 1e-99, 0.0427852, 0.03946, 0.0425368, 0.033945, 0.0381769, 0.0394737, 0.0438445, 0.0281551, 0.0508798, 0.0453609,
  1e-99, 1e-99, 0.0375865, 0.0300535, 0.0277128, 0.0305169, 0.0234267, 0.0350631, 0.0385654, 0.0338983, 0.0322301, 0.0318116,
  1e-99, 1e-99, 0.0443038, 0.0366377, 0.0288212, 0.0328482, 0.0360844, 0.0433189, 0.0453609, 0.0433189, 0.0518406, 0.0262432,
  1e-99, 1e-99, 0.0708574, 0.0557655, 0.0518406, 0.0474534, 0.0508197, 0.0597259, 0.0589256, 0.0508197, 0.0333333, 0.0277778,
  1e-99, 1e-99, 0.0502395, 0.035639, 0.0336026, 0.040089, 0.0327777, 0.0340909, 0.043373, 0.0343726, 0.0266199, 0.0183664,
  1e-99, 1e-99, 0.10176, 0.0677596, 0.0745356, 0.0769231, 0.0642824, 0.0666667, 0.0625, 0.083189, 0.0526316, 0.0526316,
  1e-99, 1e-99, 0.129023, 0.0790569, 0.0666667, 0.0645161, 0.061859, 0.0869565, 0.0744323, 0.0434783, 0.0555556, 0.0285714,
  1e-99, 1e-99, 0.11024, 0.0698771, 0.0771058, 0.0744323, 0.0785674, 0.0883883, 0.333333, 0.0714286, 0.1, 0.047619,
  1e-99, 1e-99, 0.0824786, 0.0662024, 0.0604343, 0.08, 0.047619, 0.0909091, 0.25, 0.111111, 0.1, 0.0454545,
  1e-99, 1e-99, 0.142857, 0.2, 0.101639, 0.144338, 0.2, 0.142857, 0.333333, 0.111111, 0.125, 0.166667,
  1e-99, 1e-99, 0.282843, 0.1, 0.202031, 0.333333, 0.2, 0.25, 1e-99, 0.2, 1e-99, 0.5,
  1e-99, 1e-99, 0.222222, 0.159719, 0.25, 0.176777, 0.0833333, 0.25, 0.5, 0.2, 0.333333, 0.142857,
  1e-99, 1e-99, 0.166667, 1, 0.5, 1e-99, 1e-99, 1, 1e-99, 1e-99, 1, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.390564, 0.291125, 0.225, 0.2144, 0.170314, 0.147477, 0.148214, 0.131677, 0.108865, 0.0547397
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0159981, 0.0128285, 0.0118585, 0.0130966, 0.0127299, 0.0138125, 0.0162686, 0.0127896, 0.0130118, 0.0085489
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 0.25, 0.176471, 0.2125, 0.205357, 0.139535, 0.162055, 0.196532, 0.215686, 0.209887, 0.216797, 0.218951, 0.248577, 0.250472, 0.247451, 0.228628, 0.202196, 0.142857, 0.19244, 0.129032, 0.175676, 0.21875, 0.146341, 0.191176, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.2, 0.176777, 0.0588235, 0.0515388, 0.0428199, 0.0232558, 0.0253088, 0.0337049, 0.0205649, 0.0130417, 0.0205774, 0.0141925, 0.0153571, 0.012555, 0.0151438, 0.0213197, 0.0136012, 0.026082, 0.0257159, 0.0263386, 0.0281306, 0.0477352, 0.0597437, 0.0530228, 0.0909091, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.333333, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.333333, 0.0909091, 0.166667, 0.2, 0.25, 1e-99, 1e-99, 1e-99, 0.333333, 1e-99,
  1e-99, 1e-99, 0.352941, 0.352941, 0.0769231, 0.285714, 0.25, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.35, 0.409091, 0.0555556, 0.333333, 0.166667, 0.1, 1e-99, 1e-99, 0.166667, 1e-99,
  1e-99, 1e-99, 0.228571, 0.348837, 0.0638298, 0.0952381, 0.0555556, 0.0666667, 0.153846, 0.0588235, 0.105263, 0.0333333,
  1e-99, 1e-99, 0.384615, 0.222222, 0.0909091, 0.25, 0.115385, 0.0555556, 0.117647, 0.0588235, 1e-99, 0.08,
  1e-99, 1e-99, 0.5, 0.266667, 0.285714, 0.1, 0.105263, 0.125, 0.181818, 0.1875, 1e-99, 0.047619,
  1e-99, 1e-99, 0.42029, 0.340426, 0.25641, 0.243902, 0.155556, 0.0857143, 1e-99, 0.0697674, 0.108108, 0.0444444,
  1e-99, 1e-99, 0.347826, 0.328042, 0.232804, 0.218045, 0.176923, 0.133333, 0.140625, 0.0784314, 0.120482, 0.0645161,
  1e-99, 1e-99, 0.308642, 0.282609, 0.240506, 0.225352, 0.166667, 0.121212, 0.148148, 0.2, 0.16, 0.025641,
  1e-99, 1e-99, 0.360544, 0.252809, 0.230303, 0.233083, 0.186441, 0.153846, 0.1, 0.19403, 0.145161, 0.107143,
  1e-99, 1e-99, 0.410112, 0.329412, 0.27972, 0.17037, 0.240741, 0.157895, 0.156863, 0.091954, 0.192308, 0.111111,
  1e-99, 1e-99, 0.420814, 0.305344, 0.268, 0.252475, 0.164773, 0.22314, 0.174419, 0.194915, 0.0921053, 0.0779221,
  1e-99, 1e-99, 0.417722, 0.321244, 0.194969, 0.232877, 0.166667, 0.232877, 0.222222, 0.164384, 0.196721, 0.0757576,
  1e-99, 1e-99, 0.421053, 0.302326, 0.229508, 0.191781, 0.159091, 0.103448, 0.166667, 0.204545, 0.0666667, 0.111111,
  1e-99, 1e-99, 0.407407, 0.228916, 0.245399, 0.239316, 0.190909, 0.113636, 0.147541, 0.130435, 0.0833333, 0.012987,
  1e-99, 1e-99, 0.346154, 0.181818, 0.166667, 1e-99, 0.181818, 0.0666667, 0.125, 0.176471, 1e-99, 1e-99,
  1e-99, 1e-99, 0.62069, 0.275, 0.222222, 0.193548, 0.142857, 0.173913, 0.105263, 0.0434783, 1e-99, 1e-99,
  1e-99, 1e-99, 0.333333, 0.1875, 0.172414, 0.105263, 0.0555556, 0.125, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.380952, 0.27027, 0.108108, 0.12, 0.047619, 1e-99, 0.625, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.357143, 0.3, 0.272727, 0.25, 0.2, 0.142857, 1e-99, 0.111111, 0.125, 1e-99,
  1e-99, 1e-99, 0.4, 1e-99, 0.428571, 1e-99, 0.2, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.555556, 0.285714, 1e-99, 0.25, 0.0833333, 1e-99, 1e-99, 0.2, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionShrinkingConeTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 0.5, 1, 1, 1e-99, 1e-99, 1e-99, 1,
  1e-99, 1e-99, 0.333333, 0.333333, 1e-99, 1e-99, 1, 1e-99, 1e-99, 1e-99, 1e-99, 1,
  1e-99, 1e-99, 0.19245, 0.0909091, 0.117851, 0.2, 0.25, 1, 1, 0.5, 0.333333, 0.333333,
  1e-99, 1e-99, 0.144088, 0.144088, 0.0769231, 0.202031, 0.176777, 0.333333, 0.5, 0.25, 0.166667, 0.333333,
  1e-99, 1e-99, 0.132288, 0.136364, 0.0555556, 0.19245, 0.166667, 0.1, 0.2, 0.166667, 0.166667, 0.1,
  1e-99, 1e-99, 0.0808122, 0.0900694, 0.0368521, 0.0673435, 0.0555556, 0.0666667, 0.108786, 0.0588235, 0.0744323, 0.0333333,
  1e-99, 1e-99, 0.0993073, 0.0785674, 0.0524864, 0.102062, 0.0666173, 0.0555556, 0.083189, 0.0588235, 0.0555556, 0.0565685,
  1e-99, 1e-99, 0.150756, 0.0942809, 0.142857, 0.0707107, 0.0744323, 0.125, 0.128565, 0.108253, 0.0833333, 0.047619,
  1e-99, 1e-99, 0.0780459, 0.0601793, 0.0573351, 0.0771287, 0.0587945, 0.0494872, 0.0434783, 0.0402803, 0.0540541, 0.031427,
  1e-99, 1e-99, 0.0464802, 0.0416614, 0.0350966, 0.04049, 0.036891, 0.03849, 0.046875, 0.0277297, 0.0380997, 0.0263386,
  1e-99, 1e-99, 0.0617284, 0.0554241, 0.0551759, 0.056338, 0.0745356, 0.0606061, 0.0740741, 0.0755929, 0.08, 0.025641,
  1e-99, 1e-99, 0.0495246, 0.0376865, 0.0373601, 0.0418629, 0.0397493, 0.0411171, 0.0377964, 0.0538142, 0.0483871, 0.0437409,
  1e-99, 1e-99, 0.048, 0.0440195, 0.0442277, 0.0355247, 0.0472131, 0.0455803, 0.0554594, 0.0325107, 0.060813, 0.0453609,
  1e-99, 1e-99, 0.0436364, 0.0341384, 0.0327414, 0.0353536, 0.0305975, 0.0429434, 0.0450347, 0.0406426, 0.0348125, 0.0318116,
  1e-99, 1e-99, 0.051418, 0.040798, 0.0350174, 0.039938, 0.0416667, 0.0564809, 0.06415, 0.0474534, 0.0567886, 0.0338798,
  1e-99, 1e-99, 0.0744323, 0.0592909, 0.0613386, 0.0512556, 0.0601307, 0.0597259, 0.0833333, 0.0681818, 0.0471405, 0.0555556,
  1e-99, 1e-99, 0.0549348, 0.037135, 0.038801, 0.0452265, 0.0416598, 0.035935, 0.0491803, 0.0376533, 0.031497, 0.012987,
  1e-99, 1e-99, 0.115385, 0.074227, 0.0745356, 0.0769231, 0.0909091, 0.0666667, 0.0883883, 0.101885, 0.0526316, 0.0526316,
  1e-99, 1e-99, 0.146298, 0.0829156, 0.0702728, 0.0790158, 0.0714286, 0.0869565, 0.0744323, 0.0434783, 0.0555556, 0.0285714,
  1e-99, 1e-99, 0.117851, 0.0765466, 0.0771058, 0.0744323, 0.0555556, 0.0883883, 0.333333, 0.0714286, 0.1, 0.047619,
  1e-99, 1e-99, 0.0952381, 0.085467, 0.0540541, 0.069282, 0.047619, 0.0909091, 0.279508, 0.111111, 0.1, 0.0454545,
  1e-99, 1e-99, 0.159719, 0.173205, 0.11134, 0.144338, 0.2, 0.142857, 0.333333, 0.111111, 0.125, 0.166667,
  1e-99, 1e-99, 0.282843, 0.1, 0.247436, 0.333333, 0.2, 0.25, 1e-99, 0.2, 1e-99, 0.5,
  1e-99, 1e-99, 0.248452, 0.142857, 0.25, 0.176777, 0.0833333, 0.25, 0.5, 0.2, 0.333333, 0.142857,
  1e-99, 1e-99, 0.166667, 1, 0.5, 1e-99, 1e-99, 1, 1e-99, 1e-99, 1, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.130685, 0.114809, 0.0848187, 0.104515, 0.0785463, 0.0762565, 0.0707071, 0.0821918, 0.0694087, 0.0670927
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.00790939, 0.00744199, 0.00722025, 0.00934811, 0.00959596, 0.0114961, 0.0133624, 0.0126825, 0.0133577, 0.0146408
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 0.285714, 0.155556, 0.132353, 0.128713, 0.0720339, 0.0675105, 0.0363636, 0.0851064, 0.102165, 0.106472, 0.0918972, 0.126506, 0.108652, 0.100691, 0.104701, 0.0917874, 0.0777202, 0.102564, 0.0639535, 0.0891089, 0.126437, 0.0769231, 0.169231, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.2, 0.202031, 0.0587945, 0.0441176, 0.0356985, 0.0174708, 0.0168776, 0.0148454, 0.0134565, 0.009405, 0.014909, 0.0095293, 0.0112701, 0.0085365, 0.0099699, 0.0149573, 0.00941719, 0.0200673, 0.0193828, 0.0192827, 0.0210032, 0.0381221, 0.0444116, 0.051025, 0.1, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionHPSTauBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1, 0.25, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.181818, 0.0909091, 0.125, 0.4, 1e-99, 1e-99, 1e-99, 1e-99, 0.333333, 1e-99,
  1e-99, 1e-99, 0.136364, 0.285714, 0.111111, 0.125, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.08, 0.238095, 0.117647, 0.285714, 0.0909091, 1e-99, 1e-99, 1e-99, 0.25, 1e-99,
  1e-99, 1e-99, 0.156863, 0.122449, 0.025641, 1e-99, 1e-99, 1e-99, 1e-99, 0.0769231, 0.1, 1e-99,
  1e-99, 1e-99, 0.108696, 0.1, 1e-99, 0.107143, 0.047619, 1e-99, 1e-99, 0.1, 1e-99, 0.111111,
  1e-99, 1e-99, 0.0333333, 0.0322581, 1e-99, 0.0357143, 1e-99, 0.25, 0.125, 0.1, 1e-99, 1e-99,
  1e-99, 1e-99, 0.132653, 0.115789, 0.0487805, 0.0588235, 0.0952381, 0.0357143, 0.0555556, 0.0434783, 1e-99, 0.142857,
  1e-99, 1e-99, 0.126697, 0.119342, 0.102564, 0.0967742, 0.11828, 0.0759494, 0.06, 0.0540541, 0.0540541, 0.0769231,
  1e-99, 1e-99, 0.182796, 0.101852, 0.0506329, 0.134328, 0.0740741, 0.0322581, 0.0588235, 0.115385, 0.176471, 1e-99,
  1e-99, 1e-99, 0.125628, 0.0938967, 0.0545455, 0.0900901, 0.102041, 0.0757576, 0.0217391, 0.130435, 0.0714286, 0.153846,
  1e-99, 1e-99, 0.178723, 0.116162, 0.0993377, 0.116667, 0.123711, 0.115385, 0.0909091, 0.0909091, 0.135135, 0.0689655,
  1e-99, 1e-99, 0.122977, 0.123779, 0.104, 0.101523, 0.037594, 0.117647, 0.107143, 0.0945946, 0.162162, 0.139535,
  1e-99, 1e-99, 0.104072, 0.107143, 0.0745342, 0.137931, 0.106667, 0.122449, 0.0789474, 0.0784314, 0.06, 0.107143,
  1e-99, 1e-99, 0.129032, 0.171717, 0.046875, 0.109589, 0.097561, 0.105263, 0.0526316, 0.0769231, 1e-99, 1e-99,
  1e-99, 1e-99, 0.107843, 0.115578, 0.113636, 0.0992366, 0.0434783, 0.047619, 0.0697674, 0.106383, 0.0377358, 1e-99,
  1e-99, 1e-99, 0.0888889, 0.129032, 0.12, 1e-99, 0.0909091, 1e-99, 0.0769231, 0.111111, 1e-99, 1e-99,
  1e-99, 1e-99, 0.148936, 0.0892857, 0.183673, 0.142857, 0.05, 0.105263, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.147059, 0.0606061, 0.107143, 1e-99, 1e-99, 0.1, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.122807, 0.0714286, 0.0526316, 0.133333, 1e-99, 1e-99, 0.363636, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.2, 1e-99, 0.2, 0.375, 1e-99, 1e-99, 1e-99, 0.142857, 1e-99, 1e-99,
  1e-99, 1e-99, 0.2, 1e-99, 1e-99, 1e-99, 0.2, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.166667, 0.263158, 1e-99, 0.181818, 0.333333, 1e-99, 1e-99, 0.25, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionHPSTauBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 0.5, 1, 1, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1,
  1e-99, 1e-99, 1, 0.25, 1e-99, 1e-99, 1e-99, 1, 1e-99, 1, 1e-99, 1e-99,
  1e-99, 1e-99, 0.128565, 0.0909091, 0.125, 0.282843, 0.5, 0.333333, 1e-99, 0.5, 0.333333, 1e-99,
  1e-99, 1e-99, 0.0787296, 0.142857, 0.111111, 0.125, 0.2, 1e-99, 0.25, 0.333333, 0.5, 1,
  1e-99, 1e-99, 0.0565685, 0.106479, 0.083189, 0.202031, 0.0909091, 0.166667, 0.5, 0.2, 0.25, 0.333333,
  1e-99, 1e-99, 0.0554594, 0.0499896, 0.025641, 0.0526316, 0.0454545, 0.0769231, 0.111111, 0.0769231, 0.1, 0.0909091,
  1e-99, 1e-99, 0.0486102, 0.0447214, 0.0277778, 0.061859, 0.047619, 0.0666667, 0.125, 0.1, 0.0714286, 0.111111,
  1e-99, 1e-99, 0.0333333, 0.0322581, 0.0454545, 0.0357143, 0.0769231, 0.25, 0.125, 0.1, 0.142857, 0.0833333,
  1e-99, 1e-99, 0.0367913, 0.0349118, 0.0243902, 0.0339618, 0.047619, 0.0357143, 0.0555556, 0.0434783, 0.0526316, 0.101015,
  1e-99, 1e-99, 0.0239435, 0.0221612, 0.022934, 0.0279363, 0.0356626, 0.0310062, 0.034641, 0.027027, 0.038222, 0.0444116,
  1e-99, 1e-99, 0.0443345, 0.0307095, 0.0253165, 0.0447761, 0.0523783, 0.0322581, 0.0588235, 0.0666173, 0.101885, 0.0714286,
  1e-99, 1e-99, 0.0251256, 0.0209959, 0.0181818, 0.028489, 0.0322681, 0.0338798, 0.0217391, 0.0532498, 0.0412393, 0.0769231,
  1e-99, 1e-99, 0.0275776, 0.0242214, 0.0256489, 0.0311805, 0.0357124, 0.0471056, 0.0524864, 0.0454545, 0.0604343, 0.048766,
  1e-99, 1e-99, 0.0199496, 0.0200795, 0.0203961, 0.0227012, 0.0168125, 0.0372033, 0.0437409, 0.0357534, 0.0662024, 0.0569649,
  1e-99, 1e-99, 0.0217006, 0.0218704, 0.0215162, 0.0344828, 0.0377124, 0.0499896, 0.0455803, 0.0392157, 0.034641, 0.061859,
  1e-99, 1e-99, 0.0372484, 0.0416475, 0.0270633, 0.0387456, 0.0487805, 0.0744323, 0.0526316, 0.0543928, 0.0555556, 0.0625,
  1e-99, 1e-99, 0.0229922, 0.0240997, 0.0254099, 0.0275233, 0.0217391, 0.0274929, 0.0402803, 0.0475759, 0.0266833, 0.037037,
  1e-99, 1e-99, 0.0444444, 0.0645161, 0.069282, 0.0588235, 0.0642824, 0.0833333, 0.0769231, 0.111111, 0.1, 0.111111,
  1e-99, 1e-99, 0.0562926, 0.0399298, 0.0612245, 0.0714286, 0.05, 0.0744323, 0.0769231, 0.0588235, 0.0833333, 0.0833333,
  1e-99, 1e-99, 0.0657667, 0.042855, 0.061859, 0.0357143, 0.0833333, 0.1, 0.25, 0.0769231, 0.25, 0.166667,
  1e-99, 1e-99, 0.0464167, 0.0412393, 0.0372161, 0.0942809, 0.1, 0.0833333, 0.181818, 0.25, 0.166667, 0.142857,
  1e-99, 1e-99, 0.1, 0.0625, 0.11547, 0.216506, 0.125, 0.2, 0.5, 0.142857, 0.5, 0.25,
  1e-99, 1e-99, 0.141421, 0.142857, 0.111111, 0.333333, 0.2, 0.5, 1e-99, 0.5, 1, 1e-99,
  1e-99, 1e-99, 0.117851, 0.117688, 0.166667, 0.128565, 0.333333, 0.5, 0.5, 0.25, 0.25, 0.5,
  1e-99, 1e-99, 0.2, 1, 0.5, 1e-99, 1, 1, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.463636, 0.347478, 0.268147, 0.204225, 0.179864, 0.180473, 0.145078, 0.130491, 0.0894526, 0.0407041
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.0205302, 0.0166795, 0.0151324, 0.0143338, 0.0142642, 0.0163393, 0.0158293, 0.0129843, 0.0109284, 0.00669171
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 0.25, 0.204082, 0.184211, 0.247619, 0.16875, 0.197044, 0.22561, 0.232662, 0.205128, 0.241477, 0.224615, 0.23991, 0.252323, 0.246722, 0.242762, 0.206533, 0.140845, 0.158996, 0.179348, 0.241667, 0.162791, 0.181818, 0.322581, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.166667, 0.176777, 0.0645363, 0.0492323, 0.0485621, 0.032476, 0.0311554, 0.03709, 0.0228144, 0.0139572, 0.0261919, 0.0151781, 0.0163999, 0.0134298, 0.0171484, 0.0232524, 0.0147524, 0.0257147, 0.0257925, 0.0312204, 0.0448764, 0.0435076, 0.0642824, 0.0721312, 0.0909091, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCaloTauCutBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.333333, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.428571, 0.111111, 0.222222, 0.181818, 0.333333, 1e-99, 1e-99, 1e-99, 0.5, 1e-99,
  1e-99, 1e-99, 0.333333, 0.25, 0.444444, 0.0833333, 0.142857, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.5, 0.588235, 0.0714286, 0.25, 0.307692, 1e-99, 1e-99, 1e-99, 0.166667, 1e-99,
  1e-99, 1e-99, 0.421053, 0.458333, 1e-99, 0.0555556, 0.375, 0.0769231, 0.1, 0.0769231, 0.0769231, 1e-99,
  1e-99, 1e-99, 0.454545, 0.357143, 0.181818, 0.1875, 0.2, 0.117647, 0.181818, 0.1, 1e-99, 0.0909091,
  1e-99, 1e-99, 0.777778, 0.4, 0.1875, 0.263158, 0.125, 0.0714286, 0.285714, 0.111111, 0.0714286, 0.037037,
  1e-99, 1e-99, 0.509804, 0.412698, 0.357143, 0.1875, 0.166667, 0.230769, 1e-99, 0.129032, 0.0816327, 0.0454545,
  1e-99, 1e-99, 0.438596, 0.348837, 0.246479, 0.178947, 0.165049, 0.183099, 0.164706, 0.122222, 0.0888889, 0.0447761,
  1e-99, 1e-99, 0.459459, 0.422222, 0.277778, 0.282609, 0.137931, 0.259259, 0.0714286, 0.107143, 0.16, 0.030303,
  1e-99, 1e-99, 0.435644, 0.282895, 0.257812, 0.190476, 0.211538, 0.240506, 0.155844, 0.112676, 0.127907, 0.0972222,
  1e-99, 1e-99, 0.480315, 0.323077, 0.298969, 0.177083, 0.219048, 0.179104, 0.125, 0.121951, 0.109589, 0.0847458,
  1e-99, 1e-99, 0.464088, 0.352941, 0.292079, 0.224138, 0.20438, 0.196429, 0.178571, 0.188976, 0.0568182, 0.0555556,
  1e-99, 1e-99, 0.485981, 0.378151, 0.268519, 0.231707, 0.104651, 0.216216, 0.177778, 0.217391, 0.127907, 0.047619,
  1e-99, 1e-99, 0.461538, 0.358209, 0.315789, 0.305085, 0.139535, 0.0714286, 0.190476, 0.2, 0.0263158, 0.05,
  1e-99, 1e-99, 0.456311, 0.33913, 0.251748, 0.217822, 0.162791, 0.191781, 0.121212, 0.117647, 0.0879121, 1e-99,
  1e-99, 1e-99, 0.222222, 0.333333, 0.25, 1e-99, 0.266667, 0.153846, 0.190476, 0.04, 0.0588235, 1e-99,
  1e-99, 1e-99, 0.434783, 0.263158, 0.307692, 0.166667, 0.173913, 0.0769231, 0.2, 0.04, 1e-99, 0.025,
  1e-99, 1e-99, 0.5, 0.333333, 0.208333, 0.266667, 0.0588235, 0.166667, 0.0909091, 0.0666667, 1e-99, 1e-99,
  1e-99, 1e-99, 0.588235, 0.384615, 0.315789, 0.166667, 0.2, 0.111111, 0.2, 1e-99, 0.142857, 1e-99,
  1e-99, 1e-99, 0.555556, 1e-99, 0.3, 0.111111, 0.181818, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.5, 0.0833333, 0.5, 0.25, 0.142857, 0.2, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.666667, 0.642857, 1e-99, 0.125, 0.2, 1e-99, 1e-99, 1e-99, 0.2, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCaloTauCutBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1, 1e-99, 0.333333, 1e-99, 1, 1e-99, 1,
  1e-99, 1e-99, 0.333333, 0.333333, 1e-99, 1, 1e-99, 1e-99, 1e-99, 1e-99, 1, 1e-99,
  1e-99, 1e-99, 0.247436, 0.111111, 0.157135, 0.128565, 0.333333, 0.333333, 1e-99, 0.333333, 0.5, 0.5,
  1e-99, 1e-99, 0.149071, 0.144338, 0.222222, 0.0833333, 0.142857, 0.333333, 0.2, 0.2, 0.333333, 0.2,
  1e-99, 1e-99, 0.176777, 0.186016, 0.0714286, 0.176777, 0.153846, 0.166667, 0.25, 0.1, 0.166667, 0.0909091,
  1e-99, 1e-99, 0.148865, 0.138193, 0.0588235, 0.0555556, 0.216506, 0.0769231, 0.1, 0.0769231, 0.0769231, 0.04,
  1e-99, 1e-99, 0.14374, 0.112938, 0.0909091, 0.108253, 0.1, 0.083189, 0.128565, 0.0707107, 0.0714286, 0.0524864,
  1e-99, 1e-99, 0.20787, 0.163299, 0.108253, 0.117688, 0.0883883, 0.0714286, 0.202031, 0.0785674, 0.0714286, 0.037037,
  1e-99, 1e-99, 0.0999808, 0.0809368, 0.0798596, 0.0625, 0.0680414, 0.0942111, 0.047619, 0.0645161, 0.0408163, 0.0262432,
  1e-99, 1e-99, 0.0620269, 0.0520016, 0.0416625, 0.0434011, 0.0400302, 0.0507824, 0.0440195, 0.0368514, 0.031427, 0.0182798,
  1e-99, 1e-99, 0.111435, 0.0968644, 0.0717219, 0.0783815, 0.0689655, 0.0979908, 0.0505076, 0.061859, 0.08, 0.030303,
  1e-99, 1e-99, 0.0656757, 0.043141, 0.0448794, 0.0425918, 0.0451002, 0.0551759, 0.0449883, 0.039837, 0.0385654, 0.0367465,
  1e-99, 1e-99, 0.061498, 0.0498519, 0.0555172, 0.042949, 0.0456746, 0.051703, 0.0472456, 0.0385644, 0.0387456, 0.0378995,
  1e-99, 1e-99, 0.0506362, 0.0415945, 0.0380255, 0.0358908, 0.0386241, 0.0418787, 0.0461069, 0.0385746, 0.0254099, 0.0248452,
  1e-99, 1e-99, 0.0673935, 0.0563715, 0.0498626, 0.0531573, 0.0348837, 0.0540541, 0.0628539, 0.0561302, 0.0385654, 0.0274929,
  1e-99, 1e-99, 0.084265, 0.0731191, 0.0911606, 0.0719092, 0.0569649, 0.0505076, 0.0952381, 0.0632456, 0.0263158, 0.0353553,
  1e-99, 1e-99, 0.0665598, 0.0543043, 0.041958, 0.0464398, 0.0435076, 0.0512556, 0.042855, 0.0415945, 0.0310816, 0.00970874,
  1e-99, 1e-99, 0.0907218, 0.125988, 0.111803, 0.0555556, 0.133333, 0.108786, 0.0952381, 0.04, 0.0588235, 0.0277778,
  1e-99, 1e-99, 0.13749, 0.117688, 0.108786, 0.0745356, 0.0869565, 0.0769231, 0.11547, 0.04, 0.04, 0.025,
  1e-99, 1e-99, 0.166667, 0.105409, 0.0931695, 0.133333, 0.0588235, 0.117851, 0.0909091, 0.0666667, 0.111111, 0.030303,
  1e-99, 1e-99, 0.186016, 0.172005, 0.128921, 0.096225, 0.141421, 0.111111, 0.2, 0.166667, 0.142857, 0.0625,
  1e-99, 1e-99, 0.248452, 0.142857, 0.122474, 0.111111, 0.128565, 0.2, 0.2, 0.125, 0.2, 0.142857,
  1e-99, 1e-99, 0.353553, 0.0833333, 0.353553, 0.25, 0.142857, 0.2, 1e-99, 0.333333, 1, 0.25,
  1e-99, 1e-99, 0.235702, 0.214286, 1e-99, 0.125, 0.2, 0.333333, 0.5, 0.2, 0.2, 0.125,
  1e-99, 1e-99, 0.25, 1, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99, 1, 1, 1,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.119199, 0.0991211, 0.0732018, 0.0893322, 0.0675845, 0.0733083, 0.0534759, 0.0594714, 0.0512821, 0.0161943
) ),

tauIDFactorizationByPt_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 0.00763092, 0.00695694, 0.0068261, 0.00880216, 0.00919708, 0.0117387, 0.0119576, 0.0114453, 0.0128205, 0.00809717
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 0.285714, 0.130435, 0.123077, 0.102041, 0.0684932, 0.0511628, 0.0394737, 0.0776053, 0.0840951, 0.0971302, 0.0743381, 0.106472, 0.0876217, 0.0845511, 0.0993228, 0.0882957, 0.0734463, 0.0912548, 0.0552147, 0.0837696, 0.111111, 0.0882353, 0.163934, 1e-99, 1e-99
) ),

tauIDFactorizationByEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  0.25, 0.202031, 0.0532498, 0.0435143, 0.0322681, 0.0176849, 0.0154262, 0.0161151, 0.0131177, 0.00876752, 0.0146429, 0.00870061, 0.0105423, 0.00780596, 0.00939457, 0.0149735, 0.00952117, 0.0203703, 0.0186273, 0.0184049, 0.0209424, 0.037037, 0.0509427, 0.0518406, 0.1, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_Coefficients = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.5, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.166667, 0.1, 0.1, 0.5, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.15, 0.235294, 1e-99, 0.125, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0714286, 0.222222, 0.111111, 0.25, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.142857, 0.1, 0.0277778, 1e-99, 1e-99, 1e-99, 1e-99, 0.0909091, 0.142857, 1e-99,
  1e-99, 1e-99, 0.125, 0.0652174, 1e-99, 0.0769231, 0.0454545, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.0344828, 0.0322581, 1e-99, 0.037037, 1e-99, 0.166667, 0.166667, 0.125, 1e-99, 1e-99,
  1e-99, 1e-99, 0.111111, 0.11828, 0.0519481, 0.0588235, 0.1, 0.04, 1e-99, 0.0454545, 1e-99, 1e-99,
  1e-99, 1e-99, 0.112069, 0.0952381, 0.0806452, 0.0940171, 0.0909091, 0.0724638, 0.0454545, 0.03125, 1e-99, 0.0333333,
  1e-99, 1e-99, 0.173913, 0.0900901, 0.0394737, 0.0983607, 0.0909091, 0.0344828, 0.0588235, 0.157895, 0.142857, 1e-99,
  1e-99, 1e-99, 0.113861, 0.0730594, 0.0454545, 0.0754717, 0.104167, 0.0491803, 1e-99, 0.075, 0.0512821, 0.0434783,
  1e-99, 1e-99, 0.153509, 0.097561, 0.0827586, 0.0884956, 0.116279, 0.122449, 0.0555556, 0.0731707, 0.129032, 1e-99,
  1e-99, 1e-99, 0.116279, 0.101974, 0.0978723, 0.0816327, 0.0162602, 0.0952381, 0.0754717, 0.0416667, 0.0967742, 0.025641,
  1e-99, 1e-99, 0.084507, 0.0963303, 0.060241, 0.110169, 0.0857143, 0.146341, 0.0882353, 0.0243902, 0.0512821, 0.0555556,
  1e-99, 1e-99, 0.111111, 0.161616, 0.0461538, 0.106061, 0.0833333, 0.105263, 0.0625, 0.0714286, 1e-99, 1e-99,
  1e-99, 1e-99, 0.103627, 0.101523, 0.115854, 0.092437, 0.0416667, 0.0535714, 0.0697674, 0.0869565, 0.0512821, 1e-99,
  1e-99, 1e-99, 0.0909091, 0.12, 0.0714286, 1e-99, 0.105263, 1e-99, 0.1, 0.125, 1e-99, 1e-99,
  1e-99, 1e-99, 0.108696, 0.0909091, 0.145833, 0.16, 0.047619, 0.105263, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.178571, 0.0263158, 0.0714286, 1e-99, 1e-99, 0.111111, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.127273, 0.0769231, 0.0526316, 0.125, 1e-99, 1e-99, 0.222222, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.222222, 1e-99, 0.166667, 0.181818, 1e-99, 1e-99, 1e-99, 0.142857, 1e-99, 1e-99,
  1e-99, 1e-99, 0.25, 1e-99, 1e-99, 1e-99, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.153846, 0.277778, 1e-99, 0.181818, 1e-99, 1e-99, 1e-99, 0.2, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

tauIDFactorizationByPtVsEta_signalAnalysisTauSelectionCombinedHPSTaNCBased_CoefficientUncertainty = cms.untracked.vdouble( *(
  1e-99, 1e-99, 1e-99, 1e-99, 1, 1, 1e-99, 0.5, 1e-99, 1e-99, 1e-99, 1e-99,
  1e-99, 1e-99, 0.5, 0.333333, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1, 1e-99, 1,
  1e-99, 1e-99, 0.117851, 0.1, 0.1, 0.353553, 0.5, 0.5, 1e-99, 0.333333, 0.333333, 1e-99,
  1e-99, 1e-99, 0.0866025, 0.117647, 0.125, 0.125, 0.166667, 1e-99, 1e-99, 0.25, 1, 1,
  1e-99, 1e-99, 0.0505076, 0.111111, 0.0785674, 0.176777, 0.1, 0.142857, 0.5, 0.333333, 0.333333, 1,
  1e-99, 1e-99, 0.0539949, 0.0447214, 0.0277778, 0.05, 0.0555556, 0.0909091, 0.1, 0.0909091, 0.142857, 0.142857,
  1e-99, 1e-99, 0.0559017, 0.0376533, 0.0285714, 0.0543928, 0.0454545, 0.0833333, 0.1, 0.2, 0.0833333, 0.142857,
  1e-99, 1e-99, 0.0344828, 0.0322581, 0.0416667, 0.037037, 0.111111, 0.166667, 0.166667, 0.125, 0.2, 0.142857,
  1e-99, 1e-99, 0.0335013, 0.0356626, 0.025974, 0.0339618, 0.05, 0.04, 0.0625, 0.0454545, 0.0714286, 0.0714286,
  1e-99, 1e-99, 0.0219785, 0.0203048, 0.0208225, 0.0283472, 0.0321412, 0.0324068, 0.0321412, 0.0220971, 0.030303, 0.0333333,
  1e-99, 1e-99, 0.0434783, 0.028489, 0.0227901, 0.0401556, 0.0642824, 0.0344828, 0.0588235, 0.0911606, 0.101015, 0.0833333,
  1e-99, 1e-99, 0.0237417, 0.0182648, 0.0171802, 0.0266833, 0.0329404, 0.0283943, 0.0238095, 0.0433013, 0.0362619, 0.0434783,
  1e-99, 1e-99, 0.0259477, 0.0218153, 0.0238904, 0.0279848, 0.0367707, 0.0499896, 0.0392837, 0.0422451, 0.0645161, 0.0416667,
  1e-99, 1e-99, 0.0196548, 0.018315, 0.0204078, 0.0204082, 0.0114977, 0.0336718, 0.0377358, 0.0240563, 0.0558726, 0.025641,
  1e-99, 1e-99, 0.0199185, 0.021021, 0.0190499, 0.0305555, 0.0349927, 0.0597437, 0.0509427, 0.0243902, 0.0362619, 0.0555556,
  1e-99, 1e-99, 0.0351364, 0.040404, 0.0266469, 0.0400871, 0.0481125, 0.0744323, 0.0625, 0.0505076, 0.0666667, 0.111111,
  1e-99, 1e-99, 0.0231717, 0.0227012, 0.0265787, 0.0278708, 0.0208333, 0.0309295, 0.0402803, 0.0434783, 0.0362619, 0.047619,
  1e-99, 1e-99, 0.0454545, 0.069282, 0.0505076, 0.0555556, 0.0744323, 0.0769231, 0.1, 0.125, 0.166667, 0.166667,
  1e-99, 1e-99, 0.0486102, 0.0406558, 0.0551198, 0.08, 0.047619, 0.0744323, 0.0666667, 0.0714286, 0.0909091, 0.111111,
  1e-99, 1e-99, 0.0798596, 0.0263158, 0.0505076, 0.0357143, 0.0909091, 0.111111, 0.142857, 0.125, 1e-99, 0.166667,
  1e-99, 1e-99, 0.0481046, 0.0444116, 0.0372161, 0.0883883, 0.1, 0.0909091, 0.157135, 0.5, 0.25, 0.142857,
  1e-99, 1e-99, 0.111111, 0.0714286, 0.117851, 0.128565, 0.142857, 0.5, 0.333333, 0.142857, 0.333333, 0.25,
  1e-99, 1e-99, 0.176777, 0.166667, 0.111111, 0.333333, 0.333333, 0.5, 1e-99, 0.5, 1, 1e-99,
  1e-99, 1e-99, 0.108786, 0.124226, 0.166667, 0.128565, 0.25, 0.5, 1, 0.2, 1e-99, 1,
  1e-99, 1e-99, 0.2, 1, 0.5, 1e-99, 1e-99, 1, 1e-99, 1e-99, 1, 1e-99,
  1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99, 1e-99
) ),

)
