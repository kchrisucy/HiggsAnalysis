#!/usr/bin/env python
'''
Description:
Generate all TH1 generated from the Kinematics analyzer (GEN-level info).

Usage:
./plot_1d.py -m <pseudo_mcrab> [opts]

Examples:
./plot_1d.py -m <peudo_mcrab> -o "" --url --normaliseToOne
./plot_1d.py -m Kinematics_170828_082301/ -i 'TT' --url --normaliseToOne
./plot_1d.py -m Kinematics_170828_082301/ -e "QCD_b|M_180|M_200|M_220|M_250|M_300|M_350|M_400|M_500|M_1000|M_2000" --url --mergeEWK
./plot_1d.py -m Kinematics_170828_082301/ -e "QCD_b|M_180|M_200|M_220|M_250|M_350|M_400|M_1000|M_2000|M_3000" --mergeEWK --url
./plot_1d.py -m Kinematics_170828_082301/ -e "QCD_b|M_180|M_200|M_220|M_250|M_350|M_400|M_1000|M_2000|M_3000" --url --mergeEWK --normaliseToOne

Last Used:
./plot_1d.py -m Kinematics_FullStats_170831_085353 -e "QCD_b|M_180|M_200|M_220|M_250|M_350|M_400|M_1000|M_2000|M_3000" --url --mergeEWK --normaliseToOne  
./plot_1d.py -m Kinematics_170830_060219 -e "QCD_b|M_180|M_200|M_220|M_250|M_350|M_400|M_1000|M_2000|M_3000" --url --mergeEWK --normaliseToOne

//kchristo////////////////////////////////////////////////////////////////
Triplet for bjet inside barycenter
"Inclusive"
"bJetInside"
"NObJetInside"
example : --folder "Inclusive"
//////////////////////////////////////////////////////////////////////////

'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
from optparse import OptionParser

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

import HiggsAnalysis.NtupleAnalysis.tools.dataset as dataset
import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms
import HiggsAnalysis.NtupleAnalysis.tools.counter as counter
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.styles as styles
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HiggsAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck

#================================================================================================ 
# Function Definition
#================================================================================================ 
def Print(msg, printHeader=False):
    fName = __file__.split("/")[-1]
    if printHeader==True:
        print "=== ", fName
        print "\t", msg
    else:
        print "\t", msg
    return

def Verbose(msg, printHeader=True, verbose=False):
    if not opts.verbose:
        return
    aux.Print(msg, printHeader)
    return

def GetLumi(datasetsMgr):
    Verbose("Determininig Integrated Luminosity")
    
    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            lumi += d.getLuminosity()
    Verbose("Luminosity = %s (pb)" % (lumi), True)
    return lumi


def GetListOfEwkDatasets():
    Verbose("Getting list of EWK datasets")
    return ["TT", "WJetsToQQ_HT_600ToInf", "DYJetsToQQHT", "SingleTop", "TTWJetsToQQ", "TTZToQQ", "Diboson", "TTTT"]


def GetDatasetsFromDir(opts):
    Verbose("Getting datasets")
    
    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode, 
                                                        analysisName=opts.analysisName,
                                                        optimizationMode=opts.optMode)
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        includeOnlyTasks=opts.includeOnlyTasks,
                                                        optimizationMode=opts.optMode)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        excludeTasks=opts.excludeTasks,
                                                        optimizationMode=opts.optMode)
    else:
        raise Exception("This should never be reached")
    return datasets
    

def main(opts, signalMass):

    optModes = [""]

    if opts.optMode != None:
        optModes = [opts.optMode]

    # For-loop: All optimisation modes
    for opt in optModes:
        opts.optMode = opt

        # Setup & configure the dataset manager 
        datasetsMgr = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        datasetsMgr.loadLuminosities() # from lumi.json
        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()

        # Set/Overwrite cross-sections
        for d in datasetsMgr.getAllDatasets():
            if "ChargedHiggs" in d.getName():
                datasetsMgr.getDataset(d.getName()).setCrossSection(1.0)
               
        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 
        
        # Determine integrated Lumi before removing data
        #intLumi = datasetsMgr.getDataset("Data").getLuminosity()
        intLumi = 5747.588 + 2573.399 + 4248.384 + 4008.663 + 2704.118 + 405.222 + 7539.457 + 8390.5 + 215.149

        # Remove datasets
        if 1:   #kchris, it was 0
            #datasetsMgr.remove(filter(lambda name: "Data" in name, datasetsMgr.getAllDatasetNames()))                            //it was uncomment 
            #datasetsMgr.remove(filter(lambda name: "QCD-b" in name, datasetsMgr.getAllDatasetNames()))                           //it was uncomment 
            #datasetsMgr.remove(filter(lambda name: "QCD" in name, datasetsMgr.getAllDatasetNames()))                             //it was uncomment
            #datasetsMgr.remove(filter(lambda name: "SingleTop" in name, datasetsMgr.getAllDatasetNames()))
            #datasetsMgr.remove(filter(lambda name: "DYJetsToQQHT" in name, datasetsMgr.getAllDatasetNames()))
            datasetsMgr.remove(filter(lambda name: "TTZToQQ" in name, datasetsMgr.getAllDatasetNames()))                           #it was comment
            datasetsMgr.remove(filter(lambda name: "TTWJetsToQQ" in name, datasetsMgr.getAllDatasetNames()))                       #it was comment
            #datasetsMgr.remove(filter(lambda name: "WJetsToQQ" in name, datasetsMgr.getAllDatasetNames()))
            #datasetsMgr.remove(filter(lambda name: "Diboson" in name, datasetsMgr.getAllDatasetNames()))
            datasetsMgr.remove(filter(lambda name: "TTTT" in name, datasetsMgr.getAllDatasetNames()))                              #it was comment
            #datasetsMgr.remove(filter(lambda name: "FakeBMeasurementTrijetMass" in name, datasetsMgr.getAllDatasetNames()))      //it was uncomment 
            #datasetsMgr.remove(filter(lambda name: "M_" in name and "M_" + str(opts.signalMass) not in name, datasetsMgr.getAllDatasetNames()))

        # Merge EWK samples
        if opts.mergeEWK:
            datasetsMgr.merge("EWK", GetListOfEwkDatasets())
            plots._plotStyles["EWK"] = styles.getAltEWKStyle()

        # Print dataset information
        datasetsMgr.PrintInfo()


        # Re-order datasets (different for inverted than default=baseline)
        newOrder = []
        if opts.mergeEWK:
            newOrder = ["EWK", "QCD"]
        else:
            newOrder.extend(GetListOfEwkDatasets())
            newOrder.append("QCD")
        for d in datasetsMgr.getAllDatasetNames():
            if "ChargedHiggs" in d:
                newOrder.extend([d])        
        datasetsMgr.selectAndReorder( reversed(newOrder) )


        # Apply TDR style
        style = tdrstyle.TDRStyle()
        style.setOptStat(True)
        style.setGridX(True)
        style.setGridY(False)

        # Do the topSelection histos
        folder     = opts.folder
        histoPaths = []
        histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(folder)
        hList      = [x for x in histoList if "_Vs_" not in x] # remove TH2
        histoPaths = [os.path.join(folder, h) for h in hList]
        skipMe     = ["counters", "Weighting", "config", "configInfo"]
        for h in histoPaths:
            if h in skipMe:
                continue
            PlotMC(datasetsMgr, h, intLumi)
    return


def GetHistoKwargs(histo, opts):
    '''
    Dictionary with 
    key   = histogramName
    value = kwargs
    '''
    
    if opts.normaliseToOne:
        yLabel = "Arbitrary Units"
    else:
        yLabel = "Events"
    logY       = False
    yMaxFactor = 1.2

    # Create with default values
    kwargs = {
        "xlabel"           : "x-label",
        "ylabel"           : yLabel,
        "rebinX"           : 1,
        "rebinY"           : 1,
        "ratioYlabel"      : "Data/MC",
        "ratio"            : True, 
        "stackMCHistograms": True,
        "ratioInvert"      : False, 
        "addMCUncertainty" : False, 
        "addLuminosityText": True,
        "addCmsText"       : True,
        "cmsExtraText"     : "Preliminary",
        "opts"             : {"ymin": 0.0, "ymaxfactor": yMaxFactor},
        "opts2"            : {"ymin": 0.0, "ymax": 2.0},
        "log"              : logY,
        "moveLegend"       : {"dx": -0.1, "dy": 0.0, "dh": -0.1},
        "cutBox"           : {"cutValue": 0.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True},
#        "cutLine"          : {"cutValue": 0.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True} #kchristo
        }
    
    if "eta" in histo.lower(): 
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": +5.0}#, "ymin": 1e+0, "ymaxfactor": 10}

    if "deta" in histo.lower(): 
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#eta"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +5.0}

    if "met_et" in histo.lower():
        units            = "GeV/c"
        format           = "%0.0f " + units
        kwargs["xlabel"] = "E_{T}^{miss} (%s)" % units
        kwargs["ylabel"] = yLabel + "/ %s " % format
        print kwargs["ylabel"]
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 300.0}
        kwargs["log"]    = True
        
    if "ht" in histo.lower():
        units            = "GeV"
        format           = "%0.0f " + units
        kwargs["ylabel"] = yLabel + " / %.0f " + units
        kwargs["xlabel"] = "H_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 500.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 1
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2000, "ymin": 1e+0, "ymaxfactor": 10}
        ROOT.gStyle.SetNdivisions(5, "X")

    if "MHT" in histo:
        units            = "GeV"
        format           = "%0.0f " + units
        kwargs["ylabel"] = yLabel + " / %.0f " + units
        kwargs["xlabel"] = "MHT (%s)"  % units
        kwargs["rebinX"] = 1
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 300, "ymin": 1e+0, "ymaxfactor": 10}
        ROOT.gStyle.SetNdivisions(5, "X")

    if "jt" in histo.lower():
        units            = "GeV"
        kwargs["ylabel"] = yLabel + " / %.0f " + units
        kwargs["xlabel"] = "J_{T} (%s)"  % units
        kwargs["rebinX"] = 1
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1500, "ymin": 1e+0, "ymaxfactor": 10}
        ROOT.gStyle.SetNdivisions(5, "X")

    if "alphat" in histo.lower():
        units            = ""
        kwargs["xlabel"] = "#alpha_{T}"
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0, "ymin": 1e-5, "ymaxfactor": 10}
        kwargs["log"]    = True

    if "y23" in histo.lower():
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["xlabel"] = "y_{23}"

    if "sphericity" in histo.lower():
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["xlabel"] = "Sphericity"

    if "planarity" in histo.lower():
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["xlabel"] = "Planarity"

    if "aplanarity" in histo.lower():
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["xlabel"] = "Aplanarity"

    if "centrality" in histo.lower():
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["xlabel"] = "Centrality"
        kwargs["moveLegend"] = {"dx": -0.53, "dy": 0.0, "dh": -0.1}

    if "circularity" in histo.lower():
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["xlabel"] = "Circularity"

    if "h2" in histo.lower():
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["xlabel"] = "H_{2}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.01}
        kwargs["cutBox"] = {"cutValue": 0.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "cparameter" in histo.lower():
        kwargs["ylabel"]      = yLabel + " / %.2f"
        kwargs["xlabel"]     = "C"
        kwargs["opts"]       = {"xmin": 0.0, "xmax": 1.01}
        kwargs["moveLegend"] = {"dx": -0.53, "dy": 0.0, "dh": -0.1}

    if "dparameter" in histo.lower():
        kwargs["ylabel"] = yLabel + " / %.2f"
        kwargs["xlabel"] = "D"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.01}
        ROOT.gStyle.SetNdivisions(10, "X")

    if  histo == "Y":
        units            = ""
        kwargs["xlabel"] = "Y"
        kwargs["ylabel"] = yLabel + " / %.1f"

        # kchristo +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if "DiJet_BJet_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 8.0}

    if "GenTopPt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
        # ----------Study Boosted Topologies----------------------------
    if "bfromH_Higgs_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "objectsfromHiggstop_maxdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{max}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                   
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "objectsfromHiggstop_mindR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "objectsfromHiggstop_Prob_mindR_lt_p8" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}
        
    if "bfromHiggstop_Higgstop_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                               
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "b_top_fromhiggs_underdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

        # try with fat jets
    #if "Hs_topQuark_fatjet_mindR" in histo:
     #   kwargs["ylabel"] = yLabel + " / %.1f"
      #  kwargs["xlabel"] = "#DeltaR_{min}"
       # kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        #kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
    
        
        
        #

    if "b_higgs_underdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

        
    if "WfromHiggstop_closestb_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "WfromHiggstop_samedecayb_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "for_WbfromH_intofatjet_W_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "for_WbfromH_intofatjet_Top_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Probdecayproductsintofatjet_Hs" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}
        
    if "Hs_QuarksintofatjetMultiplicity" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 6.0}
        
    if "Hs_isbQuarkintofatjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Pairsinbarycenter_enoughdeltaR_Hs" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}
        
    if "WNOTfromHiggstop_closestb_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "WNOTfromHiggstop_samedecayb_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}


    if "objectsNOTfromHiggstop_maxdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{max}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                               
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "objectsNOTfromHiggstop_mindR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "objectsNOTfromHiggstop_Prob_mindR_lt_p8" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}
        
        #for both
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}
        
        # 12.01
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7_tf" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "  
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
        #

    if "bNOTfromHiggstop_top_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                               
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "b_topNOTfromhiggs_underdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "for_WbNOTfromH_intofatjet_W_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "for_WbNOTfromH_intofatjet_Top_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Probdecayproductsintofatjet_NOTHs" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}

    if "NotHs_QuarksintofatjetMultiplicity" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 6.0}
        
    if "NotHs_isbQuarkintofatjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "Pairsinbarycenter_enoughdeltaR_NOTHs" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}

        # try with gen-jets++++++++++++++++++++++++++++++++
    if "Hs_deltaR_bQuark_closestbJet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}

    if "Hs_deltaR_bQuark_fromH_closestbJet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}
        
    if "Hs_deltaR_QfromW_closestJet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}
        
    if "Hs_numof_obj_fromtop_matchedwith_genjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 3.5}
        
    if "Hs_ProbdecayJetsintofatjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}
 
    ##############################################try for both Hs and NotHs######################################
    if "JetsintofatjetMultiplicity" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 6.0}
        
    if "isbJetintofatjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}        
    ##############################################################################################################

    if "Hs_BjetInsideFatJet_Top_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 2
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "Hs_BjetInsideFatJet_W_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 2
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_BjetInsideFatJet_H_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 2
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "NotHs_deltaR_bQuark_closestbJet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}

    if "NotHs_deltaR_QfromW_closestJet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}

    if "NotHs_numof_obj_fromtop_matchedwith_genjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 3.5}

    if "NotHs_ProbdecayJetsintofatjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}
        
    if "NotHs_BjetInsideFatJet_Top_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 2
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "NotHs_BjetInsideFatJet_W_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 2
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        #++++++++++++++++++++++++++++++++++++++++++++++++++

        # Update on boosted Topologies
    if "Hs_objectsfromtop_top_maxdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{max}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "Hs_objectsfromtop_top_mindR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "Hs_which_objectfromtop_maxdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "NotHs_objectsfromtop_top_maxdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{max}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "NotHs_objectsfromtop_top_mindR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "NotHs_which_objectfromtop_maxdR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
       

    if "Iffatjet_Hs_Top_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Iffatjet_Hs_W_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Iffatjet_Hs_bfromH_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Iffatjet_Hs_EventsWithHighTop_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Iffatjet_Hs_EventsWithHighW_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Iffatjet_Hs_HighWandTop_pT_dRmin" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Iffatjet_Hs_bfromT_Top_dEta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#eta"
        #kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}
        
    if "Iffatjet_Hs_bfromT_Top_dPhi" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#phi"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}
        
    if "Iffatjet_Hs_bfromT_Top_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}
        
    if "Iffatjet_Hs_bfromT_W_dEta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#eta"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_Hs_bfromT_W_dPhi" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#phi"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}
        
    if "Iffatjet_Hs_bfromT_W_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}
        
    if "Iffatjet_Hs_bfromT_bfromH_dEta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#eta"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_Hs_bfromT_bfromH_dPhi" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#phi"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}
        
    if "Iffatjet_Hs_bfromT_bfromH_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}


    if "Iffatjet_NotHs_Top_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "Iffatjet_NotHs_W_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "Iffatjet_NotHs_bfromH_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "Iffatjet_NotHs_EventsWithHighTop_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Iffatjet_NotHs_EventsWithHighW_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "Iffatjet_NotHs_HighWandTop_pT_dRmin" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Iffatjet_NotHs_bfromT_Top_dEta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#eta"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_NotHs_bfromT_Top_dPhi" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#phi"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_NotHs_bfromT_Top_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_NotHs_bfromT_W_dEta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#eta"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_NotHs_bfromT_W_dPhi" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#phi"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}
        
    if "Iffatjet_NotHs_bfromT_W_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_NotHs_bfromT_bfromH_dEta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#eta"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_NotHs_bfromT_bfromH_dPhi" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#Delta#phi"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

    if "Iffatjet_NotHs_bfromT_bfromH_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0}

        # Update 4/12/2017 on Boosted Topologies
    if "mostdistantfromtop_isb__top_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "mostdistantfromtop_isb__W_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "mostdistantfromtop_isb__b_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 200.0}
        
    if "mostdistantfromtop_isb__dRmax_b_top" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "mostdistantfromtop_isb_dRmin_b_objfromW" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
        # ----Study Boosted Topologies With fat Jets -------------------------
        # try for both Hs and not
    if "Hs_QuarksFromW_deltaR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "Hs_QuarksFromW_Prob_deltaR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_QuarksintoBaryCenterMultiplicity" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "Hs_isbQuarkintoBaryCenter" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}
        
    if "Hs_OnlyQQ_dR_less_p7" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_Prob_Diquark_match_with_fj" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_MatchedWithDiquark_fj_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_MatchedWithDiquark_Prob_fj_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_MatchedWithDiquark_Prob_fj_pT_pTcuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "Hs_objectsfromtop_dRmax_pTcuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{max}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "Hs_objectsfromtop_dRmin_pTcuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "Hs_objectsfromtop_Prob_dRmax_pTcuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{max}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_objectsfromtop_Prob_dRmin_pTcuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_QuarksFromTop_Passed_pTcuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
        
    if "Hs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "Hs_objectsfromtop_mindR_ltp8_matchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 500.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 500.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 200.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 300.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
        # -------------Update 24.01.2018 after pt cuts---------------------------------------
    if "Hs_Bc_MassOf_bq" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "m_{bq} (%s)" % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 350.0}
        
    if "Hs_Bc_bqq_Top_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_Bc_qq_Top_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_Bc_qq_bq_Top_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_Bc_bq_Top_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_Bc_bqq_W_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 300.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "Hs_Bc_qq_W_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 300.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "Hs_Bc_qq_bq_W_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 300.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_Bc_bq_W_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 300.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "Hs_bqcase_deltaR_bq" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{bq}"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}
        
    if "Hs_bqcase_Prob_deltaR_lt" in histo:
         kwargs["ylabel"] = yLabel + " / %.1f"
         kwargs["xlabel"] = "#DeltaR_{bq}"
         kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
         
    if "Hs_BoostedW_deltaR_qq" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{qq}"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}

    if "Hs_BoostedW_Prob_deltaR_lt" in histo:
         kwargs["ylabel"] = yLabel + " / %.1f"
         kwargs["xlabel"] = "#DeltaR_{qq}"
         kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "Hs_BoostedWcase_matchWithAk4" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}

    if "Hs_otherBcloseToTopProd" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "Other b close"
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 3.5}
        
    if "Hs_otherBclose_BoostedTop" in histo:
         kwargs["ylabel"] = yLabel + " / %.1f"
         kwargs["xlabel"] = " "
         kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
         
    if "Hs_otherBclose_BoostedW" in histo:
         kwargs["ylabel"] = yLabel + " / %.1f"
         kwargs["xlabel"] = " "
         kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Hs_MatchedWithDiquark_fj_NumOf_Subjets" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "Subjets" 
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 3.5}
        
    if "Hs_MatchedWithDiquark_fj_hasBsubjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "Has b-subjet"
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 1.5}
        
    if "Hs_MatchedWithDiquark_fj_CSV" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "CSV"
        kwargs["cutBox"] = {"cutValue": 0.8484, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}
        
    if "Hs_QuarksintoFatJetMultiplicity" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}
        
    if "Hs_BoostedWinFatJet_dR_qq" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{qq}"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "Hs_BoostedWinFatJet_Prob_deltaR_lt" in histo:
         kwargs["ylabel"] = yLabel + " / %.1f"
         kwargs["xlabel"] = "#DeltaR_{qq}"
         kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
             
        # n jettiness --------------------------------------------------------
    if "InFatJet_fatjet_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "InFatJet_Higgs_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}
        
    if "InFatJet_hasBsubjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "Has b-subjet"
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 1.5}
        
    if "InFatJet_Njettinesstau1" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#tau_{1}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}
        
    if "InFatJet_Njettinesstau2" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#tau_{2}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}
        
    if "InFatJet_Njettinesstau3" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#tau_{3}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}
        
    if "InFatJet_Njettinesstau4" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#tau_{4}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}
        
    if "InFatJet_tau2DIVtau1" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#tau_{2} / #tau_{1}"
        kwargs["cutBox"] = {"cutValue": 0.6, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}
        
    if "InFatJet_tau3DIVtau2" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#tau_{3} / #tau_{2}"
        kwargs["cutBox"] = {"cutValue": 0.67, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0}
        
    if "InFatJet_ldgORsubldg" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}

    if "InFatJet_TFtau21cut_fatjet_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T,jet} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

    if "InFatJet_TFtau32cut_fatjet_pT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T,jet} (%s)"  % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}


        # -----TTsample, top-pt Reweighting ----------------------------------
    if "ttsample_Top_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        # kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0}

  
        # -----TT sample, bs not from the tops --------------------------------
    if "ttsample_bfromtop_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 150.0}
        
    if "ttsample_bNOTfromtop_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 150.0}

    if "ttsample_bfromtop_eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "ttsample_bNOTfromtop_eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "ttsample_cfromtop_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 150.0}

    if "ttsample_cNOTfromtop_pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 150.0}

    if "ttsample_cfromtop_eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "ttsample_cNOTfromtop_eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}
                
    if "ttsample_extra_b_or_c" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "ttsample_intocuts_extra_b_or_c" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}


    if "tt_extra_b_overptCut" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "tt_extra_b_underetaCut" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "tt_extra_b_intobothCuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "tt_extra_c_overptCut" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "tt_extra_c_underetaCut" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "tt_extra_c_intobothCuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}


    # NV on TT sample, bs not from the tops              

    if "ttsample_intocuts_extra_b" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "Extra b "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 6.0}

    if "ttsample_extrab_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "ttsample_Massof_extradib" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "m_{bb} (%s)" % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 250.0}

    if "ttsample_top_extrab_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "ttsample_top_extrab_mindR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "ttsample_bfromtop_extrab_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "ttsample_bfromtop_extrab_mindR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}


    # ---------For higgsMC  b from HIggs as the extra b from above tt analysis------------------
    if "bfromH_as_extra_b_overptCut" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "bfromH_as_extra_b_underetaCut" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "bfromH_as_extra_b_intobothCuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "bfromH_as_extra_b_intocuts" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "bfromH_as_extrab_top_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "bfromH_as_extrab_top_mindR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}
        
    if "bfromH_as_extrab_bfromtop_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

    if "bfromH_as_extrab_bfromtop_mindR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0}

        # -------------soft-b analysis --------------------
        # For M_200, deltaR b from Higgs with other objects from H
    if "bfromH_topfromH_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}
        
    if "bfromH_bfromtopfromH_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}
        
    if "bfromH_objfromW_dR" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR"
        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0}
        
    if "bfromH_Closer_to" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = " "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "Massof_dib" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "m_{bb} (%s)" % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 120.0}


        # Gen particle lvl
    if "bAssociated_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 150.0}
        
    if "bAssociated_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        #kwargs["cutLine"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "softbAssociated_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True} 
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

        
    if "bfromH_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 500.0}
        
    if "bfromH_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "bfromtopfromH_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 500.0}

    if "bfromtopfromH_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}
        
    if "bfromAssociatedTop_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 500.0}
        
    if "bfromAssociatedTop_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

#    if "bOther_Pt" in histo:
#        kwargs["ylabel"] = yLabel + " / %.1f"
#        units            = "GeV/c"
#        kwargs["xlabel"] = "p_{T} (%s)"  % units
#        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
#        kwargs["opts"]   = {"xmin": 0.0, "xmax": 120.0}

#    if "bOther_Eta" in histo:
#        kwargs["ylabel"] = yLabel + " / %.1f"
#        kwargs["xlabel"] = "#eta"
#        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
#        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
#        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}


#    if "softbOther_Eta" in histo:
#        kwargs["ylabel"] = yLabel + " / %.1f"
#        kwargs["xlabel"] = "#eta"
#        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
#        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
#        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

# Gen jet lvl                                                                                                                     
    if "bjetAssociated_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 150.0}

    if "bjetAssociated_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "softbjetAssociated_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}


    if "bjetfromH_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 30.0}

    if "bjetfromH_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "bjetfromtopfromH_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 500.0}

    if "bjetfromtopfromH_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "bjetfromAssociatedTop_Pt" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 500.0}

    if "bjetfromAssociatedTop_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}


        # ---- Histos for dRmin and soft b
    #if "deltaRMin_b_bjet" in histo:
    #    kwargs["ylabel"] = yLabel + " / %.1f"
    #    kwargs["xlabel"] = "#DeltaR_{min}"
    #    kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
    #    kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    #if "deltaRMinb_bjet_proportion_undercut" in histo:
    #    kwargs["ylabel"] = yLabel + " / %.1f"
        #kwargs["xlabel"] = ""
    #    kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
   # if "deltaRMinb_bjet_proportion_undercut_point3" in histo:
   #     kwargs["ylabel"] = yLabel + " / %.1f"
        #kwargs["xlabel"] = ""
   #     kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

   # if "deltaRMin_bassoc_bjet" in histo:                                                 
   #     kwargs["ylabel"] = yLabel + " / %.1f"                                                      
   #     kwargs["xlabel"] = "#DeltaR_{min}"                                          
#        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                 
   #     kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0} 

    if "deltapTMin_bassoc_bjet" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "GeV"
        kwargs["xlabel"] = "#Deltap_{T}_{min} (%s)" % units
#        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 300.0}



    if "associated_b_issoft" in histo:                                                
        kwargs["ylabel"] = yLabel + " / %.1f"                                         
        kwargs["xlabel"] = "p_{T}"                                                                  
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}                                                                                         

    if "soft_associated_b_underetacut" in histo:                                       
        kwargs["ylabel"] = yLabel + " / %.1f"                                                               
        kwargs["xlabel"] = "|#eta|"                                                                                
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}  



    if "soft_nonmatched_associated_b_deltaRMin" in histo:                                                      
        kwargs["ylabel"] = yLabel + " / %.1f"                                                         
        kwargs["xlabel"] = "#DeltaR_{min}"                                                          
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}              
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0} 
        
    if "associated_bjet_issoft" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "p_{T}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "numofmatchedjetswith_assocb" in histo: #after cut deltaR<0.4       
        kwargs["ylabel"] = yLabel + " / %.0f"
        kwargs["xlabel"] = "Number of jets "
       #kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 3.5}

        # ----------------------------- closer look on b from Higgs -----------------------------------------


    if "bfromH_issoft" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "p_{T}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "soft_bfromH_underetacut" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "|#eta|"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}
        
    if "numofmatchedjetswith_bfromH" in histo: #after cut deltaR<0.4                   
        kwargs["ylabel"] = yLabel + " / %.0f"
        kwargs["xlabel"] = "Number of jets "
       #kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 3.5}
        
    if "numofmatchedjetswithpT_bfromH" in histo: #after cut deltaR<0.4                               
        kwargs["ylabel"] = yLabel + " / %.0f"
        kwargs["xlabel"] = "Number of jets "
       #kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                    
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 3.5}

    if "softbjetfromH_Eta" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": -2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["cutBox"] = {"cutValue": 2.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0}

    if "bjetfromH_issoft" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "p_{T}"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 2.0}

    if "soft_nonmatched_bfromH_deltaRMin" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        kwargs["xlabel"] = "#DeltaR_{min}"
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}

    if "soft_nonmatched_bfromH_deltapT" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "ratio "
        kwargs["xlabel"] = "#Deltap_{T}_{min} (%s)" % units
#        kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}                  
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.0}


     # --------------------------------- soft bjet multiplicity- --------------------------------------------- 
    if "SoftBJets" in histo:
        kwargs["ylabel"] = yLabel + " / %.0f"
        kwargs["xlabel"] = "Number of Soft b-jets"
        #kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}          
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 5.5}
        
    if "SoftBJetsMultiplicity" in histo:
        kwargs["ylabel"] = yLabel + " / %.0f"
        kwargs["xlabel"] = "   "
        #kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 5.5}

    if "SoftBJetsCombo" in histo:
        kwargs["ylabel"] = yLabel + " / %.0f"
        kwargs["xlabel"] = "   "
        #kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}              
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 16.5}


        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #if "_Pt" in histo:
    #    kwargs["rebinX"] = 1
    #    units            = "GeV/c"
    #    kwargs["ylabel"] = yLabel + " / %.0f " + units
    #    kwargs["xlabel"] = "p_{T} (%s)"  % units
    #    kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
    #    if "GenJet1" in histo:
    #        kwargs["xlabel"] = "jet_{1} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 1000, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "GenJet2" in histo:
    #        kwargs["xlabel"] = "jet_{2} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 600, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "GenJet3" in histo:
    #        kwargs["xlabel"] = "jet_{3} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 400, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "GenJet4" in histo:
    #        kwargs["xlabel"] = "jet_{4} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 300, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "GenJet5" in histo:
    #        kwargs["xlabel"] = "jet_{5} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 200, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "GenJet6" in histo:
    #        kwargs["xlabel"] = "jet_{6} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 150, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "BQuark1" in histo:
    #        kwargs["xlabel"] = "b_{1} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 600, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "BQuark2" in histo:
    #        kwargs["xlabel"] = "b_{2} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 300, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "BQuark3" in histo:
    #        kwargs["xlabel"] = "b_{3} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 200, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "BQuark4" in histo:
    #        kwargs["xlabel"] = "b_{4} p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 100, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "MaxPt" in histo:
    #        kwargs["log"]    = True
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 1000}#, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "MaxMass" in histo:
    #        kwargs["log"]    = True
    #        ROOT.gStyle.SetNdivisions(5, "X")
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 1000}
    #    if "BQuarkPair_dRMin" in histo:
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 500}#, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "MaxDiJetMass" in histo:
    #        kwargs["rebinX"] = 2
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 600}#, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "MaxTriJetPt" in histo:
    #        units            = "GeV/c"
    #        kwargs["ylabel"] = yLabel + " / %.0f " + units
    #        kwargs["xlabel"] = "p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0.0, "xmax": +1000.0}
    #    if "dRMinDiJet_NoBJets" in histo:
    #        units            = "GeV/c"
    #        kwargs["ylabel"] = yLabel + " / %.0f " + units
    #        kwargs["xlabel"] = "p_{T} (%s)"  % units
    #        kwargs["opts"]   = {"xmin": 0, "xmax": 800, "ymin": 1e+0, "ymaxfactor": 10}
    #    if "BQuarkPair_dR" in histo:
    #        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

#    if "eta" in histo.lower(): 
#        kwargs["ylabel"] = yLabel + " / %.1f"
#        kwargs["xlabel"] = "#eta"
#        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
#        kwargs["opts"]   = {"xmin": -5.0, "xmax": +5.0}#, "ymin": 1e+0, "ymaxfactor": 10}
#
#        if "dEta" in histo:
#            #kwargs["ylabel"] = yLabel + " / %.2f"
#            kwargs["xlabel"] = "#Delta#eta"
#            kwargs["opts"]   = {"xmin": 0.0, "xmax": +5.0}
                
    #if "_Mass" in histo:                        #commented out by kchristo
     #   kwargs["rebinX"] = 2
     #   units            = "GeV/c^{2}"
     #   kwargs["ylabel"] = yLabel + " / %.0f " + units
     #   kwargs["xlabel"] = "M (%s)"  % units
     #   kwargs["log"]    = True
     #   kwargs["opts"]   = {"xmin": 0, "xmax": 800}

    if "BQuarkPair_dRMin_jet2_dPhi" in histo:
        kwargs["ylabel"] = yLabel + " / %.1f"
        units            = "rads"
        kwargs["xlabel"] = "#Delta#phi (%s)" % units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +3.2}

    if opts.normaliseToOne:
        yMin = 1e-5
    else:
        yMin = 1e0

    if kwargs["log"] == True:
        yMaxFactor = 2.0
    else:
        yMaxFactor = 1.2

    # Finalise and return
    kwargs["opts"]["ymaxfactor"] = yMaxFactor
    kwargs["opts"]["ymin"]       = yMin
    return kwargs
  

def PlotMC(datasetsMgr, histo, intLumi):

    kwargs = {}
    if opts.normaliseToOne:
        p = plots.MCPlot(datasetsMgr, histo, normalizeToOne=True, saveFormats=[], **kwargs)
    else:
        p = plots.MCPlot(datasetsMgr, histo, normalizeToLumi=intLumi, saveFormats=[], **kwargs)

    # Get histogram<->kwargs dictionary 
    kwargs = GetHistoKwargs(histo, opts)

    # Customise styling
    if 0:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kSolid))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerSize(0))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerStyle(6))

    # Customise signal
        for d in datasetsMgr.getAllDatasetNames():
            if "Charged" not in d:
                continue
            else:
                #p.histoMgr.setHistoDrawStyle(d, "HIST")
                p.histoMgr.setHistoLegendStyle(d, "L")
                p.histoMgr.setHistoLegendStyle(d, "L")
        
    #  Customise QCD style
    p.histoMgr.setHistoDrawStyle("QCD", "P")
    p.histoMgr.setHistoLegendStyle("QCD", "P")
    p.histoMgr.setHistoLegendLabelMany({
            "QCD": "QCD (MC)",
            })

    # Customise style
    signalM = []
    for m in signalMass:
        signalM.append(m.rsplit("M_")[-1])
    for m in signalM:
        dName = "ChargedHiggs_HplusTB_HplusToTB_M_%s" %m
        if dName in datasetsMgr.getAllDatasetNames():
            p.histoMgr.forHisto(dName, styles.getSignalStyleHToTB_M(m))
            p.histoMgr.setHistoLegendStyle(dName, "LP")

    # Plot customised histogram
    plots.drawPlot(p, 
                   histo,  
                   xlabel       = kwargs.get("xlabel"),
                   ylabel       = kwargs.get("ylabel"),
                   log          = kwargs.get("log"),
                   rebinX       = kwargs.get("rebinX"), 
                   cmsExtraText = "Preliminary", 
                   #createLegend = {"x1": 0.62, "y1": 0.75, "x2": 0.92, "y2": 0.92},
                   moveLegend   = kwargs.get("moveLegend"),
                   opts         = kwargs.get("opts"),
                   opts2        = {"ymin": 0.6, "ymax": 1.4},
                   cutBox       = kwargs.get("cutBox"),
                   #cutLine      = kwargs.get("cutLine") #kchristo
                   )
    
    # Save plot in all formats    
    saveName = histo.split("/")[-1]
    if opts.folder == "":
        savePath = os.path.join(opts.saveDir, opts.optMode)
    else:
        savePath = os.path.join(opts.saveDir, histo.split("/")[0], opts.optMode)
    SavePlot(p, saveName, savePath, [".png", ".pdf"]) 
    return


def SavePlot(plot, saveName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(saveFormats), ", ".join(saveFormats) ) )
    
    # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
        
    savePath = os.path.join(saveDir, saveName)

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = savePath + ext
        saveNameURL = saveNameURL.replace("/publicweb/k/kchristo/", "http://home.fnal.gov/~kchristo/")
        if opts.url:
            Print(saveNameURL, i==0)
        else:
            Print(savePath + ext, i==0)
        plot.saveAs(savePath, formats=saveFormats)
    return


#================================================================================================ 
# Main
#================================================================================================ 
if __name__ == "__main__":
    '''
    https://docs.python.org/3/library/argparse.html
 
    name or flags...: Either a name or a list of option strings, e.g. foo or -f, --foo.
    action..........: The basic type of action to be taken when this argument is encountered at the command line.
    nargs...........: The number of command-line arguments that should be consumed.
    const...........: A constant value required by some action and nargs selections.
    default.........: The value produced if the argument is absent from the command line.
    type............: The type to which the command-line argument should be converted.
    choices.........: A container of the allowable values for the argument.
    required........: Whether or not the command-line option may be omitted (optionals only).
    help............: A brief description of what the argument does.
    metavar.........: A name for the argument in usage messages.
    dest............: The name of the attribute to be added to the object returned by parse_args().
    '''
    
    # Default Settings
    ANALYSISNAME = "kcKinematics"
    SEARCHMODE   = "80to1000"
    DATAERA      = "Run2016"
    OPTMODE      = ""
    BATCHMODE    = True
    PRECISION    = 3
    #SIGNALMASS   = [200, 500, 800, 2000]
    SIGNALMASS   = [300, 500, 800, 1000]
    #SIGNALMASS   = [200, 500, 800, 1000, 2000, 3000]
    INTLUMI      = -1.0
    SUBCOUNTERS  = False
    LATEX        = False
    MERGEEWK     = False
    URL          = False
    NOERROR      = True
    SAVEDIR      = "/publicweb/k/kchristo/" + ANALYSISNAME
    VERBOSE      = False
    HISTOLEVEL   = "Vital" # 'Vital' , 'Informative' , 'Debug'
    NORMALISE    = False
    #FOLDER       = "Inclusive"
    FOLDER       = "TH1"

    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--analysisName", dest="analysisName", type="string", default=ANALYSISNAME,
                      help="Override default analysisName [default: %s]" % ANALYSISNAME)

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--mergeEWK", dest="mergeEWK", action="store_true", default=MERGEEWK, 
                      help="Merge all EWK samples into a single sample called \"EWK\" [default: %s]" % MERGEEWK)

    #parser.add_option("--signalMass", dest="signalMass", type=float, default=SIGNALMASS, 
                      #help="Mass value of signal to use [default: %s]" % SIGNALMASS)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("--histoLevel", dest="histoLevel", action="store", default = HISTOLEVEL,
                      help="Histogram ambient level (default: %s)" % (HISTOLEVEL))

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("-n", "--normaliseToOne", dest="normaliseToOne", action="store_true", 
                      help="Normalise the histograms to one? [default: %s]" % (NORMALISE) )

    parser.add_option("--folder", dest="folder", type="string", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )

    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if opts.mcrab == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        #print __doc__
        sys.exit(1)

    # Sanity check
    if opts.mergeEWK:
        Print("Merging EWK samples into a single Datasets \"EWK\"", True)

    # Sanity check
    allowedMass = [180, 200, 220, 250, 300, 350, 400, 500, 800, 1000] # 2000,3000
    signalMass = []
    for m in sorted(SIGNALMASS, reverse=True):
        signalMass.append("ChargedHiggs_HplusTB_HplusToTB_M_%.f" % m)

    # Call the main function
    main(opts, signalMass)

    if not opts.batchMode:
        raw_input("=== plot_1d.py: Press any key to quit ROOT ...")
