#!/usr/bin/env python
'''
Description:
Script for plotting TH2 plots

Useful Link:
https://nixtricks.wordpress.com/2011/03/03/simple-loops-in-csh-ortcsh/

Usage:
./plot_2d.py -m <pseudo_mcrab_directory> [opts]

Examples:
./plot_2d.py -m BoostedTopologies_StdSelections_TopCut100_AllSelections_NoTrgMatch__H2Cut0p5_NoTopMassCut_170831_085353/ --url -i "M_500"

Last Used:


Obsolete:
foreach x ( 180 200 220 250 300 350 400 500 800 1000 2000 3000 )
./plot_2d.py -m Hplus2tbAnalysis_StdSelections_TopCut100_AllSelections_NoTrgMatch_TopCut10_H2Cut0p5_170826_073257/ -i M_$x --url
end

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
import array
ROOT.gROOT.SetBatch(True)
from ROOT import *

import HiggsAnalysis.NtupleAnalysis.tools.dataset as dataset
import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms
import HiggsAnalysis.NtupleAnalysis.tools.counter as counter
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.styles as styles
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.crosssection as xsect

#================================================================================================ 
# Function Definition
#================================================================================================ 
def GetLumi(datasetsMgr):
    '''
    '''
    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            #print "dataset = %s, lumi = %s" % (d.getName(), lumi)
            lumi += d.getLuminosity()
    #print "Luminosity = %s (pb)" % (lumi), True)
    return lumi

def GetListOfEwkDatasets():
    return  ["TT", "WJetsToQQ_HT_600ToInf", "SingleTop", "DYJetsToQQHT", "TTZToQQ",  "TTWJetsToQQ", "Diboson", "TTTT"]

def GetDatasetsFromDir(opts):
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

def GetAllDatasetsFromDir(opts):
    datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                    dataEra=opts.dataEra,
                                                    searchMode=opts.searchMode, 
                                                    analysisName=opts.analysisName,
                                                    optimizationMode=opts.optMode)
    return datasets
    
def main(opts):

    optModes = [""]
    if opts.optMode != None:
        optModes = [opts.optMode]
        
    # For-loop: All optimisation modes
    for opt in optModes:
        opts.optMode = opt

        # Quick and dirty way to get total int lumi
        allDatasetsMgr = GetAllDatasetsFromDir(opts)
        allDatasetsMgr.updateNAllEventsToPUWeighted()
        allDatasetsMgr.loadLuminosities() # from lumi.json
        plots.mergeRenameReorderForDataMC(allDatasetsMgr) 

        # Determine integrated Lumi before removing data
        # opts.intLumi = GetLumi(allDatasetsMgr) //kchristo                               
        #opts.intLumi = 5747.588 + 2573.399 + 4248.384 + 4008.663 + 2704.118 + 405.222 + 7539.457 + 8390.5 + 215.149

        # Setup & configure the dataset manager 
        datasetsMgr    = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        datasetsMgr.loadLuminosities() # from lumi.json

        # Set/Overwrite cross-sections
        for d in datasetsMgr.getAllDatasets():
            if "ChargedHiggs" in d.getName():
                datasetsMgr.getDataset(d.getName()).setCrossSection(1.0) # ATLAS 13 TeV H->tb exclusion limits
            
        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()

        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 

        # Print dataset information
        datasetsMgr.PrintInfo()

        # Merge EWK samples
        if opts.mergeEWK:
            datasetsMgr.merge("EWK", GetListOfEwkDatasets())
            plots._plotStyles["EWK"] = styles.getAltEWKStyle()

        # Apply TDR style
        style = tdrstyle.TDRStyle()
        style.setOptStat(True)
        style.setWide(True, 0.15)
        # style.setPadRightMargin()#0.13)

#######################################################################################################  
#######################################################################################################      
#here Efficiency plots for pT and tau32 tau21 cuts#####################################################
        datasets = datasetsMgr.getAllDatasets()
        signals = []
        for dataset in datasets:
            if "M_300" in dataset.getName():
                signal300_dataset = dataset
                signals.append(signal300_dataset)
            if "M_200" in dataset.getName():
                signal200_dataset = dataset
                signals.append(signal200_dataset)
            if "M_500" in dataset.getName():
                signal500_dataset = dataset
                signals.append(signal500_dataset)
            if "M_1000" in dataset.getName():
                signal1000_dataset = dataset
                signals.append(signal1000_dataset)
            if "M_800" in dataset.getName():
                signal800_dataset = dataset
                signals.append(signal800_dataset)           

        for hName in [opts.folder+"/"+"Hs_TopProdInFatJet_tau3DIVtau2_Vs_fatjet_pT"]:
            myKwargs = GetHistoKwargs(hName)
            #            for signal_dataset in signals:                 
            GetEfficiencyBeforeTauCut(signals, hName, opts.intLumi, 32, **myKwargs)
            GetEfficiencyAfterTauCut(signals, hName, opts.intLumi, 32, **myKwargs)
            
        for hName in [opts.folder+"/"+"Hs_WProdInFatJet_tau3DIVtau2_Vs_fatjet_pT"]:
            myKwargs = GetHistoKwargs(hName)
            #            for signal_dataset in signals:   
            GetEfficiencyBeforeTauCut(signals, hName, opts.intLumi, 32, **myKwargs)
            GetEfficiencyAfterTauCut(signals, hName, opts.intLumi, 32, **myKwargs)

        for hName in [opts.folder+"/"+"Hs_bqInFatJet_tau3DIVtau2_Vs_fatjet_pT"]:
            myKwargs = GetHistoKwargs(hName)
            #            for signal_dataset in signals:                                              
            GetEfficiencyBeforeTauCut(signals, hName, opts.intLumi, 32, **myKwargs)
            GetEfficiencyAfterTauCut(signals, hName, opts.intLumi, 32, **myKwargs)
            
        for hName in [opts.folder+"/"+"Hs_TopProdInFatJet_tau2DIVtau1_Vs_fatjet_pT"]:
            myKwargs = GetHistoKwargs(hName)
            #            for signal_dataset in signals: 
            GetEfficiencyBeforeTauCut(signals, hName, opts.intLumi, 21, **myKwargs)
            GetEfficiencyAfterTauCut(signals, hName, opts.intLumi, 21, **myKwargs)

        for hName in [opts.folder+"/"+"Hs_WProdInFatJet_tau2DIVtau1_Vs_fatjet_pT"]:
            myKwargs = GetHistoKwargs(hName)
            #            for signal_dataset in signals:       
            GetEfficiencyBeforeTauCut(signals, hName, opts.intLumi, 21, **myKwargs)
            GetEfficiencyAfterTauCut(signals, hName, opts.intLumi, 21, **myKwargs)

        for hName in [opts.folder+"/"+"Hs_bqInFatJet_tau2DIVtau1_Vs_fatjet_pT"]:
            myKwargs = GetHistoKwargs(hName)
            #            for signal_dataset in signals:   
            GetEfficiencyBeforeTauCut(signals, hName, opts.intLumi, 21, **myKwargs)
            GetEfficiencyAfterTauCut(signals, hName, opts.intLumi, 21, **myKwargs)



#######################################################################################################   
#######################################################################################################   
#######################################################################################################


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Do 2D histos
        histoNames  = []
        saveFormats = [".png",".pdf"] #[".C", ".png", ".pdf"]
        histoList   = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(opts.folder)
        histoPaths  = [opts.folder +"/" + h for h in histoList]
        #histoKwargs = GetHistoKwargs(histoPaths)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Axes divisions!
        ROOT.gStyle.SetNdivisions(5, "X")
        ROOT.gStyle.SetNdivisions(5, "Y")

        # For-loop: All histogram
        for histoName in histoPaths:

            myKwargs = GetHistoKwargs(histoName) 
            
            # For-loop: All datasets
            for d in datasetsMgr.getAllDatasetNames():
                #if d != "TT":
                #    continue

                Plot2d(datasetsMgr, d, histoName, myKwargs, opts)
                
                # Avoid replacing canvas with same name warning
                for o in gROOT.GetListOfCanvases():
                    # print o.GetName()
                    o.SetName(o.GetName() + "_" + d)
    return

def GetHistoKwargs(histoName): 
    '''
    Dictionary with 
    key   = histogramName
    value = kwargs
    '''
    
    histoKwargs = {}
    
    # Defaults
    logX        = False
    logY        = False
    logZ        = False
    rebinX      = 1
    rebinY      = 1
    labelX      = None
    labelY      = None
    gridX       = True
    gridY       = True
    _moveLegend = {"dx": -0.1, "dy": -0.01, "dh": 0.1}
    xMin        = 0
    xMax        = 5
    yMin        = 0
    yMax        = 5
    zMin        = 1e-10
    zMax        = 1e+5
    normToOne   = False


    kwargs = {
        "xlabel"           : labelX,
        "ylabel"           : labelY,
        "rebinX"           : 1,
        "rebinY"           : 1,
        "normaliseToOne"   : normToOne, 
        "stackMCHistograms": False,
        "addMCUncertainty" : False,
        "addLuminosityText": False, # you add True if you run on data kchristo
        "addCmsText"       : True,
        "cmsExtraText"     : "Preliminary",
        "logX"             : logX,
        "logY"             : logY,
        "logZ"             : logZ,
        "moveLegend"       : _moveLegend,
        "gridX"            : gridX,
        "gridY"            : gridY,
        "xmin"             : xMin,
        "xmax"             : xMax,
        "ymin"             : yMin,
        "ymax"             : yMax,
        "zmin"             : zMin,
        "zmax"             : zMax,
        "verbose"          : False, #here
        "drawStyle"        : "CPE",  
        "legStyle"         : "LP",  #to here
        "cutBoxY"          : {"cutValue": 0.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True},
        }
    
    h = histoName.replace("TH2/", "")
    kwargs["name"] = h

# kchristo +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# --- Boosted Topologies-----------------------------------------
    if "Iffatjet_Hs_Top_pT_Vs_W_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "Iffatjet_NotHs_Top_pT_Vs_W_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+1
        kwargs["logZ"] = True
    
    elif "Hs_mostdistantfromtop_isb__dRqq_Vs_W_pT" in h:
        kwargs["xmax"] = 2.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "Hs_mostdistantfromtop_isb__dRqq_Vs_top_pT" in h:
        kwargs["xmax"] = 2.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "NotHs_mostdistantfromtop_isb__dRqq_Vs_W_pT" in h:
        kwargs["xmax"] = 2.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "NotHs_mostdistantfromtop_isb__dRqq_Vs_top_pT" in h:
        kwargs["xmax"] = 2.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

        # update 24 /01 
    elif "Hs_BoostedW_W_pT_Vs_Fatjet_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

        # AK8 & n jettiness
    elif "InFatJet_Top_pT_Vs_fatjet_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 3e+1
        kwargs["logZ"] = True
        
    elif "InFatJet_W_pT_Vs_fatjet_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 3e+1
        kwargs["logZ"] = True
        
    elif "InFatJet_tau2DIVtau1_Vs_fatjet_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 3e+1
        kwargs["cutBoxY"] = {"cutValue": 0.6, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["logZ"] = True
        
    elif "InFatJet_tau3DIVtau2_Vs_fatjet_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 3e+1
        kwargs["cutBoxY"] = {"cutValue": 0.67, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["logZ"] = True

    elif "Hs_TopProdInFatJet_tau2DIVtau1_Vs_fatjet_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 3e+1
        kwargs["cutBoxY"] = {"cutValue": 0.6, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["logZ"] = True
        
    elif "InFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 3e+1
        kwargs["logZ"] = True
        
    elif "InFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT" in h:
        kwargs["xmax"] = 1000.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 3e+1
        kwargs["logZ"] = True
        
        # ----- MCHiggs, tops' and their objects dR-------------------------------------------------------------
        
    elif "bfromH_Higgs_dR_Vs_Higgs_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "objectsfromHiggstop_maxdR_Vs_Higgstop_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "objectsfromHiggstop_maxdR_Vs_Higgs_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "objectsfromHiggstop_mindR_Vs_Higgstop_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

    elif "objectsfromHiggstop_mindR_Vs_Higgs_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

    elif "bfromHiggstop_Higgstop_dR_Vs_Higgstop_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True


    elif "WfromHiggstop_closestb_dR_Vs_W_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

    elif "WfromHiggstop_samedecayb_dR_Vs_W_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "WfromHiggstop_closestb_dR_Vs_samedecaytop_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "WfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

    elif "WNOTfromHiggstop_closestb_dR_Vs_W_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

    elif "WNOTfromHiggstop_samedecayb_dR_Vs_W_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

    elif "WNOTfromHiggstop_closestb_dR_Vs_samedecaytop_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

    elif "WNOTfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True


    elif "objectsNOTfromHiggstop_maxdR_Vs_top_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True

    elif "objectsNOTfromHiggstop_mindR_Vs_top_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
    elif "bNOTfromHiggstop_top_dR_Vs_top_Pt" in h:
        kwargs["xmax"] = 5.0
        kwargs["ymax"] = 1000.0
        kwargs["zmin"] = 1e-1
        kwargs["zmax"] = 1e+3
        kwargs["logZ"] = True
        
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
        
    else:
        raise Exception("Unexpected histogram with name \"%s\"" % h)

    if 0:
        for k in kwargs:
            print "%s = %s" % (k, kwargs[k])
        sys.exit()

    return kwargs
    
def GetHisto(datasetsMgr, dataset, histoName):
    h = datasetsMgr.getDataset(dataset).getDatasetRootHisto(histoName)
    return h

def Plot2d(datasetsMgr, dataset, histoName, kwargs, opts):
    '''
    '''
    # Definitions
    saveName = histoName.replace("/", "_")    
    
    #kwargs = {}
    #print opts.normaliseToOne, "+++++++++++++++++++++"
    #if opts.normaliseToOne:
    #    p = plots.MCPlot(datasetsMgr, histoName, normalizeToOne=True , saveFormats=[], **kwargs) 
    #else:
    #    p = plots.MCPlot(datasetsMgr, histoName, normalizeToLumi=opts.intLumi, saveFormats=[], **kwargs)

    
    # Get histogram<->kwargs dictionary                                                                                                              
    #kwargs = GetHistoKwargs(histoName) 

    # Get Histos for the plotter object
    refHisto = GetHisto(datasetsMgr, dataset, histoName)

    # Create the plotting object
    p = plots.PlotBase(datasetRootHistos=[refHisto], saveFormats=kwargs.get("saveFormats"))
    # p.histoMgr.normalizeMCToLuminosity(opts.intLumi) #it was uncomment
    
    # Customise axes (before drawing histo)    
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().RebinX(kwargs.get("rebinX")))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().RebinY(kwargs.get("rebinY")))
    if 0:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetXaxis().SetTitle(kwargs.get("xlabel")))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetYaxis().SetTitle(kwargs.get("ylabel")))
        # p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetXaxis().SetTitle( h.getRootHisto().GetXaxis().GetTitle() + "%0.1f" ))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetXaxis().SetRangeUser(kwargs.get("xmin"), kwargs.get("xmax")))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetYaxis().SetRangeUser(kwargs.get("ymin"), kwargs.get("ymax")))

    #p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetYaxis().SetRangeUser(kwargs.get("ymin"), kwargs.get("ymax")))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetZaxis().SetRangeUser(kwargs.get("zmin"), kwargs.get("zmax")))

    # Create a frame                                                                                                                                                             
    fOpts = {"ymin": 0.0, "ymaxfactor": 1.0}
    p.createFrame(saveName, opts=fOpts)
        
    # SetLog
    SetLogAndGrid(p, **kwargs)

    # Add cut line/box
    if 0:
        _kwargs = { "lessThan": True}
        #p.addCutBoxAndLine(cutValue=400.0, fillColor=ROOT.kGray, box=False, line=True, **_kwargs)
        p.addCutBoxAndLineY(cutValue=0.6, fillColor=ROOT.kGray, box=False, line=True, **_kwargs)

    # Customise Legend
    moveLegend = {"dx": -0.1, "dy": +0.0, "dh": -0.1}
    p.setLegend(histograms.moveLegend(histograms.createLegend(), **kwargs.get("moveLegend")))
    p.removeLegend()

    # Customise histogram (after frame is created)

    #kchristo -----------------------------------------           
    if opts.normaliseToOne:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetZaxis().SetTitle("Arbitrary Units"))
    else:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetZaxis().SetTitle("Events"))
    # -------------------------------------------------- 

    #p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetZaxis().SetTitle("Events"))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetZaxis().SetTitleOffset(1.4))
    #ROOT.gPad.RedrawAxis()

    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMinimum(kwargs.get("zmin")))    
    #p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMaximum(kwargs.get("zmax")))

    # Drawing style    
    p.histoMgr.setHistoDrawStyle(dataset, "COLZ")
    p.histoMgr.setHistoDrawStyleAll("COLZ")

    # The lines below do nothing since the Legend is disabled
    if 0:
        p.histoMgr.setHistoLegendStyle(dataset, "F")
        if "FakeB" in dataset:
            p.histoMgr.setHistoLegendLabelMany({
                    dataset: "QCD (Data)",
                    })
        if dataset == "QCD":
            p.histoMgr.setHistoLegendLabelMany({
                    dataset: "QCD (MC)",
                    })
            
    # Draw the plot
    p.draw()

    # Add canvas text
    if kwargs.get("addLuminosityText"):
        histograms.addStandardTexts(lumi=opts.intLumi)
    else:
        histograms.addStandardTexts()

    # Add dataset name on the canvas
    if dataset == "QCD":
        histograms.addText(0.17, 0.88, "QCD (H_{T} binned)", 27)
    else:
        histograms.addText(0.17, 0.88, plots._legendLabels[dataset], 27)

    # Save the plots
    SavePlot(p, saveName, os.path.join(opts.saveDir, opts.analysisName, opts.folder, dataset, opts.optMode) )
    return

def HasKeys(keyList, **kwargs):
    for key in keyList:
        if key not in kwargs:
            raise Exception("Could not find the keyword \"%s\" in kwargs" % (key) )
    return 

def SetLogAndGrid(p, **kwargs):
    '''
    '''
    HasKeys(["logX", "logY", "logZ", "gridX", "gridY"], **kwargs)
    logX  = kwargs.get("logX")
    logY  = kwargs.get("logY")
    logZ  = kwargs.get("logZ")
    gridX = kwargs.get("gridX")
    gridY = kwargs.get("gridY")

    p.getPad().SetLogx(logX)
    p.getPad().SetLogy(logY)
    p.getPad().SetGridx(gridX)
    p.getPad().SetGridy(gridY)
    if logZ != None:
        p.getPad().SetLogz(logZ)

    return
#######################################################################################################
#######################################################################################################
#here for the Efficiency plots#########################################################################
def set_range(start, end, step):
     while start <= end:
          yield start
          start += step

def GetEfficiencyBeforeTauCut(signals, histoName, intLumi, tauwhat, **mykwargs):
    efficiencyPlots   = []

    taucut = 0
    if tauwhat == 32:
        taucut = 0.67
    elif tauwhat == 21:
        taucut = 0.6

    for signal_dataset in signals:
        HasKeys(["verbose", "drawStyle", "legStyle"], **mykwargs)
        drawStyle    = mykwargs.get("drawStyle")
        legStyle     = mykwargs.get("legStyle")
        legName      = plots._legendLabels[signal_dataset.getName()]

        SignalName = signal_dataset.getName()
        sig = signal_dataset.getDatasetRootHisto(histoName)
        #sig.normalizeToLuminosity(intLumi)       
        h_Signal = sig.getHistogram()

        windowOpeningValues = []

        xeff=[]
        xeffPLACEforDot=[]
        yeff=[]

        my_xrange = []
        my_yrange = []
        my_xstep  = []

        x_bins = h_Signal.GetXaxis().GetNbins()
        y_bins = h_Signal.GetYaxis().GetNbins()

        i = 1
        while i < x_bins+1:
            my_xrange.append(i)
            i = i + 1

        j = 1
        while j < y_bins+1:
            my_yrange.append(j)
            j = j + 1
            
        
        k = 17
        my_xstep.append(k)
        k = 20
        while k < x_bins+1:
            my_xstep.append(k)
            k = k + 5
        
        xmin = h_Signal.GetXaxis().GetBinLowEdge(1)
        xmax = h_Signal.GetXaxis().GetBinUpEdge(x_bins)
        ymin = h_Signal.GetYaxis().GetBinLowEdge(1)
        ymax = h_Signal.GetYaxis().GetBinUpEdge(y_bins)

        #xstepPreviousValue = 0
        
        signal_total = 1e-6 
        
        for xcut in my_xrange:
            for ycut in my_yrange:
                sigContent = h_Signal.GetBinContent(xcut,ycut)
                if sigContent < 0 :
                    sigContent = 0
                signal_total  +=  sigContent
        signal_total -= 1e-6

        for xstep in my_xstep:
            xstepvalue = h_Signal.GetXaxis().GetBinUpEdge(xstep)
            #xstepvalueFORPLOT = xstepvalue - 25
            
            signal_pass = 0
            for xcut in my_xrange:
                #xvalue = h_Signal.GetXaxis().GetBinCenter(xcut)
                xvalue = h_Signal.GetXaxis().GetBinLowEdge(xcut)
                if xvalue < xstepvalue:
                    continue
                for ycut in my_yrange:
                    #    yvalue = h_Signal.GetYaxis().GetBinCenter(ycut)
                        
                        
                    #if(R==1):
                    #    sigContent = h_Signal.GetBinContent(xcut)
                    #    if (sigContent < 0):
                    #        sigContent = 0
                    #    signal_total  +=  sigContent

                        #elif(yvalue < taucut):
                        #    sigContent = h_Signal.GetBinContent(xcut, ycut)
                        #    if (sigContent < 0):
                        #        sigContent = 0
                        #    signal_pass  +=  sigContent
                            
                    #elif(R==2):
                    #    if xvalue < xstepvalue:
                    #        continue
                    #    sigContent = h_Signal.GetBinContent(xcut)
                    #    if (sigContent < 0):
                    #        sigContent = 0
                    #    signal_pass  +=  sigContent
                     
                    sigContent = h_Signal.GetBinContent(xcut,ycut) 
                    if sigContent < 0 :                                                                                                                                                 
                        sigContent = 0
                    
                    #if (xvalue > xstepvalue) :
                    signal_pass  +=  sigContent                        
                
            #xstepPreviousValue = xstepvalue
            efficiency = signal_pass/signal_total
            xeff.append(xstepvalue)
            #xeffPLACEforDot.append(xstepvalueFORPLOT)
            yeff.append(efficiency)
            
            tGraph_eff = ROOT.TGraph(len(xeff), array.array("d",xeff), array.array("d", yeff))
            
            if "M_500" in SignalName:
                styles.signalStyleHToTB500.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 500 GeV"
            elif "M_300" in SignalName:
                styles.signalStyleHToTB300.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 300 GeV"
            elif "M_1000" in SignalName:
                styles.signalStyleHToTB1000.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 1000 GeV"
            elif "M_800" in SignalName:
                styles.signalStyleHToTB800.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 800 GeV"
            elif "M_200" in SignalName:
                styles.signalStyleHToTB200.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 200 GeV"
            elif "TT" in SignalName:
                styles.ttStyle.apply(tGraph_eff)
                legend = "TT"
            else:
                styles.ttStyle.apply(tGraph_eff)
                legend = "other"
                
            if "M_" in SignalName:
                tGraph_eff.SetLineWidth(1)
                tGraph_eff.SetMarkerSize(1)
            else:
                tGraph_eff.SetLineWidth(2)
                tGraph_eff.SetMarkerSize(1)
                

        efficiencyGraph = histograms.HistoGraph(tGraph_eff, legName, "lp", drawStyle)

        efficiencyPlots.append(efficiencyGraph)
        
    if (1):
        #xMin = xmin - 0.1     
        xMin = xmin
        #rebinX = 2  
        xMax = xmax
        xTitle = "p_{T}"
        units = ""
        _format = "%0.1f" + units
        #binwidth = h_Signal.GetXaxis().GetBinWidth(1)          
        #yTitle = "N_{S}/#sqrt{N_{S}+N_{B}} / "   + str(round(binwidth, 2)) + units         
        yTitle_eff = "Efficiency"  + units # + str(round(binwidth, 2)) + units      
        yMin = 0.0
        yMax_eff = 1.0 + 0.1

    options_eff = {"ymin": yMin  , "ymax": yMax_eff, "xmin":xMin, "xMax":xMax}

    SignalName = SignalName.replace("ChargedHiggs_HplusTB_HplusToTB_M_", "_")
    p_eff = plots.PlotBase(datasetRootHistos=efficiencyPlots, saveFormats=mykwargs.get("saveFormats"))

    saveName_eff = "Eff_BeforeTauCut_"+histoName.split("/")[-1]

    p_eff.createFrame(saveName_eff, opts=options_eff)
    p_eff.getFrame().GetXaxis().SetTitle(xTitle)
    p_eff.getFrame().GetYaxis().SetTitle(yTitle_eff)

    # Set range    
    p_eff.getFrame().GetXaxis().SetRangeUser(xMin, xMax)
        
    #moveLegend = {"dx": -0.10, "dy": -0.01, "dh": -0.1}
    moveLegend = {"dx": 0, "dy": 0, "dh": 0}
    p_eff.setLegend(histograms.moveLegend(histograms.createLegend(), **moveLegend))

    # Add Standard Texts to plot                                                                                                                                                                                   
    histograms.addStandardTexts()

    p_eff.draw()

    savePath = os.path.join(opts.saveDir, "HplusMasses", histoName.split("/")[0], opts.optMode)
    SavePlot(p_eff, saveName_eff, savePath)
    return


#also here#####################################################################################
def GetEfficiencyAfterTauCut(signals, histoName, intLumi, tauwhat, **mykwargs):
    efficiencyPlots   = []

    taucut = 0
    if tauwhat == 32:
        taucut = 0.67
    elif tauwhat == 21:
        taucut = 0.6

    for signal_dataset in signals:
        HasKeys(["verbose", "drawStyle", "legStyle"], **mykwargs)
        drawStyle    = mykwargs.get("drawStyle")
        legStyle     = mykwargs.get("legStyle")
        legName      = plots._legendLabels[signal_dataset.getName()]

        SignalName = signal_dataset.getName()
        sig = signal_dataset.getDatasetRootHisto(histoName)
        h_Signal = sig.getHistogram()

        windowOpeningValues = []

        xeff=[]
        xeffPLACEforDot=[]
        yeff=[]
        
        my_xrange = []
        my_yrange = []
        my_xstep  = []

        x_bins = h_Signal.GetXaxis().GetNbins()
        y_bins = h_Signal.GetYaxis().GetNbins()

        i = 1
        while i < x_bins+1:
            my_xrange.append(i)
            i = i + 1

        j = 1
        while j < y_bins+1:
            my_yrange.append(j)
            j = j + 1


        k = 17
        my_xstep.append(k)
        k = 20
        while k < x_bins+1:
            my_xstep.append(k)
            k = k + 5

        xmin = h_Signal.GetXaxis().GetBinLowEdge(1)
        xmax = h_Signal.GetXaxis().GetBinUpEdge(x_bins)
        ymin = h_Signal.GetYaxis().GetBinLowEdge(1)
        ymax = h_Signal.GetYaxis().GetBinUpEdge(y_bins)
        
        signal_total = 1e-6

        for xcut in my_xrange:
            for ycut in my_yrange:
                sigContent = h_Signal.GetBinContent(xcut,ycut)
                if sigContent < 0 :
                    sigContent = 0
                signal_total  +=  sigContent
        signal_total -= 1e-6

        for xstep in my_xstep:
            xstepvalue = h_Signal.GetXaxis().GetBinUpEdge(xstep)
            
            signal_pass = 0
            for xcut in my_xrange:
                xvalue = h_Signal.GetXaxis().GetBinLowEdge(xcut)
                if xvalue < xstepvalue:
                    continue
                for ycut in my_yrange:
                    yvalue = h_Signal.GetYaxis().GetBinLowEdge(ycut)
                    
                    if yvalue > taucut:
                        continue                   

                    sigContent = h_Signal.GetBinContent(xcut,ycut)
                    if sigContent < 0 :
                        sigContent = 0
                    signal_pass  +=  sigContent
                    
            efficiency = signal_pass/signal_total
            xeff.append(xstepvalue)
            yeff.append(efficiency)

            tGraph_eff = ROOT.TGraph(len(xeff), array.array("d",xeff), array.array("d", yeff))

            if "M_500" in SignalName:
                styles.signalStyleHToTB500.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 500 GeV"
            elif "M_300" in SignalName:
                styles.signalStyleHToTB300.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 300 GeV"
            elif "M_1000" in SignalName:
                styles.signalStyleHToTB1000.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 1000 GeV"
            elif "M_800" in SignalName:
                styles.signalStyleHToTB800.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 800 GeV"
            elif "M_200" in SignalName:
                styles.signalStyleHToTB200.apply(tGraph_eff)
                legend = "H^{+} m_{H^{+}} = 200 GeV"
            elif "TT" in SignalName:
                styles.ttStyle.apply(tGraph_eff)
                legend = "TT"
            else:
                styles.ttStyle.apply(tGraph_eff)
                legend = "other"

            if "M_" in SignalName:
                tGraph_eff.SetLineWidth(1)
                tGraph_eff.SetMarkerSize(1)
            else:
                tGraph_eff.SetLineWidth(2)
                tGraph_eff.SetMarkerSize(1)
                
        efficiencyGraph = histograms.HistoGraph(tGraph_eff, legName, "lp", drawStyle)

        efficiencyPlots.append(efficiencyGraph)
        
    if (1):
        #xMin = xmin - 0.1             
        xMin = xmin
        #rebinX = 2            
        xMax = xmax
        xTitle = "p_{T}"
        units = ""
        _format = "%0.1f" + units
        #binwidth = h_Signal.GetXaxis().GetBinWidth(1)      
        #yTitle = "N_{S}/#sqrt{N_{S}+N_{B}} / "   + str(round(binwidth, 2)) + units      
        yTitle_eff = "Efficiency"  + units # + str(round(binwidth, 2)) + units           
        yMin = 0.0
        yMax_eff = 1.0 + 0.1

    options_eff = {"ymin": yMin  , "ymax": yMax_eff, "xmin":xMin, "xMax":xMax}

    SignalName = SignalName.replace("ChargedHiggs_HplusTB_HplusToTB_M_", "_")
    p_eff = plots.PlotBase(datasetRootHistos=efficiencyPlots, saveFormats=mykwargs.get("saveFormats"))

    saveName_eff = "Eff_AfterTauCut_"+histoName.split("/")[-1]
                
    p_eff.createFrame(saveName_eff, opts=options_eff)
    p_eff.getFrame().GetXaxis().SetTitle(xTitle)
    p_eff.getFrame().GetYaxis().SetTitle(yTitle_eff)

    # Set range                   
    p_eff.getFrame().GetXaxis().SetRangeUser(xMin, xMax)
    #moveLegend = {"dx": -0.10, "dy": -0.01, "dh": -0.1}
    moveLegend = {"dx": 0, "dy": 0, "dh": 0}
    p_eff.setLegend(histograms.moveLegend(histograms.createLegend(), **moveLegend))

    # Add Standard Texts to plot  
    histograms.addStandardTexts()

    p_eff.draw()

    savePath = os.path.join(opts.saveDir, "HplusMasses", histoName.split("/")[0], opts.optMode)
    SavePlot(p_eff, saveName_eff, savePath)
    return


########################################################################################################
########################################################################################################
########################################################################################################
        
    #save histos
def SavePlot(plot, plotName, saveDir, saveFormats = [".png",".pdf"]):
    # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = saveName + ext
        saveNameURL = saveNameURL.replace("/publicweb/k/kchristo/", "http://home.fnal.gov/~kchristo/")
        if i==0:
            print "=== plot_2d.py:"
     # print histos
        if opts.url:
            print "\t", saveNameURL
        else:
            print "\t", saveName + ext
        plot.saveAs(saveName, formats=saveFormats)
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
    ANALYSISNAME = "BoostedTopologies"
    SEARCHMODE   = "80to1000"
    DATAERA      = "Run2016"
    OPTMODE      = None
    BATCHMODE    = True
    PRECISION    = 3
    INTLUMI      = -1.0
    SUBCOUNTERS  = False
    LATEX        = False
    MERGEEWK     = False
    URL          = False
    NOERROR      = True
    SAVEDIR      = "/publicweb/k/kchristo/"
    VERBOSE      = False
    HISTOLEVEL   = "Vital" # 'Vital' , 'Informative' , 'Debug' 
    FOLDER       = "TH2"
    NORMALISE    = False

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

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("--histoLevel", dest="histoLevel", action="store", default = HISTOLEVEL,
                      help="Histogram ambient level (default: %s)" % (HISTOLEVEL))

    parser.add_option("--folder", dest="folder", action="store", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("-n", "--normaliseToOne", dest="normaliseToOne", action="store_true",
                      help="Normalise the histograms to one? [default: %s]" % (NORMALISE) )
    #"normalizeToOne", "normalizeByCrossSection", "normalizeToLumi"

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

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plot_2d.py: Press any key to quit ROOT ...")
