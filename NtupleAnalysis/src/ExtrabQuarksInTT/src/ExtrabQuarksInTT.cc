// -*- c++ -*-                                                                                                 
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"
// User                                                                                                            
#include "Auxiliary/interface/Table.h"
#include "Auxiliary/interface/Tools.h"
#include "Tools/interface/MCTools.h"
#include "Tools/interface/DirectionalCut.h"
#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"

//kchristo
#include "DataFormat/interface/AK8Jet.h"
#include "DataFormat/interface/AK8JetsSoftDrop.h"

// ROOT                                                                                                             
#include "TDirectory.h"
#include "Math/VectorUtil.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

typedef Particle<ParticleCollection<double> > genParticle;

struct PtComparator
{
  bool operator() (const genParticle p1, const genParticle p2) const { return ( p1.pt() > p2.pt() ); }
  bool operator() (const math::XYZTLorentzVector p1, const math::XYZTLorentzVector p2) const { return ( p1.pt() > p2.pt() ); }
};


class ExtrabQuarksInTT: public BaseSelector {
public:
  explicit ExtrabQuarksInTT(const ParameterSet& config, const TH1* skimCounters);
  virtual ~ExtrabQuarksInTT() {}

  /// Books histograms                                                                                                   
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree                                          
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event                                                      
  virtual void process(Long64_t entry) override;
  virtual vector<genParticle> GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy=true, const bool hasNoDaughters=false);

  virtual vector<GenJet> GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCut, std::vector<float> etaCut);
  virtual vector<GenJet> GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCut, std::vector<float> etaCut, vector<genParticle> genParticlesToMatch);
  
  virtual TMatrixDSym ComputeMomentumTensor(std::vector<math::XYZTLorentzVector> jets, double r = 2.0);
  virtual TMatrixDSym ComputeMomentumTensor2D(std::vector<math::XYZTLorentzVector> jets);
  virtual vector<float> GetMomentumTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
						     float &C,
						     float &D,
						     float &H2);
  virtual vector<float> GetMomentumTensorEigenValues2D(std::vector<math::XYZTLorentzVector> jets,
						       float &Circularity);
  virtual vector<float> GetSphericityTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
						       float &y23, float &Sphericity, float &SphericityT, float &Aplanarity, float &Planarity, float &Y);
    virtual double GetAlphaT(std::vector<math::XYZTLorentzVector> jets,
			     float &HT,
			     float &JT,
			     float &MHT,
			     float &Centrality);

  //kchristo -----------------------------------------------------------------------------------------------------------------------
  virtual unsigned int GetTheLastCopy(unsigned int firstcopy_index);
  virtual unsigned int GetTheFirstCopy(unsigned int lastcopy_index);
  //--------------------------------------------------------------------------------------------------------------------------------

private:
  // Input parameters 
  const double cfg_Verbose;
  const ParameterSet PSet_ElectronSelection;
  const double cfg_ElectronPtCut;
  const double cfg_ElectronEtaCut;
  const ParameterSet PSet_MuonSelection;
  const double cfg_MuonPtCut;
  const double cfg_MuonEtaCut;
  const ParameterSet PSet_TauSelection;
  const double cfg_TauPtCut;
  const double cfg_TauEtaCut;
  const ParameterSet PSet_JetSelection;
  const std::vector<float> cfg_JetPtCuts;
  const std::vector<float> cfg_JetEtaCuts;
  const DirectionalCut<int> cfg_JetNumberCut;
  const ParameterSet PSet_HtSelection;
  const DirectionalCut<float> cfg_HtCut;
  const ParameterSet PSet_BJetSelection;
  const std::vector<float> cfg_BJetPtCuts;
  const std::vector<float> cfg_BJetEtaCuts;
  const DirectionalCut<int> cfg_BJetNumberCut;
  // METSelection PSet_METSelection;                                                                                                                 
  //TopologySelection PSet_TopologySelection;
  //const DirectionalCut<double> cfg_SphericityCut;
  //const DirectionalCut<double> cfg_AplanarityCut;
  //const DirectionalCut<double> cfg_PlanarityCut;
  //const DirectionalCut<double> cfg_CircularityCut;
  //const DirectionalCut<double> cfg_Y23Cut;
  //const DirectionalCut<double> cfg_CparameterCut;
  //const DirectionalCut<double> cfg_DparameterCut;
  //const DirectionalCut<double> cfg_FoxWolframMomentCut;
  //const DirectionalCut<double> cfg_AlphaTCut;
  //const DirectionalCut<double> cfg_CentralityCut;
  // TopSelection PSet_TopSelection;                                                                        
  const HistogramSettings cfg_PtBinSetting;
  const HistogramSettings cfg_EtaBinSetting;
  const HistogramSettings cfg_PhiBinSetting;
  const HistogramSettings cfg_MassBinSetting;
  const HistogramSettings cfg_DeltaEtaBinSetting;
  const HistogramSettings cfg_DeltaPhiBinSetting;
  const HistogramSettings cfg_DeltaRBinSetting;

  Tools auxTools;

  // Counters                                                                                                                                         
  Count cAllEvents;
  Count cTrigger;
  Count cElectronVeto;
  Count cMuonVeto;
  Count cTauVeto;
  Count cJetSelection;
  Count cBJetSelection;
  //Count cTopologySelection;
  Count cTopSelection;
  Count cSelected;

  // Event Variables   
  WrappedTH1 *h_genMET_Et;
  WrappedTH1 *h_genMET_Phi;
  WrappedTH1 *h_genHT_GenJets;

  // Event-Shape Variables                                                                                                                            
  //WrappedTH1 *h_y23;
  //WrappedTH1 *h_Sphericity;
  //WrappedTH1 *h_SphericityT;
  //WrappedTH1 *h_Y;
  //WrappedTH2 *h_S_Vs_Y;
  //WrappedTH1 *h_Aplanarity;
  //WrappedTH1 *h_Planarity;
  //WrappedTH1 *h_CParameter;
  //WrappedTH1 *h_DParameter;
  //WrappedTH1 *h_H2;
  //WrappedTH1 *h_Circularity;
  //WrappedTH1 *h_Centrality;
  //WrappedTH1 *h_HT;
  //WrappedTH1 *h_JT;
  //WrappedTH1 *h_MHT;
  //WrappedTH1 *h_AlphaT;

  // GenParticles: BQuarks                                                                                                                            
  //  vector<WrappedTH1*> vh_BQuarks_Eta;                                                                                                             

  // GenParticles: BQuarks pair closest together                                                                                                      

  // GenJets                                                                                                                                          
  WrappedTH1 *h_GenJets_N;
  
  //kchristo ---------------------------------------------------------------------------------------                                                                                                              
  // --OV----TT sample, bs not from the tops---------------
  WrappedTH1 *h_ttsample_bfromtop_pt;
  WrappedTH1 *h_ttsample_bNOTfromtop_pt;
  WrappedTH1 *h_ttsample_bfromtop_eta;
  WrappedTH1 *h_ttsample_bNOTfromtop_eta;

  WrappedTH1 *h_ttsample_cfromtop_pt;
  WrappedTH1 *h_ttsample_cNOTfromtop_pt;
  WrappedTH1 *h_ttsample_cfromtop_eta;
  WrappedTH1 *h_ttsample_cNOTfromtop_eta;

  WrappedTH1 *h_ttsample_extra_b_or_c;
  WrappedTH1 *h_ttsample_intocuts_extra_b_or_c;

  WrappedTH1 *h_tt_extra_b_overptCut;
  WrappedTH1 *h_tt_extra_b_underetaCut;
  WrappedTH1 *h_tt_extra_b_intobothCuts;

  WrappedTH1 *h_tt_extra_c_overptCut;
  WrappedTH1 *h_tt_extra_c_underetaCut;
  WrappedTH1 *h_tt_extra_c_intobothCuts;

  // --NV----TT sample, bs not from the tops---------------
  WrappedTH1 *h_ttsample_intocuts_extra_b;
  WrappedTH1 *h_ttsample_extrab_dR;
  WrappedTH1 *h_ttsample_Massof_extradib;

  WrappedTH1 *h_ttsample_top_extrab_dR;
  WrappedTH1 *h_ttsample_top_extrab_mindR;
  WrappedTH1 *h_ttsample_bfromtop_extrab_dR;
  WrappedTH1 *h_ttsample_bfromtop_extrab_mindR;

  WrappedTH2 *h_ttsample_ldgtop_extrab_dR_Vs_subldgtop_extrab_dR;


  // ----For higgsMC  b from HIggs as the extra b from above tt analysis---
  WrappedTH1 *h_bfromH_as_extra_b_overptCut;
  WrappedTH1 *h_bfromH_as_extra_b_underetaCut;
  WrappedTH1 *h_bfromH_as_extra_b_intobothCuts;

  WrappedTH1 *h_bfromH_as_extra_b_intocuts;
  WrappedTH1 *h_bfromH_as_extrab_top_dR;
  WrappedTH1 *h_bfromH_as_extrab_top_mindR;
  WrappedTH1 *h_bfromH_as_extrab_bfromtop_dR;
  WrappedTH1 *h_bfromH_as_extrab_bfromtop_mindR;

  WrappedTH2 *h_bfromH_as_extrab_ldgtop_dR_Vs_subldgtop_extrab_dR;
 
  // ---------------------------------------------------------------------------------------------------
                                                  
};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(ExtrabQuarksInTT);

ExtrabQuarksInTT::ExtrabQuarksInTT(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    cfg_Verbose(config.getParameter<bool>("verbose")),
    PSet_ElectronSelection(config.getParameter<ParameterSet>("ElectronSelection")),
    cfg_ElectronPtCut(config.getParameter<float>("ElectronSelection.electronPtCut")),
    cfg_ElectronEtaCut(config.getParameter<float>("ElectronSelection.electronEtaCut")),
    PSet_MuonSelection(config.getParameter<ParameterSet>("MuonSelection")),
    cfg_MuonPtCut(config.getParameter<float>("MuonSelection.muonPtCut")),
    cfg_MuonEtaCut(config.getParameter<float>("MuonSelection.muonEtaCut")),
    PSet_TauSelection(config.getParameter<ParameterSet>("TauSelection")),
    cfg_TauPtCut(config.getParameter<float>("TauSelection.tauPtCut")),
    cfg_TauEtaCut(config.getParameter<float>("TauSelection.tauEtaCut")),
    PSet_JetSelection(config.getParameter<ParameterSet>("JetSelection")),
    cfg_JetPtCuts(config.getParameter<std::vector<float>>("JetSelection.jetPtCuts")),
    cfg_JetEtaCuts(config.getParameter<std::vector<float>>("JetSelection.jetEtaCuts")),
    cfg_JetNumberCut(config, "JetSelection.numberOfJetsCut"),
    PSet_HtSelection(config.getParameter<ParameterSet>("JetSelection")),
    cfg_HtCut(config, "JetSelection.HTCut"),
    PSet_BJetSelection(config.getParameter<ParameterSet>("BJetSelection")),
    cfg_BJetPtCuts(config.getParameter<std::vector<float>>("BJetSelection.jetPtCuts")),
    cfg_BJetEtaCuts(config.getParameter<std::vector<float>>("BJetSelection.jetEtaCuts")),
    cfg_BJetNumberCut(config, "BJetSelection.numberOfBJetsCut"),
    // PSet_METSelection(config.getParameter<ParameterSet>("METSelection")),                        
    //PSet_TopologySelection(config.getParameter<ParameterSet>("TopologySelection")),
    //cfg_SphericityCut(config, "TopologySelection.SphericityCut"),
    //cfg_AplanarityCut(config, "TopologySelection.AplanarityCut"),
    //cfg_PlanarityCut(config, "TopologySelection.PlanarityCut"),
    //cfg_CircularityCut(config, "TopologySelection.CircularityCut"),
    //cfg_Y23Cut(config, "TopologySelection.Y23Cut"),
    //cfg_CparameterCut(config, "TopologySelection.CparameterCut"),
    //cfg_DparameterCut(config, "TopologySelection.DparameterCut"),
    //cfg_FoxWolframMomentCut(config, "TopologySelection.FoxWolframMomentCut"),
    //cfg_AlphaTCut(config, "TopologySelection.AlphaTCut"),
    //cfg_CentralityCut(config, "TopologySelection.CentralityCut"),
    // PSet_TopSelection(config.getParameter<ParameterSet>("TopSelection")),                                                                         
    cfg_PtBinSetting(config.getParameter<ParameterSet>("CommonPlots.ptBins")),
    cfg_EtaBinSetting(config.getParameter<ParameterSet>("CommonPlots.etaBins")),
    cfg_PhiBinSetting(config.getParameter<ParameterSet>("CommonPlots.phiBins")),
    cfg_MassBinSetting(config.getParameter<ParameterSet>("CommonPlots.invMassBins")),
    cfg_DeltaEtaBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaEtaBins")),
    cfg_DeltaPhiBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaPhiBins")),
    cfg_DeltaRBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaRBins")),
    cAllEvents(fEventCounter.addCounter("All events")),
    cTrigger(fEventCounter.addCounter("Trigger")),
    cElectronVeto(fEventCounter.addCounter("e-veto")),
    cMuonVeto(fEventCounter.addCounter("#mu-veto")),
    cTauVeto(fEventCounter.addCounter("#tau-veto")),
    cJetSelection(fEventCounter.addCounter("Jets + H_{T}")),
    cBJetSelection(fEventCounter.addCounter("b-jets")),
    //cTopologySelection(fEventCounter.addCounter("Topology")),
    cTopSelection(fEventCounter.addCounter("Top")),
    cSelected(fEventCounter.addCounter("All Selections"))
{ }

void ExtrabQuarksInTT::book(TDirectory *dir) {

  // Fixed binning                                                                                                                                   
  const int nBinsPt   = 4*cfg_PtBinSetting.bins();
  const double minPt  = cfg_PtBinSetting.min();
  const double maxPt  = 4*cfg_PtBinSetting.max();

  //const int nBinsEta  = 2*cfg_EtaBinSetting.bins();                                                                                               
  //  const double minEta = 2*cfg_EtaBinSetting.min();                                                                                              
  //  const double maxEta = 2*cfg_EtaBinSetting.max();                                                                                               

  //  const int nBinsRap  = cfg_EtaBinSetting.bins();                                   
  //  const double minRap = cfg_EtaBinSetting.min();                                                                                                
  //  const double maxRap = cfg_EtaBinSetting.max();                                                                                               
 
  const int nBinsPhi  = cfg_PhiBinSetting.bins();
  const double minPhi = cfg_PhiBinSetting.min();
  const double maxPhi = cfg_PhiBinSetting.max();

  const int nBinsM  = cfg_MassBinSetting.bins();
  const double minM = cfg_MassBinSetting.min();
  const double maxM = cfg_MassBinSetting.max();

  //const int nBinsdEta  = 2*cfg_DeltaEtaBinSetting.bins();                                                                  
  //  const double mindEta = cfg_DeltaEtaBinSetting.min();                                                                                           
  //  const double maxdEta = cfg_DeltaEtaBinSetting.max();                                                                              

  //  const int nBinsdRap  = 2*cfg_DeltaEtaBinSetting.bins();                                                                                        
  //  const double mindRap = cfg_DeltaEtaBinSetting.min();                                                                                         
  //  const double maxdRap = cfg_DeltaEtaBinSetting.max();                                                                                          

  //  const int nBinsdPhi  = 2*cfg_DeltaPhiBinSetting.bins();                                                                                        
  //  const double mindPhi = cfg_DeltaPhiBinSetting.min();                                                                            
  //  const double maxdPhi = cfg_DeltaPhiBinSetting.max();                                                                         

  const int nBinsdR  = 2*cfg_DeltaRBinSetting.bins();
  const double mindR = cfg_DeltaRBinSetting.min();
  const double maxdR = cfg_DeltaRBinSetting.max();

  TDirectory* th1 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH1");
  TDirectory* th2 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH2");

  // kchristo/////////////////////////////////////////////////////////////////////////////////////////////
  std::string myInclusiveLabel = "Inclusive";
  std::string myTrueLabel      = "BoolIsTrue";
  std::string myFalseLabel     = "BoolIsFalse";
  
  // Create directories  
  TDirectory* myInclusiveDir    = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myInclusiveLabel);
  TDirectory* myBoolIsFalseDir  = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myFalseLabel);
  TDirectory* myBoolIsTrueDir   = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myTrueLabel);
  std::vector<TDirectory*> myTripletDirs = {myInclusiveDir, myBoolIsFalseDir, myBoolIsTrueDir};
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Event Variables                                                                                      
  h_genMET_Et         =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital      , th1, "genMET_Et"    , ";Gen E_{T}^{miss} (GeV)"       , 100,  0.0,   +500.0);
  h_genMET_Phi        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "genMET_Phi"   , ";Gen E_{T}^{miss} #phi (rads)" , nBinsPhi, minPhi, maxPhi);
  h_genHT_GenJets     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital      , th1, "genHT_GenJets", ";GenJ H_{T} (GeV)"             , nBinsM  , minM  , maxM   );

  // Event-Shape Variables                                  
  //h_y23         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "y23"        , ";y_{23}"        , 25, 0.0,    0.25);
  //h_Sphericity  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Sphericity" , ";Sphericity"    , 20, 0.0,    1.00);
  //h_SphericityT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "SphericityT", ";Sphericity_{T}", 20, 0.0,    1.00);
  //h_Y           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Y"          , ";Y"             , 50, 0.0,    0.50);
  //h_S_Vs_Y      = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "S_Vs_Y"     , ";Sphericity;Y=#frac{#sqrt{3}}{2}x(Q1-Q2)", 100, 0.0, 1.0, 50, 0.0, 0.5);
  //h_Aplanarity  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Aplanarity" , ";Aplanarity" , 25, 0.0, 0.5);
  //h_Planarity   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Planarity"  , ";Planarity"  , 25, 0.0, 0.5);
  //h_CParameter  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "CParameter" , ";C"          , 20, 0.0, 1.0);
  //h_DParameter  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "DParameter" , ";D"          , 20, 0.0, 1.0);
  //h_H2          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "H2"         , ";H_{2}"      , 20, 0.0, 1.0);
  //h_Circularity = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Circularity", ";Circularity", 20, 0.0, 1.0);
  //h_Centrality  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Centrality" , ";Centrality" , 20, 0.0, 1.0);
  //h_HT          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HT"         , ";H_{T}"      , 30, 0.0, 4000.0);
  //h_JT          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "JT"         , ";J_{T}"      , 30, 0.0, 4000.0);
  //h_MHT         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MHT"        , ";MHT"        , 50, 0.0,  500.0);
  //h_AlphaT      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "AlphaT"     , ";#alpha_{T}" , 20, 0.0,    1.0);

  // GenParticles: B-quarks                                                                                                                        

  // GenParticles: BQuarks pairs                                                                        

  // GenParticles: BQuarks pairs with maximum pT                                          

  // GenParticles: BQuarks pairs with maximum mass                        

  // GenParticles: BQuarks pair closest together                    

  // Leading Jets                                                                                 

  // GenJets                                                                                       
  h_GenJets_N   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJets_N" , ";genJet multiplicity" , 30, 0.0, 30.0);

  // GenJets: Dijet with largest mass                                                                                                                

  // GenJets: Untagged jet pair with min dR                                                                                                          

  // GenJets: Trijet with largest pT                                                                                                                

  // Correlations                                           

  //kchristo----------------------------------------------------------------------------------------------------------------------------         
  
  // --OV-------TT sample, bs not from the tops---------------    
  h_ttsample_bfromtop_pt       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_bfromtop_pt",";p_{T} (GeV)",30, 0.0, 150.0);
  h_ttsample_bNOTfromtop_pt    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_bNOTfromtop_pt",";p_{T} (GeV)",30, 0.0, 150.0);
  h_ttsample_bfromtop_eta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_bfromtop_eta",";#eta",100 , -5.0, 5.0);
  h_ttsample_bNOTfromtop_eta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_bNOTfromtop_eta",";#eta",100 , -5.0, 5.0);

  h_ttsample_cfromtop_pt       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_cfromtop_pt",";p_{T} (GeV)",30, 0.0, 150.0);
  h_ttsample_cNOTfromtop_pt    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_cNOTfromtop_pt",";p_{T} (GeV)",30, 0.0, 150.0);
  h_ttsample_cfromtop_eta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_cfromtop_eta",";#eta",100 , -5.0, 5.0);
  h_ttsample_cNOTfromtop_eta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_cNOTfromtop_eta",";#eta",100 , -5.0, 5.0);
  h_ttsample_extra_b_or_c      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_extra_b_or_c"," ",5 , 0.0, 5.0);
  h_ttsample_intocuts_extra_b_or_c      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_intocuts_extra_b_or_c"," ",5 , 0.0, 5.0);

  h_tt_extra_b_overptCut       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "tt_extra_b_overptCut"," ", 2 , 0.0, 2.0);
  h_tt_extra_b_underetaCut     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "tt_extra_b_underetaCut"," ", 2 , 0.0, 2.0);
  h_tt_extra_b_intobothCuts    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "tt_extra_b_intobothCuts"," ", 2 , 0.0, 2.0);

  h_tt_extra_c_overptCut       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "tt_extra_c_overptCut"," ", 2 , 0.0, 2.0);
  h_tt_extra_c_underetaCut     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "tt_extra_c_underetaCut"," ", 2 , 0.0, 2.0);
  h_tt_extra_c_intobothCuts    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "tt_extra_c_intobothCuts"," ", 2 , 0.0, 2.0);

  // --NV-----TT sample, bs not from the tops---------------
  h_ttsample_intocuts_extra_b = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_intocuts_extra_b","Extra b",6 , 0.0, 6.0);
  h_ttsample_extrab_dR        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_extrab_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_ttsample_Massof_extradib  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_Massof_extradib",";m_{bb} (GeV)", 50, 0.0, 250.0);
  
  h_ttsample_top_extrab_dR         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_top_extrab_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_ttsample_top_extrab_mindR      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_top_extrab_mindR", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );
  h_ttsample_bfromtop_extrab_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_bfromtop_extrab_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_ttsample_bfromtop_extrab_mindR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "ttsample_bfromtop_extrab_mindR", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );

  h_ttsample_ldgtop_extrab_dR_Vs_subldgtop_extrab_dR = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "ttsample_ldgtop_extrab_dR_Vs_subldgtop_extrab_dR", ";#DeltaR_{ldgtop,b};#DeltaR_{subldgtop,b}", nBinsdR, 0.0, 5.0, nBinsdR, 0.0, 5.0);

  //---For higgsMC  b from HIggs as the extra b from above tt analysis-----
  h_bfromH_as_extra_b_overptCut    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_as_extra_b_overptCut"," ", 2 , 0.0, 2.0);
  h_bfromH_as_extra_b_underetaCut  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_as_extra_b_underetaCut"," ", 2 , 0.0, 2.0);
  h_bfromH_as_extra_b_intobothCuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_as_extra_b_intobothCuts"," ", 2 , 0.0, 2.0);
  h_bfromH_as_extra_b_intocuts     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_as_extra_b_intocuts"," ", 2 , 0.0, 2.0);
  
  h_bfromH_as_extrab_top_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_as_extrab_top_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_bfromH_as_extrab_top_mindR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_as_extrab_top_mindR", ";#DeltaR_{min}", 50, 0.0 , 5.0 );
  h_bfromH_as_extrab_bfromtop_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_as_extrab_bfromtop_dR", ";#DeltaR", 50 ,0.0 , 5.0 );
  h_bfromH_as_extrab_bfromtop_mindR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_as_extrab_bfromtop_mindR", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );

  h_bfromH_as_extrab_ldgtop_dR_Vs_subldgtop_extrab_dR = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bfromH_as_extrab_ldgtop_dR_Vs_subldgtop_extrab_dR", ";#DeltaR_{ldgtop,b};#DeltaR_{subldgtop,b}", nBinsdR, 0.0, 5.0, nBinsdR, 0.0, 5.0);

    // -----------------------------------------------------------------------------------------------------------------------------------
  // end of plots

  return;
}

void ExtrabQuarksInTT::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
}


void ExtrabQuarksInTT::process(Long64_t entry) {

  if ( !fEvent.isMC() ) return;

  // Increment Counter                                                                                               
  cAllEvents.increment();

  // Initialise MCTools object                 
  MCTools mcTools(fEvent);

  //================================================================================================    
  // 1) Apply trigger                                                                             
  //================================================================================================          
  if (cfg_Verbose) std::cout << "=== Trigger" << std::endl;
  if ( !(fEvent.passTriggerDecision()) ) return;
  cTrigger.increment();


  //================================================================================================              
  // 2) MET filters (to remove events with spurious sources of fake MET)                                      
  //================================================================================================         
  // nothing to do                                                                                                               

  //================================================================================================    
  // 3) Primarty Vertex (Check that a PV exists)                                                                  
  //================================================================================================                
  // nothing to do                                                                                                                                   

  //================================================================================================                      
  // 4) Electron veto (fully hadronic + orthogonality)                                                                      
  //================================================================================================                          
  if (cfg_Verbose) std::cout << "=== Electron veto" << std::endl;
  vector<genParticle> selectedElectrons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_ElectronPtCut, cfg_ElectronEtaCut, 11, true, false);
  if (0)
    {
      std::cout << "\nnElectrons = " << selectedElectrons.size() << std::endl;
      for (auto& p: selectedElectrons) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
    }
  if ( selectedElectrons.size() > 0 ) return;
  cElectronVeto.increment();


  //================================================================================================
  // 5) Muon veto (fully hadronic + orthogonality)                                                                                
  //================================================================================================                   
  if (cfg_Verbose) std::cout << "=== Muon veto" << std::endl;
  vector<genParticle> selectedMuons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_MuonPtCut, cfg_MuonEtaCut, 13, true, false);
  if (cfg_Verbose)
    {
      std::cout << "\nnMuons = " << selectedMuons.size() << std::endl;
      for (auto& p: selectedMuons) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
    }
  if ( selectedMuons.size() > 0 ) return;
  cMuonVeto.increment();


  //================================================================================================            
  // 6) Tau veto (HToTauNu orthogonality)                                                                              
  //================================================================================================       
  if (cfg_Verbose) std::cout << "=== Tau veto" << std::endl;
  vector<genParticle> selectedTaus = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_TauPtCut, cfg_TauEtaCut, 15, true, false);
  if (0)
    {
      std::cout << "\nnTaus = " << selectedTaus.size() << std::endl;
      for (auto& p: selectedTaus) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
      std::cout << "" << std::endl;
    }
  if ( selectedTaus.size() > 0 ) return;
  cTauVeto.increment();


  //================================================================================================           
  // 7) Jet Selection                                                                                        
  //================================================================================================                   
  if (cfg_Verbose) std::cout << "=== Jet Selection" << std::endl;
  vector<GenJet> selectedJets = GetGenJets(fEvent.genjets(), cfg_JetPtCuts, cfg_JetEtaCuts);
  // for (auto& p: selectedJets) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;            
  if (!cfg_JetNumberCut.passedCut(selectedJets.size())) return;

  // HT Selection                                                                                              
  double genJ_HT = 0.0;
  std::vector<math::XYZTLorentzVector> selJets_p4;
  math::XYZTLorentzVector jet_p4, bjet_p4 ; // I put bjet_p4 here for the p4 of the B jets
  for(auto jet: selectedJets)
    {
      jet_p4 = jet.p4();
      genJ_HT += jet.pt();
      selJets_p4.push_back( jet_p4 );
    }

  if ( !cfg_HtCut.passedCut(genJ_HT) ) return;
  cJetSelection.increment();

  //================================================================================================                                                  
  // 8) BJet Selection                                                                                                                                
  //================================================================================================                                                  
  if (cfg_Verbose) std::cout << "=== BJet Selection" << std::endl;
  vector<genParticle> selectedBQuarks = GetGenParticles(fEvent.genparticles().getGenParticles(), 10, 3, 5, true, false);
  std::sort( selectedBQuarks.begin(), selectedBQuarks.end(), PtComparator() );
  if (0) for (auto& p: selectedBQuarks) mcTools.PrintGenParticle(p);

  // Match b-quarks with GenJets                                                                                                                      
  vector<GenJet> selectedBJets = GetGenJets(selectedJets, cfg_BJetPtCuts, cfg_BJetEtaCuts, selectedBQuarks);
  // for (auto& p: selectedBJets) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;                                   
  
  // kchristo ----------------------------------------------------------------------------------------
  // useful to a vector with the p4 of the B jets
  
  std::vector<math::XYZTLorentzVector> selectedBJets_p4;
  for (auto& bjet: selectedBJets)
    {
      bjet_p4 = bjet.p4();
      selectedBJets_p4.push_back(bjet_p4);
    }
  
  // -------------------------------------------------------------------------------------------------

  // Get selected jets excluding the matched bjets                                       
  bool isBJet = false;
  std::vector<math::XYZTLorentzVector> selJets_NoBJets_p4;

  // For-loop: Selected jets                                                   
  for (auto& jet: selectedJets)
    {
      isBJet = false;

      // For-loop: Selected bjets                                                                                     
      for (auto& bjet: selectedBJets)
        {
          double dR = ROOT::Math::VectorUtil::DeltaR(jet.p4(), bjet.p4());
          if (dR < 0.01) isBJet = true;
        }
      if (isBJet) continue;
      jet_p4 = jet.p4();
      selJets_NoBJets_p4.push_back(jet_p4);
    }

  if (!cfg_BJetNumberCut.passedCut(selectedBJets.size())) return;
  cBJetSelection.increment();


  //================================================================================================                           
  // 9) All Soft Jets                                                                                                                               
  //================================================================================================                                                

  //================================================================================================                                  
  // 10) MET selection                                                                                                                               
  //================================================================================================                                                 
  // nothing to do                                                                                                                                 

  //================================================================================================                    
  // 11) Topology selection                                                                                  
  //================================================================================================                                  
  //float C, D, H2;
  //float Circularity;
  //float y23, Sphericity, SphericityT, Aplanarity, Planarity, Y; // functions to return values when properly implemented           
  //float HT, JT, MHT, Centrality;
  //vector<float> a = GetMomentumTensorEigenValues(selJets_p4, C, D, H2);
  //vector<float> b = GetMomentumTensorEigenValues2D(selJets_p4, Circularity);
  //vector<float> c = GetSphericityTensorEigenValues(selJets_p4, y23, Sphericity, SphericityT, Aplanarity, Planarity, Y);
  //double alphaT   = GetAlphaT(selJets_p4, HT, JT, MHT, Centrality);

  // Apply cuts                                                                                           
  //if ( !cfg_CparameterCut.passedCut(C) ) return;
  //if ( !cfg_DparameterCut.passedCut(D) ) return;
  //if ( !cfg_FoxWolframMomentCut.passedCut(H2) ) return;
  //if ( !cfg_CircularityCut.passedCut(Circularity) ) return;
  //if ( !cfg_Y23Cut.passedCut(y23) ) return;
  //if ( !cfg_SphericityCut.passedCut(Sphericity) ) return;
  //if ( !cfg_AplanarityCut.passedCut(Aplanarity) ) return;
  //if ( !cfg_PlanarityCut.passedCut(Planarity) ) return;
  //if ( !cfg_CentralityCut.passedCut(Centrality) ) return;
  //if ( !cfg_AlphaTCut.passedCut(alphaT) ) return;
  //cTopologySelection.increment();


  //================================================================================================                                    
  // 12) Top selection                                                                                                                             
  //================================================================================================      
  cTopSelection.increment();


  //================================================================================================                    
  // Standard Selections                                                                                                         
  //================================================================================================                                                
 

  //================================================================================================                                                 
  // All Selections                                                                                                       
  //================================================================================================            
  cSelected.increment();
  //HERE                                                                                         
  // Fill Histograms     
  h_GenJets_N    ->Fill(selectedJets.size());
  h_genMET_Et    ->Fill(fEvent.genMET().et());
  h_genMET_Phi   ->Fill(fEvent.genMET().Phi());
  h_genHT_GenJets->Fill(genJ_HT);
  //h_y23          ->Fill(y23);
  //h_Sphericity   ->Fill(Sphericity);
  //h_SphericityT  ->Fill(SphericityT);
  //h_Y            ->Fill(Y);
  //h_S_Vs_Y       ->Fill(Sphericity, Y);
  //h_Aplanarity   ->Fill(Aplanarity);
  //h_Planarity    ->Fill(Planarity);
  //h_CParameter   ->Fill(C);
  //h_DParameter   ->Fill(D);
  //h_H2           ->Fill(H2);
  //h_Circularity  ->Fill(Circularity);
  //h_Centrality   ->Fill(Centrality);
  //h_HT           ->Fill(HT);
  //h_JT           ->Fill(JT);
  //h_MHT          ->Fill(MHT);
  //h_AlphaT       ->Fill(alphaT);

  ///////////////////////////////////////////////////////////////////////////                                            
  // GenParticles                                                                                                           
  ///////////////////////////////////////////////////////////////////////////                                      
  if (0) std::cout << "=== GenParticles" << std::endl;
  std::vector<math::XYZTLorentzVector> bQuarks_p4;
  math::XYZTLorentzVector tbWPlus_BQuark_p4;
  math::XYZTLorentzVector tbWPlus_Wqq_Quark_p4;
  math::XYZTLorentzVector tbWPlus_Wqq_AntiQuark_p4;
  math::XYZTLorentzVector tbWMinus_BQuark_p4;
  math::XYZTLorentzVector tbWMinus_Wqq_Quark_p4;
  math::XYZTLorentzVector tbWMinus_Wqq_AntiQuark_p4;

  // Define the table                                                                                                                                
  Table table("Evt | Index | PdgId | Status | Charge | Pt | Eta | Phi | E | Vertex (mm) | Lxy (mm) | d0 (mm) | Mothers | Daughters |", "Text"); //LaTeX or Text                                                                                                            
  int row = 0;

// For-loop: GenParticles                                                                                                                           
 for (auto& p: fEvent.genparticles().getGenParticles()) {

   //    mcTools.PrintGenParticle(p);                                                                                                                
   // mcTools.PrintGenDaughters(p);                                                                                                                  

   // Particle properties                                                                                                                            
   short genP_index   = p.index();
   int genP_pdgId     = p.pdgId();
   int genP_status    = p.status();
   double genP_pt     = p.pt();
   double genP_eta    = p.eta();
   double genP_phi    = p.phi();
   double genP_energy = p.e();
   int genP_charge    = p.charge();
   //ROOT::Math::XYZPoint genP_vtx(p.vtxX()*10, p.vtxY()*10, p.vtxZ()*10); // in mm                                                                 
   math::XYZTLorentzVector genP_p4;
   genP_p4 = p.p4();

   // Get vectors for mom/daus                                                                                                                       
   std::vector<short> genP_mothers   = p.mothers();
   std::vector<short> genP_daughters = p.daughters();

   // Assign mother/daughers for Lxy & d0 calculation                                                                                                
   // genParticle m;                                                                                                                                 
   // genParticle d;                                                                                                                                 
   // if (genP_mothers.size() > 0  ) m = fEvent.genparticles().getGenParticles()[genP_mothers.at(0)];                                                
   // if (genP_daughters.size() > 0) d = fEvent.genparticles().getGenParticles()[genP_daughters.at(0)];                                              

   // Filtering                                                                                                                                      
   //if (!p.isLastCopy()) continue;                                                                                                                  
   //if ( !mcTools.IsQuark(genP_pdgId) ) continue;                                                                                                   
   // if (!p.isPrompt()) continue;                                                                                                                   
   // if (!p.isPromptDecayed()) continue;                                                                                                            
   // if (!p.isPromptFinalState()) continue;                                                                                                         
   // if (!p.isDecayedLeptonHadron()) continue;                                                                                                      
   // if (!p.isTauDecayProduct()) continue;                                                                                                          
   // if (!p.isPromptTauDecayProduct()) continue;                                                                                                    
   // if (!p.isDirectTauDecayProduct()) continue;                                                                                                    
   // if (!p.isDirectPromptTauDecayProduct()) continue;                                                                                              
   // if (!p.isDirectPromptTauDecayProductFinalState()) continue;                                                                                    
   // if (!p.isDirectHadronDecayProduct()) continue;                                                                                                 
   // if (!p.isDirectHardProcessTauDecayProductFinalState()) continue;                                                                               
   //// if (!p.isHardProcess()) continue;                                                                                                            
   // if (!p.fromHardProcess()) continue;                       
   // if (!p.fromHardProcessDecayed()) continue;                                                                                                     
   // if (!p.fromHardProcessFinalState()) continue;                                                                                                  
   // if (!p.isHardProcessTauDecayProduct()) continue;                                                                                               
   // if (!p.isDirectHardProcessTauDecayProduct()) continue;                                                                                         
   // if (!p.fromHardProcessBeforeFSR()) continue;                                                                                                   
   // if (!p.isFirstCopy()) continue;                                                                                                                
   // if (!p.isLastCopyBeforeFSR()) continue;                                                                                                        

   // Add table rows                                                                                                                                 
   if(0) //to print the table with mothers/daugthers etc  if (1) //if(cfg_Verbose)
     {
       table.AddRowColumn(row, auxTools.ToString(entry)           );
       table.AddRowColumn(row, auxTools.ToString(genP_index)      );
       table.AddRowColumn(row, auxTools.ToString(genP_pdgId)      );
       table.AddRowColumn(row, auxTools.ToString(genP_status)     );
       table.AddRowColumn(row, auxTools.ToString(genP_charge)     );
       table.AddRowColumn(row, auxTools.ToString(genP_pt , 3)     );
       table.AddRowColumn(row, auxTools.ToString(genP_eta, 4)     );
       table.AddRowColumn(row, auxTools.ToString(genP_phi, 3)     );
       table.AddRowColumn(row, auxTools.ToString(genP_energy, 3)  );
       // table.AddRowColumn(row, "("+auxTools.ToString(genP_vtx.x(), 3)+", "+auxTools.ToString(genP_vtx.y(), 3) +", "+auxTools.ToString(genP_vtx.z(\), 3) + ")" );                                                                                                       
       table.AddRowColumn(row, "N/A" );
       table.AddRowColumn(row, "N/A" ); // auxTools.ToString(genP_Lxy, 3)                                                                            
       table.AddRowColumn(row, "N/A" ); // auxTools.ToString(genP_d0,  3)                                                                            
       if (genP_mothers.size() < 6)
	 {
	   table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genP_mothers) );
	 }
       else table.AddRowColumn(row, ".. Too many .." );
       if (genP_daughters.size() < 6)
	 {
	   table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genP_daughters) );
	 }
       else table.AddRowColumn(row, ".. Too many .." );
       row++;
     }

   // b-quarks                                                                                                                                       
   if(std::abs(genP_pdgId) == 5)
     {
       bQuarks_p4.push_back( genP_p4 );
       if ( mcTools.HasMother(p, +6) ) tbWPlus_BQuark_p4  = genP_p4;
       if ( mcTools.HasMother(p, -6) ) tbWMinus_BQuark_p4 = genP_p4;
     }// b-quarks                                                                                                                                    


   // W-                                                                                                                                             
   if (mcTools.HasMother(p, -24) )
     {
       if (genP_pdgId > 0) tbWPlus_Wqq_Quark_p4 = genP_p4;
       else tbWPlus_Wqq_AntiQuark_p4 = genP_p4;
     }

   // W+                                                                                                                                             
   if (mcTools.HasMother(p, +24) )
     {
       //std::cout << auxTools.ConvertIntVectorToString(genP_mothers) << std::endl;                                                                  
       if (genP_pdgId > 0) tbWMinus_Wqq_Quark_p4 = genP_p4;
       else tbWMinus_Wqq_AntiQuark_p4 = genP_p4;
     }

 }//for-loop: genParticles                                                                                                                           

 if (0) table.Print();// to print the table with mothers/daugthers etc if(1) //if (cfg_Verbose)


 //////////////////////////////////////////////////////////////////////////////////////////////////////   
 // GenJets                                                                                                    
 //////////////////////////////////////////////////////////////////////////////////////////////////////         

 //////////////////////////////////////////////////////////////////////////////////////////////////////        
 // MaxDiJet System (DiJet combination with largest mass)                                                      
 ////////////////////////////////////////////////////////////////////////////////////////////////////// 

 //////////////////////////////////////////////////////////////////////////////////////////////////////         
 // GenJets Correlations                                                                                     
 //////////////////////////////////////////////////////////////////////////////////////////////////////           

 ////////////////////////////////////////////////////////////////////////////////////////////////////// 
 // GenJets: Trijet with largest pT                                                        
 //////////////////////////////////////////////////////////////////////////////////////////////////////

 //////////////////////////////////////////////////////////////////////////////////////////////////////       
 // B-quarks (pT sorted)                                                                                  
 //////////////////////////////////////////////////////////////////////////////////////////////////////   

 //////////////////////////////////////////////////////////////////////////////////////////////////////            
 // BquarkPair - Ldg Jets Correlations                                                                         
 //////////////////////////////////////////////////////////////////////////////////////////////////////   

 //////kchristo/////////////////////////////////////////Old Version////////////////////////////////////////////////////////////////////////// 
 ////////////////////////////////////////// TT sample, bs & cs not from the tops  ///////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 if(0)
   {
     bool wehavesignal   = false;
     unsigned int wehave_extra_b = 0;
     unsigned int wehave_intocuts_extra_b = 0;
     unsigned int wehave_extra_c = 0;
     unsigned int wehave_intocuts_extra_c = 0;
     
     for (auto& p: fEvent.genparticles().getGenParticles())
       {
	 if(std::abs(p.pdgId()) == 5 && p.isFirstCopy()) // find the b
	   {
	     std::vector<short> b_mothers = p.mothers();
	     for (unsigned int i = 0; i < b_mothers.size() ; i++)
	       {
		 genParticle m = fEvent.genparticles().getGenParticles()[b_mothers.at(i)]; //create the particle object
		 if(std::abs(m.pdgId()) == 6)//if the b comes from top
		   {
		     wehavesignal   = true;
		     h_ttsample_bfromtop_pt     -> Fill (p.pt());
		     h_ttsample_bfromtop_eta    -> Fill (p.eta());
		   }
		 else 
		   {
		     wehave_extra_b++;
		     h_ttsample_bNOTfromtop_pt  -> Fill (p.pt());
		     h_ttsample_bNOTfromtop_eta -> Fill(p.eta());

		     h_tt_extra_b_overptCut    -> Fill("under 30 GeV",0);
		     h_tt_extra_b_underetaCut  -> Fill("|#eta| under 2.4 ",0);
		     h_tt_extra_b_intobothCuts -> Fill("into Cuts",0);

		     if(p.pt() < 30)             h_tt_extra_b_overptCut   -> Fill("under 30 GeV",1);
		     else                        h_tt_extra_b_overptCut   -> Fill("over 30 GeV",1);
		     if(std::abs(p.eta()) < 2.4) h_tt_extra_b_underetaCut -> Fill("|#eta| under 2.4 ",1);
		     else                        h_tt_extra_b_underetaCut -> Fill("|#eta| over 2.4 ",1);
		     if(std::abs(p.eta()) < 2.4 && p.pt() > 30) 
		       {
			 h_tt_extra_b_intobothCuts -> Fill("into Cuts",1);
			 wehave_intocuts_extra_b++;
		       }
		     else h_tt_extra_b_intobothCuts -> Fill("out of Cuts",1);

		   }
	       }//for mothers of bs
	   }    //if it is b

	 if(std::abs(p.pdgId()) == 4 && p.isFirstCopy()) // find the c
	   {
	     std::vector<short> c_mothers = p.mothers();
	     //std::cout << "------------------" << std::endl;
             for (unsigned int i = 0; i < c_mothers.size() ; i++)
	       {
		 unsigned int c_firstcopy_mother = GetTheFirstCopy(c_mothers.at(i));
		 genParticle m = fEvent.genparticles().getGenParticles()[c_firstcopy_mother]; //create the particle object
		 std::vector<short> c_grand_mothers = m.mothers();
		 for (unsigned int j = 0; j < c_grand_mothers.size() ; j++)
		   {
		     genParticle grandm = fEvent.genparticles().getGenParticles()[c_grand_mothers.at(j)]; //create the particle object
		     if (std::abs(grandm.pdgId()) == 6 || std::abs(m.pdgId()) == 6)        // if c comes from top 
		       {
			 h_ttsample_cfromtop_pt     -> Fill (p.pt());
			 h_ttsample_cfromtop_eta    -> Fill (p.eta()); 
		       }
		     else
		       {
			 wehave_extra_c++;
			 h_ttsample_cNOTfromtop_pt  -> Fill (p.pt());
			 //std::cout <<"mother:" << m.pdgId() << ", grandm:" << grandm.pdgId() << std::endl;
			 h_ttsample_cNOTfromtop_eta -> Fill (p.eta());

			 h_tt_extra_c_overptCut    -> Fill("under 30 GeV",0);
			 h_tt_extra_c_underetaCut  -> Fill("|#eta| under 2.4 ",0);
			 h_tt_extra_c_intobothCuts -> Fill("into Cuts",0);

			 if(p.pt() < 30)             h_tt_extra_c_overptCut   -> Fill("under 30 GeV",1);
			 else                        h_tt_extra_c_overptCut   -> Fill("over 30 GeV",1);
			 if(std::abs(p.eta()) < 2.4) h_tt_extra_c_underetaCut -> Fill("|#eta| under 2.4 ",1);
			 else                        h_tt_extra_c_underetaCut -> Fill("|#eta| over 2.4 ",1);
			 if(std::abs(p.eta()) < 2.4 && p.pt() > 30) 
			   {
			     h_tt_extra_c_intobothCuts -> Fill("into Cuts",1);
			     wehave_intocuts_extra_c++;
			   }
			 else h_tt_extra_c_intobothCuts -> Fill("out of Cuts",1);

		       }
		   } //for charm grand_mothers
	       }     //for charm mothers
	     //std::cout << "------------------" << std::endl;
	   }         //if it is c

       }             //for loop particle
     
     if (wehavesignal) 
       {
	 h_ttsample_extra_b_or_c       ->Fill("0 extra",0);
	 h_ttsample_extra_b_or_c       ->Fill("extra b",0);
	 h_ttsample_extra_b_or_c       ->Fill("2=<extra b",0);
	 h_ttsample_extra_b_or_c       ->Fill("extra c",0);
	 h_ttsample_extra_b_or_c       ->Fill("2=<extra c",0); //just to determine the name of the first bin, zero stands for 0 entries

	 if     (wehave_extra_b == 0 && wehave_extra_c ==0) h_ttsample_extra_b_or_c ->Fill("0 extra",1); 
	 if     (wehave_extra_b == 1) h_ttsample_extra_b_or_c                       ->Fill("extra b",1);
	 else if(wehave_extra_b >= 2) h_ttsample_extra_b_or_c                       ->Fill("2=<extra b",1);
	 if     (wehave_extra_c == 1) h_ttsample_extra_b_or_c                       ->Fill("extra c",1);
         else if(wehave_extra_c >= 2) h_ttsample_extra_b_or_c                       ->Fill("2=<extra c",1);
	 
	 h_ttsample_intocuts_extra_b_or_c       ->Fill("0 extra",0);
         h_ttsample_intocuts_extra_b_or_c       ->Fill("extra b",0);
         h_ttsample_intocuts_extra_b_or_c       ->Fill("2=<extra b",0);
         h_ttsample_intocuts_extra_b_or_c       ->Fill("extra c",0);
         h_ttsample_intocuts_extra_b_or_c       ->Fill("2=<extra c",0); //just to determine the name of the first bin, zero stands for 0 entries 
	 if     (wehave_intocuts_extra_b == 0 && wehave_intocuts_extra_c ==0) h_ttsample_intocuts_extra_b_or_c ->Fill("0 extra",1);
         if     (wehave_intocuts_extra_b == 1) h_ttsample_intocuts_extra_b_or_c                                ->Fill("extra b",1);
         else if(wehave_intocuts_extra_b >= 2) h_ttsample_intocuts_extra_b_or_c                                ->Fill("2=<extra b",1);
         if     (wehave_intocuts_extra_c == 1) h_ttsample_intocuts_extra_b_or_c                                ->Fill("extra c",1);
         else if(wehave_intocuts_extra_c >= 2) h_ttsample_intocuts_extra_b_or_c                                ->Fill("2=<extra c",1);
       }
     
   }


 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 //////kchristo/////////////////////////////////////////New Version////////////////////////////////////////////////////////////////////////// 
 /////////////////////////////////////////////// TT sample, bs not from the tops  /////////////////////////////////////////////////////////// 
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 if(1)
   {
     math::XYZTLorentzVector ldg_pT_extrab_p4(0,0,0,0), subldg_pT_extrab_p4(0,0,0,0) , ldg_pT_top_p4(0,0,0,0), subldg_pT_top_p4(0,0,0,0) ;
     bool wehavesignal   = false, wehaveldg_b = false, wehavesubldg_b = false, wehaveldg_top = false ;
     unsigned int wehave_intocuts_extra_b = 0;
     unsigned int Numoftop = 0, Numofbfrom_top =0;

     std::vector<math::XYZTLorentzVector> top_p4(2), bfromtop_p4(2);
     math::XYZTLorentzVector p4_initializer(0,0,0,0);
     for (unsigned int i = 0; i < top_p4.size() ; i++) //initialize                                 
       {
	 top_p4.at(i) = p4_initializer;
	 bfromtop_p4.at(i) = p4_initializer;
       }

     for (auto& p: fEvent.genparticles().getGenParticles())
       {
         if(std::abs(p.pdgId()) == 5 && p.isFirstCopy()) // find the b                                                                        
           {
	     std::vector<short> b_mothers = p.mothers();
             for (unsigned int i = 0; i < b_mothers.size() ; i++)
               {
                 genParticle m = fEvent.genparticles().getGenParticles()[b_mothers.at(i)]; //create the particle object          
                 if(std::abs(m.pdgId()) == 6)//if the b comes from top                                                 
                   {
		     if(Numoftop>=2) 
		       {
			 Numoftop++;
			 continue; //if we have more than 2 tops
		       }
                     wehavesignal   = true;
		     top_p4.at(Numoftop)            = m.p4();
		     bfromtop_p4.at(Numofbfrom_top) = p.p4();
		     Numoftop++;
		     Numofbfrom_top++;
		   }
                 else
		   {
                     h_tt_extra_b_overptCut    -> Fill("under 30 GeV",0);      //Just to determinate the label of the first bin
                     h_tt_extra_b_underetaCut  -> Fill("|#eta| under 2.4 ",0); //same as above
                     h_tt_extra_b_intobothCuts -> Fill("into Cuts",0);         //same as above

                     if(p.pt() < 30)             h_tt_extra_b_overptCut   -> Fill("under 30 GeV",1);
                     else                        h_tt_extra_b_overptCut   -> Fill("over 30 GeV",1);
                     if(std::abs(p.eta()) < 2.4) h_tt_extra_b_underetaCut -> Fill("|#eta| under 2.4 ",1);
                     else                        h_tt_extra_b_underetaCut -> Fill("|#eta| over 2.4 ",1);

                     if(std::abs(p.eta()) < 2.4 && p.pt() > 30)
                       {
			 if(!wehaveldg_b)
			   {
			     ldg_pT_extrab_p4 = p.p4();
			     wehaveldg_b = true; 
			   }
			 else if(!wehavesubldg_b)          // if we have ldg in pt but not subldg
			   {
			     if(p.pt() > ldg_pT_extrab_p4.pt()) // check if this one is the ldg
			       {
				 subldg_pT_extrab_p4 = ldg_pT_extrab_p4;
				 ldg_pT_extrab_p4    = p.p4();
			       }
			     else subldg_pT_extrab_p4 = p.p4();
			     wehavesubldg_b= true;
			   }
			 else
			   {
			     if(p.pt() > ldg_pT_extrab_p4.pt()) // check if this one is the ldg          
			       {
                                 subldg_pT_extrab_p4 = ldg_pT_extrab_p4;
                                 ldg_pT_extrab_p4    = p.p4();
                               }
			     else if (p.pt() > subldg_pT_extrab_p4.pt()) // check if this one is the subldg
			       {
				 subldg_pT_extrab_p4 =p.p4();
			       }
			   }

                         h_tt_extra_b_intobothCuts -> Fill("into Cuts",1);
                         wehave_intocuts_extra_b++;
                       }
                     else h_tt_extra_b_intobothCuts -> Fill("out of Cuts",1);

                   }
               }//for mothers of bs                                                                          
           }    //if it is b   

	 if (wehavesignal)
	   {
	     h_ttsample_intocuts_extra_b  ->Fill("0",0);
	     h_ttsample_intocuts_extra_b  ->Fill("1",0);
	     h_ttsample_intocuts_extra_b  ->Fill("2",0);
	     h_ttsample_intocuts_extra_b  ->Fill("3",0);
	     h_ttsample_intocuts_extra_b  -> Fill("4",0);
	     h_ttsample_intocuts_extra_b  ->Fill("5 =<",0); //just to determine the labels of the bins, 0 stands for zero entries 

	     if      (wehave_intocuts_extra_b == 0) h_ttsample_intocuts_extra_b  ->Fill("0",1);
	     else if (wehave_intocuts_extra_b == 1) h_ttsample_intocuts_extra_b  ->Fill("1",1);
	     else if (wehave_intocuts_extra_b == 2) h_ttsample_intocuts_extra_b  ->Fill("2",1);
	     else if (wehave_intocuts_extra_b == 3) h_ttsample_intocuts_extra_b  ->Fill("3",1);
	     else if (wehave_intocuts_extra_b == 4) h_ttsample_intocuts_extra_b  ->Fill("4",1);
	     else if (wehave_intocuts_extra_b >= 5) h_ttsample_intocuts_extra_b  ->Fill("5 =< ",1);
	   }
       }// for gen particles
     
     if(wehaveldg_b && Numoftop < 3)
       {
	 for(unsigned int i = 0; i < top_p4.size() ; i++)
           {
             if(!wehaveldg_top)
               {
                 ldg_pT_top_p4 = top_p4.at(i);
                 wehaveldg_top = true;
               }
             else
               {
                 if(top_p4.at(i).pt() > ldg_pT_top_p4.pt()) // check if this one is the ldg    
                   {
                     subldg_pT_top_p4 = ldg_pT_top_p4;
                     ldg_pT_top_p4    = top_p4.at(i);
                   }
                 else subldg_pT_top_p4 = top_p4.at(i);
               }
           }
       }
     
     if(wehaveldg_b && wehavesubldg_b && Numoftop < 3)
       {	 
	 double deltaR_twoextrab = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,subldg_pT_extrab_p4);
	 h_ttsample_extrab_dR       -> Fill(deltaR_twoextrab);
	 math::XYZTLorentzVector extra_di_b_sump4 = ldg_pT_extrab_p4 + subldg_pT_extrab_p4;
	 h_ttsample_Massof_extradib -> Fill(extra_di_b_sump4.M());
	 
	 for(unsigned int i = 0; i < top_p4.size() ; i++) 
	   {
	     double deltaR_top_ldgextrab    = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,top_p4.at(i)); 
	     double deltaR_top_subldgextrab = ROOT::Math::VectorUtil::DeltaR(subldg_pT_extrab_p4,top_p4.at(i));
	     h_ttsample_top_extrab_dR -> Fill(deltaR_top_ldgextrab);
	     h_ttsample_top_extrab_dR -> Fill(deltaR_top_subldgextrab);
	     
	     // Changed now
	     //if(deltaR_top_ldgextrab < deltaR_top_subldgextrab) h_ttsample_top_extrab_mindR -> Fill(deltaR_top_ldgextrab);
	     //else                                               h_ttsample_top_extrab_mindR -> Fill(deltaR_top_subldgextrab);
	     h_ttsample_top_extrab_mindR -> Fill(deltaR_top_ldgextrab);
	     h_ttsample_top_extrab_mindR -> Fill(deltaR_top_subldgextrab); 
	     //

	     if(bfromtop_p4.at(i) == p4_initializer) continue;
	     double deltaR_bfromtop_ldgextrab    = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,bfromtop_p4.at(i));
             double deltaR_bfromtop_subldgextrab = ROOT::Math::VectorUtil::DeltaR(subldg_pT_extrab_p4,bfromtop_p4.at(i));
             h_ttsample_bfromtop_extrab_dR -> Fill(deltaR_bfromtop_ldgextrab);
             h_ttsample_bfromtop_extrab_dR -> Fill(deltaR_bfromtop_subldgextrab);
	     
	     // Changed now                                                                                                                                                                             
             //if(deltaR_bfromtop_ldgextrab < deltaR_bfromtop_subldgextrab) h_ttsample_bfromtop_extrab_mindR -> Fill(deltaR_bfromtop_ldgextrab);
             //else                                                      h_ttsample_bfromtop_extrab_mindR -> Fill(deltaR_bfromtop_subldgextrab);
	     h_ttsample_bfromtop_extrab_mindR -> Fill(deltaR_bfromtop_ldgextrab);
	     h_ttsample_bfromtop_extrab_mindR -> Fill(deltaR_bfromtop_subldgextrab);
	     //
	   }

	 double deltaR_ldgtop_ldgextrab         = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,ldg_pT_top_p4);
	 double deltaR_subldgtop_ldgextrab      = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,subldg_pT_top_p4);
	 double deltaR_ldgtop_subldgextrab      = ROOT::Math::VectorUtil::DeltaR(subldg_pT_extrab_p4,ldg_pT_top_p4);
	 double deltaR_subldgtop_subldgextrab   = ROOT::Math::VectorUtil::DeltaR(subldg_pT_extrab_p4,subldg_pT_top_p4);
	 
	 h_ttsample_ldgtop_extrab_dR_Vs_subldgtop_extrab_dR   -> Fill(deltaR_ldgtop_ldgextrab,deltaR_subldgtop_ldgextrab);
	 h_ttsample_ldgtop_extrab_dR_Vs_subldgtop_extrab_dR   -> Fill(deltaR_ldgtop_subldgextrab,deltaR_subldgtop_subldgextrab);
       }

     else if(wehaveldg_b && !wehavesubldg_b && Numoftop < 3)
       {
	 for(unsigned int i = 0; i < top_p4.size() ; i++)
           {
             double deltaR_top_ldgextrab    = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,top_p4.at(i));
             h_ttsample_top_extrab_dR    -> Fill(deltaR_top_ldgextrab);
             h_ttsample_top_extrab_mindR -> Fill(deltaR_top_ldgextrab);

             if(bfromtop_p4.at(i) == p4_initializer) continue;
	     double deltaR_bfromtop_ldgextrab    = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,bfromtop_p4.at(i));
             h_ttsample_bfromtop_extrab_dR    -> Fill(deltaR_bfromtop_ldgextrab);
             h_ttsample_bfromtop_extrab_mindR -> Fill(deltaR_bfromtop_ldgextrab);
           }
	 
	 double deltaR_ldgtop_ldgextrab         = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,ldg_pT_top_p4);
         double deltaR_subldgtop_ldgextrab      = ROOT::Math::VectorUtil::DeltaR(ldg_pT_extrab_p4,subldg_pT_top_p4);
	 h_ttsample_ldgtop_extrab_dR_Vs_subldgtop_extrab_dR   -> Fill(deltaR_ldgtop_ldgextrab,deltaR_subldgtop_ldgextrab);

       }
     
   }

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 //////kchristo/////////////////////////////////////////New Version////////////////////////////////////////////////////////////////////////// 
 //////////////////////////Applied on HiggsMC the TT sample, bs not from the tops as the b from HIggs////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

 if(0)
   {
     math::XYZTLorentzVector bfromH_as_extrab_p4(0,0,0,0), ldg_pT_top_p4(0,0,0,0), subldg_pT_top_p4(0,0,0,0) ;
     bool wehavesignal   = false, wehave_bfromH = false, wehaveldg_top = false ;
     unsigned int wehave_intocuts_bfromH_as_extra_b = 0;
     unsigned int Numoftop = 0, Numofbfrom_top =0;
     
     std::vector<math::XYZTLorentzVector> top_p4(2), bfromtop_p4(2);
     math::XYZTLorentzVector p4_initializer(0,0,0,0);
     for (unsigned int i = 0; i < top_p4.size() ; i++) //initialize                                                 
       {
         top_p4.at(i) = p4_initializer;
         bfromtop_p4.at(i) = p4_initializer;
       }
     
     for (auto& p: fEvent.genparticles().getGenParticles())
       {
         if(std::abs(p.pdgId()) == 5 && p.isFirstCopy()) // find the b                                                                        
           {
	     std::vector<short> b_mothers = p.mothers();
             for (unsigned int i = 0; i < b_mothers.size() ; i++)
               {
                 genParticle m = fEvent.genparticles().getGenParticles()[b_mothers.at(i)]; //create the particle object                    
                 if(std::abs(m.pdgId()) == 6)//if the b comes from top                                                          
                   {
                     if(Numoftop>=2)
                       {
                         Numoftop++;
                         continue; //if we have more than 2 tops                                                            
                       }
                     wehavesignal   = true;
                     top_p4.at(Numoftop)            = m.p4();
                     bfromtop_p4.at(Numofbfrom_top) = p.p4();
                     Numoftop++;
                     Numofbfrom_top++;
		   }
                 else if(std::abs(m.pdgId()) == 37)
                   {
		     h_bfromH_as_extra_b_overptCut    -> Fill("under 30 GeV",0);      //Just to determinate the label of the first bin       
		     h_bfromH_as_extra_b_underetaCut  -> Fill("|#eta| under 2.4 ",0); //same as above                                  
		     h_bfromH_as_extra_b_intobothCuts -> Fill("into Cuts",0);         //same as above                               
		     
		     if(p.pt() < 30)             h_bfromH_as_extra_b_overptCut   -> Fill("under 30 GeV",1);
                     else                        h_bfromH_as_extra_b_overptCut   -> Fill("over 30 GeV",1);
                     if(std::abs(p.eta()) < 2.4) h_bfromH_as_extra_b_underetaCut -> Fill("|#eta| under 2.4 ",1);
                     else                        h_bfromH_as_extra_b_underetaCut -> Fill("|#eta| over 2.4 ",1);
		     
                     if(std::abs(p.eta()) < 2.4 && p.pt() > 30)
                       {
                         if(!wehave_bfromH)
                           {
                             bfromH_as_extrab_p4 = p.p4();
                             wehave_bfromH = true;
                           }
                      			 
                         h_bfromH_as_extra_b_intobothCuts -> Fill("into Cuts",1);
                         wehave_intocuts_bfromH_as_extra_b++;
                       }
                     else h_bfromH_as_extra_b_intobothCuts -> Fill("out of Cuts",1);
		   }
               }//for mothers of bs                                                                              
           }    //if it is b                                                                          
	 
         if (wehavesignal)
	 {
	   h_bfromH_as_extra_b_intocuts  ->Fill("0",0);
           h_bfromH_as_extra_b_intocuts  ->Fill("1",0);
	   
           if      (wehave_intocuts_bfromH_as_extra_b == 0) h_bfromH_as_extra_b_intocuts  ->Fill("0",1);
           else if (wehave_intocuts_bfromH_as_extra_b == 1) h_bfromH_as_extra_b_intocuts  ->Fill("1",1);
	 }
       }// for gen particles                                                                                      
     
     if(wehave_bfromH && Numoftop < 3)
       {
         for(unsigned int i = 0; i < top_p4.size() ; i++)
           {
             if(!wehaveldg_top)
	       {
                 ldg_pT_top_p4 = top_p4.at(i);
                 wehaveldg_top = true;
               }
             else
               {
                 if(top_p4.at(i).pt() > ldg_pT_top_p4.pt()) // check if this one is the ldg                       
                   {
                     subldg_pT_top_p4 = ldg_pT_top_p4;
                     ldg_pT_top_p4    = top_p4.at(i);
                   }
                 else subldg_pT_top_p4 = top_p4.at(i);
               }
           }
       }
     
     if(wehave_bfromH && Numoftop < 3)
       {
	 for(unsigned int i = 0; i < top_p4.size() ; i++)
	   {
	     double deltaR_top_bfromH_as_extrab    = ROOT::Math::VectorUtil::DeltaR(bfromH_as_extrab_p4,top_p4.at(i));
	     h_bfromH_as_extrab_top_dR    -> Fill(deltaR_top_bfromH_as_extrab);
	     h_bfromH_as_extrab_top_mindR -> Fill(deltaR_top_bfromH_as_extrab);
	     
	     if(bfromtop_p4.at(i) == p4_initializer) continue;
	     double deltaR_bfromtop_bfromH_as_extrab    = ROOT::Math::VectorUtil::DeltaR(bfromH_as_extrab_p4,bfromtop_p4.at(i));
	     h_bfromH_as_extrab_bfromtop_dR    -> Fill(deltaR_bfromtop_bfromH_as_extrab);
	     h_bfromH_as_extrab_bfromtop_mindR -> Fill(deltaR_bfromtop_bfromH_as_extrab);
	   }

	 double deltaR_ldgtop_bfromH_as_extrab     = ROOT::Math::VectorUtil::DeltaR(bfromH_as_extrab_p4,ldg_pT_top_p4);
         double deltaR_subldgtop_bfromH_as_extrab  = ROOT::Math::VectorUtil::DeltaR(bfromH_as_extrab_p4,subldg_pT_top_p4);
         h_bfromH_as_extrab_ldgtop_dR_Vs_subldgtop_extrab_dR  -> Fill(deltaR_ldgtop_bfromH_as_extrab,deltaR_subldgtop_bfromH_as_extrab);
       }
     
   }

 return;
}


///////kchristo/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////FInd the Last Copy Of a Particle////////////////////////////////////////////////
// by index/////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int ExtrabQuarksInTT::GetTheLastCopy(unsigned int firstcopy_index)
{
  unsigned int lastcopy_index = firstcopy_index; //in case if the first and last copy is the same index

  //For-loop: GenParticles                      
  for (auto& p: fEvent.genparticles().getGenParticles())
    {
      if (p.index() == firstcopy_index) // find the particle with the same index
	{
	  std::vector<short> genP_daughters = p.daughters(); //take the daughters of the particle
	  for (unsigned int i = 0; i < genP_daughters.size() ; i++) //loop over the daughters
	    {
	      genParticle d; //create the particle object                                                                  
	      d =  fEvent.genparticles().getGenParticles()[genP_daughters.at(i)];
	      if (d.pdgId() == p.pdgId()) // if the daughter is the same particle with the mother then
		firstcopy_index = d.index(); //find the daughter in the particle loop and take her's daughter
	      else
		lastcopy_index = p.index(); // if its not, it means that we found last copy
	    }
	}
    }
  return lastcopy_index;
}
		

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////kchristo/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////FInd the First Copy Of a Particle////////////////////////////////////////////////  
// by index/////////////////////////////////////////////////////////////////////////////////////////////////////            

unsigned int ExtrabQuarksInTT::GetTheFirstCopy(unsigned int lastcopy_index)
{
  unsigned int firstcopy_index = lastcopy_index; //in case if the first and last copy is the same index      
  //For-loop: GenParticles                                                                      
  for (auto& p: fEvent.genparticles().getGenParticles())
    {
      if (p.index() == lastcopy_index) // find the particle with the same index      
        {
	  std::vector<short> genP_mothers = p.mothers(); //take the mothers of the particle          
          for (unsigned int i = 0; i < genP_mothers.size() ; i++) //loop over the mothers   
            {
              genParticle m; //create the particle object           
              m =  fEvent.genparticles().getGenParticles()[genP_mothers.at(i)];
              if (m.pdgId() == p.pdgId()) // if the daughter is the same particle with the mother then           
                firstcopy_index = m.index(); //find the daughter in the particle loop and take her's daughter                 
              else
                lastcopy_index = p.index(); // if its not, it means that we found last copy                     
            }
        }
    }
  return firstcopy_index;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<float> ExtrabQuarksInTT::GetMomentumTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
							 float &C,
							 float &D,
							 float &H2){

  // Tensor required for calculation of: Sphericity, Aplanarity, Planarity                                              
  // Need all particles in event to calculate kinematic variables. Use all tracks (ch. particles) instead.                 
  // Links:                                                                                                           
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_24/doc/html/d5/d29/EventShapeVariables_8h_source.html                  
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_24/doc/html/dd/d99/classEventShapeVariables.html            
  // Get the Linear Momentum Tensor                                                                              
                                     
  TMatrixDSym MomentumTensor = ComputeMomentumTensor(jets, 1.0);

  // Find the Momentum-Tensor EigenValues (Q1, Q2, Q3)                                                           
  TMatrixDSymEigen eigen(MomentumTensor);
  TVectorD eigenvals = eigen.GetEigenValues();

  // Store & Sort the eigenvalues                                                                                         
  vector<float> eigenvalues(3);
  eigenvalues.at(0) = eigenvals(0); // Q1                                                                                     
  eigenvalues.at(1) = eigenvals(1); // Q2                                                                               
  eigenvalues.at(2) = eigenvals(2); // Q3                                                                   
  sort( eigenvalues.begin(), eigenvalues.end(), std::greater<float>() );

  // Calculate the eigenvalues sum (Requirement: Q1 + Q2 + Q3 = 1)         
  float eigenSum = std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0);
  if ( (eigenSum - 1.0) > 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Q1+Q2+Q3=1. Found that Q1+Q2+Q3 = " << eigenSum << ", instead.";
    }

  // Save the final eigenvalues                                                                
  float Q1 = eigenvalues.at(0);
  float Q2 = eigenvalues.at(1);
  float Q3 = eigenvalues.at(2);

  // Sanity check on eigenvalues: Q1 >= Q2 >= Q3 (Q1 >= 0)                                                                 
  bool bQ1Zero = (Q1 >= 0.0);
  bool bQ1Q2   = (Q1 >= Q2);
  bool bQ2Q3   = (Q2 >= Q3);
  bool bInequality = bQ1Zero * bQ1Q2 * bQ2Q3;

  if ( !(bInequality) )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that eigenvalues are ordered as Q1 >= Q2 >= Q3 (Q1 >= 0). Q1 = " << Q1 << ", Q2 = " << Q2 << ", Q3 = " << Q3;
    }

  // Calculate the linear combinations C and D                                                   
  C  = 3*(Q1*Q2 + Q1*Q3 + Q2*Q3); // Used to measure the 3-jet structure. Vanishes for perfece 2-jet event. Related to the 2nd Fox-Wolfram Moment (H2)                                                                                                                          
  D  = 27*Q1*Q2*Q3; // Used to measure the 4-jet structure. Vanishes for a planar event                          
  H2 = 1-C; // The C-measure is related to the second Fox-Wolfram moment (see below), $C = 1 - H_2$.          

  // C                                                                                    
  if ( abs(C-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the quantity C satisfies the inequality: 0.0 <= C <= 1.0. Found that C = " << C;
    }

  // D                                                                     
  if ( abs(C-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the quantity C satisfies the inequality: 0.0 <= D <= 1.0. Found that D = " << D;
    }

// 2nd Fox-Wolfram Moment                                                                                                                           
if ( abs(H2-1.0) > 1 + 1e-4 )
  {
    throw hplus::Exception("LogicError") << "Failure of requirement that the 2nd Fox-Wolfram Moment (H2) satisfies the inequality: 0.0 <= H2 <= 1.0\
. Found that H2 = " << H2;
  }

 if (0)
   {

     Table vars("Variable | Value | Allowed Range | Definition", "Text"); //LaTeX or Text                                                            
     vars.AddRowColumn(0, "C");
     vars.AddRowColumn(0, auxTools.ToString(C) );
     vars.AddRowColumn(0, "0.0 <= C <= 1.0");
     vars.AddRowColumn(0, "C = 3 x (Q1Q2 + Q1Q3 + Q2Q3");
     //                                                                                                                                              
     vars.AddRowColumn(1, "D");
     vars.AddRowColumn(1, auxTools.ToString(D) );
     vars.AddRowColumn(1, "0.0 <= D <= 1.0");
     vars.AddRowColumn(1, "D = 27 x Q1 x Q2 x Q3");
     //                                                                                                                                              
     vars.AddRowColumn(2, "2nd F-W Moment");
     vars.AddRowColumn(2, auxTools.ToString(H2) );
     vars.AddRowColumn(2, "0.0 <= H2 <= 1.0 ");
     vars.AddRowColumn(2, "H2 = 1-C");
     vars.Print();
   }

 return eigenvalues;
}


vector<float> ExtrabQuarksInTT::GetMomentumTensorEigenValues2D(std::vector<math::XYZTLorentzVector> jets,
							   float &Circularity){

  // For Circularity, the momentum tensor is the 2¡Á2 submatrix of Mjk, normalized by the sum of pT instead by the sum of       
  // This matrix has two eigenvalues Qi with 0 < Q1 < Q2. The following definition for the circularity C has been used:      
  //  C = 2 ¡Á min (Q1,Q2) / (Q1 +Q2)                                                                                         
  // The circularity C is especially interesting for hadron colliders because it only uses the momentum values in x and y direction transverse
  // to the beam line. So C is a two dimensional event shape variable and is therefore independent from a boost along z.                  
  // In addition, the normalization by the sum of the particle momenta makes C highly independent from energy calibration effects  (systematic uncertainty).                                                                                                                                              
  // C takes small values for linear and high values for circular events.                                       

  // Find the Momentum-Tensor EigenValues (E1, E2)                                                 
  TMatrixDSym MomentumTensor = ComputeMomentumTensor2D(jets);
  TMatrixDSymEigen eigen(MomentumTensor);
  TVectorD eigenvals = eigen.GetEigenValues();

  // Store & Sort the eigenvalues                                                                                       
  vector<float> eigenvalues(2);
  eigenvalues.at(0) = eigenvals[0]; // Q1                                                                
  eigenvalues.at(1) = eigenvals[1]; // Q2                                                   
  sort( eigenvalues.begin(), eigenvalues.end(), std::greater<float>() );
  
  // Save the final eigenvalues      
  float Q1 = eigenvalues.at(0);
  float Q2 = eigenvalues.at(1);

  // Sanity check on eigenvalues: (Q1 > Q2)                                                                        
  if (Q1 == 0) return eigenvalues;

  bool bInequality = (Q1 > 0 && Q1 > Q2);
  if ( !(bInequality) )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that eigenvalues are ordered as Q1 >= Q2. Found Q1 = " << Q1 << ", Q2 " << Q2;
    }

  // Calculate circularity                                                                           
  Circularity = 2*std::min(Q1, Q2)/(Q1+Q2); // is this definition correct?                  

  if (0)
    {
      Table vars("Variable | Value | Allowed Range | Definition", "Text"); //LaTeX or Text        
      vars.AddRowColumn(0, "Circularity");
      vars.AddRowColumn(0, auxTools.ToString(Circularity) );
      vars.AddRowColumn(0, "0.0 <= C <= 1.0 ");
      vars.AddRowColumn(0, "C = 2 ¡Á min (Q1,Q2)/(Q1+Q2)");
      vars.Print();
    }

  return eigenvalues;
}



vector<float> ExtrabQuarksInTT::GetSphericityTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
							   float &y23, float &Sphericity, float &SphericityT, float &Aplanarity,
							   float &Planarity, float &Y){

  // C, D parameters        
  // Need all particles in event to calculate kinematic variables. Use all tracks (ch. particles) instead.                
  // Links:                                                                                                                  
  // http://home.fnal.gov/~mrenna/lutp0613man2/node234.html                                                                                          

  // Sanity check: at least 3 jets (for 3rd-jet resolution)                                                                       
  if( (jets.size()) < 3 )
    {
      vector<float> zeros(3, 0);
      return zeros;
    }

  // Sort the jets by pT (leading jet first)                              
  std::sort( jets.begin(), jets.end(), PtComparator() );

  // Get the Sphericity Tensor                                         
  TMatrixDSym SphericityTensor = ComputeMomentumTensor(jets, 2.0);

  // Find the Momentum-Tensor EigenValues (Q1, Q2, Q3)                                   
  TMatrixDSymEigen eigen(SphericityTensor);
  TVectorD eigenvals = eigen.GetEigenValues();

  // Store & Sort the eigenvalues                                
  vector<float> eigenvalues(3);
  eigenvalues.at(0) = eigenvals(0); // Q1                                      
  eigenvalues.at(1) = eigenvals(1); // Q2                                 
  eigenvalues.at(2) = eigenvals(2); // Q3                                                                                   
  sort( eigenvalues.begin(), eigenvalues.end(), std::greater<float>() );

  // Calculate the eigenvalues sum (Requirement: Q1 + Q2 + Q3 = 1)                                              
  float eigenSum = std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0);
  if ( (eigenSum - 1.0) > 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Q1+Q2+Q3=1. Found that Q1+Q2+Q3 = " << eigenSum << ", instead.";
    }

  // Save the final eigenvalues                                                                                                                       
  float Q1 = eigenvalues.at(0);
  float Q2 = eigenvalues.at(1);
  float Q3 = eigenvalues.at(2);

  // Sanity check on eigenvalues: Q1 >= Q2 >= Q3 (Q1 >= 0)                                           
  bool bQ1Zero = (Q1 >= 0.0);
  bool bQ1Q2   = (Q1 >= Q2);
  bool bQ2Q3   = (Q2 >= Q3);
  bool bInequality = bQ1Zero * bQ1Q2 * bQ2Q3;

  if ( !(bInequality) )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that eigenvalues are ordered as Q1 >= Q2 >= Q3 (Q1 >= 0)";
    }


  // Calculate the event-shape variables                                                            
  float pT3Squared  = pow(jets.at(2).Pt(), 2);
  float HT2Squared  = pow(jets.at(0).Pt() + jets.at(1).Pt(), 2);
  y23         = pT3Squared/HT2Squared;
  Sphericity  = -1.0;
  SphericityT = -1.0;
  Aplanarity  = -1.0;
  Planarity   = -1.0;
  Y = (sqrt(3.0)/2.0)*(Q2-Q3); // (Since Q1>Q2, then my Q1 corresponds to Q3 when Q's are reversed-ordered). Calculate the Y (for Y-S plane)

  // Check the value of the third-jet resolution                                                                                   
  if (abs(y23-0.25) > 0.25 + 1e-4)
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that y23 satisfies the inequality: 0.0 <= y23 <= 0.25. Found that y23 = " << y23;
    }

  // Calculate the Sphericity (0 <= S <= 1). S~0 for a 2-jet event, and S~1 for an isotropic one          
  Sphericity = 1.5*(Q2 + Q3);
  if ( abs(Sphericity-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Sphericity (S) satisfies the inequality: 0.0 <= S <= 1.0. Found that S = " << Sphericity;
    }

  // Calculate the Sphericity (0 <= S <= 1). S~0 for a 2-jet event, and S~1 for an isotropic one                          
  SphericityT = 2.0*Q2/(Q1 + Q2);
  if ( abs(SphericityT-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Transverse Sphericity (ST) satisfies the inequality: 0.0 <= ST <= 1.0. Found that ST = " << SphericityT;
    }

  // Calculate the Aplanarity (0 <= A <= 0.5).  It measures the transverse momentum component out of the event plane             
  // A~0 for a planar event, A~0.5 for an isotropic one                                                                   
  Aplanarity = 1.5*(Q3);
  if ( abs(Aplanarity-0.5) > 0.5 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Aplanarity (A) satisfies the inequality: 0.0 <= A <= 0.5";
    }

  // Calculate the Aplanarity (0 <= P <= 0.5)                                                                 
  Planarity  = (2.0/3.0)*(Sphericity-2*Aplanarity);
  if ( abs(Planarity-0.5) > 0.5 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Planarity (P) satisfies the inequality: 0.0 <= P <= 0.5";
    }

  if (0)
    {

      Table vars("Variable | Value | Allowed Range | Definition", "Text"); //LaTeX or Text                        
      vars.AddRowColumn(0, "y23");
      vars.AddRowColumn(0, auxTools.ToString(y23) );
      vars.AddRowColumn(0, "0.0 <= y23 <= 0.25");
      vars.AddRowColumn(0, "y23 = pow(jet3_Pt, 2) / pow(jet1_Pt + jet2_Pt, 2)" );
      //                              
      vars.AddRowColumn(1, "Sphericity");
      vars.AddRowColumn(1, auxTools.ToString(Sphericity) );
      vars.AddRowColumn(1, "0.0 <= S <= 1.0");
      vars.AddRowColumn(1, "S = 1.5 x (Q2 + Q3)");
      //                                                                                     
      vars.AddRowColumn(2, "Sphericity (T)");
      vars.AddRowColumn(2, auxTools.ToString(SphericityT) );
      vars.AddRowColumn(2, "0.0 <= S (T) <= 1.0");
      vars.AddRowColumn(2, "S (T) = (2 x Q2)/(Q1 + Q2)");
      //                                                                                                               
      vars.AddRowColumn(3, "Aplanarity");
      vars.AddRowColumn(3, auxTools.ToString(Aplanarity) );
      vars.AddRowColumn(3, "0.0 <= A <= 0.5 ");
      //                                 
      vars.AddRowColumn(4, "Planarity");
      vars.AddRowColumn(4, auxTools.ToString(Planarity) );
      vars.AddRowColumn(4, "0.0 <= P <= 0.5 ");
      vars.AddRowColumn(4, "P (2/3) x (S - 2A)");
      //                                                                                                
      vars.AddRowColumn(5, "Y");
      vars.AddRowColumn(5, auxTools.ToString(Y) );
      vars.AddRowColumn(5, "");
      vars.AddRowColumn(5, "Y = sqrt(3)/2 x (Q1 - Q2)");
      vars.Print();
    }

  return eigenvalues;

}

TMatrixDSym ExtrabQuarksInTT::ComputeMomentumTensor(std::vector<math::XYZTLorentzVector> jets, double r)
{

  // r = 2: Corresponds to sphericity tensor (Get: Sphericity, Aplanarity, Planarity, ..)                      
  // r = 1: Corresponds to linear measures (Get: C, D, Second Fox-Wolfram moment, ...)                              
  TMatrixDSym momentumTensor(3);
  momentumTensor.Zero();

  if (r!=1.0 && r!=2.0)
    {
      throw hplus::Exception("LogicError") << "Invalid value r-value in computing the Momentum Tensor (r=" << r << ").  Supported valued are r=2.0 and r=1.0.";
    }

  // Sanity Check                                                                          
  if ( jets.size() < 2 )
    {
      return momentumTensor;
    }

  // Declare the Matrix normalisation (sum of momentum magnitutes to power r). That is: sum(|p|^{r})         
  double normalisation = 0.0;
  double trace = 0.0;

  // For-loop: Jets                                                                        
  for (auto& jet: jets){

    // Get the |p|^2 of the jet                                                                   
    double p2 = pow(jet.P(), 2); // jet.P();                                                      

    // For r=2, use |p|^{2}, for r=1 use |p| as the momentum weight                 
    double pR = ( r == 2.0 ) ? p2 : TMath::Power(p2, 0.5*r);

    // For r=2, use |1|, for r=1 use   (|p|^{2})^{-0.5} = |p|^{2 (-1/2)} = |p|^{-1} = 1.0/|p| 
    double pRminus2 = ( r == 2.0 ) ? 1.0 : TMath::Power(p2, 0.5*r - 1.0); //    

    // Add pR to the matrix normalisation factor                       
    normalisation += pR;

    // Fill the momentum (r=1) or  sphericity (r=2) tensor (Must be symmetric: Mij = Mji)                                      
    momentumTensor(0,0) += pRminus2*jet.px()*jet.px(); // xx     
    momentumTensor(0,1) += pRminus2*jet.px()*jet.py(); // xy                   
    momentumTensor(0,2) += pRminus2*jet.px()*jet.pz(); // xz              

    momentumTensor(1,0) += pRminus2*jet.py()*jet.px(); // yx                                            
    momentumTensor(1,1) += pRminus2*jet.py()*jet.py(); // yy                                                
    momentumTensor(1,2) += pRminus2*jet.py()*jet.pz(); // yz                                         

    momentumTensor(2,0) += pRminus2*jet.pz()*jet.px(); // zx                                        
    momentumTensor(2,1) += pRminus2*jet.pz()*jet.py(); // zy                                                              
    momentumTensor(2,2) += pRminus2*jet.pz()*jet.pz(); // zz     

  }// for (auto& jet: jets){                                                                                      

  // Normalise the tensors to have unit trace (Mxx + Myy + Mzz = 1)                                           
  momentumTensor *= (1.0/normalisation);
  trace = momentumTensor(0,0) + momentumTensor(1,1) + momentumTensor(2,2);

  // Print the tensor                                                           
  if (0)
    {
      std::cout << "\nMomentum Tensor (r = " << r << "):" << std::endl;
      Table tensor(" |  | ", "Text"); //LaTeX or Text                                                        
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,0) ) );
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,1) ) );
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,2) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,0) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,1) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,2) ) );
      tensor.AddRowColumn(2, auxTools.ToString( momentumTensor(2,0) ) );
      tensor.AddRowColumn(2, auxTools.ToString( momentumTensor(2,1) ) );
      tensor.AddRowColumn(2, auxTools.ToString( momentumTensor(2,2) ) );
      tensor.AddRowColumn(3, "");
      tensor.AddRowColumn(4, "Normalisation");
      tensor.AddRowColumn(4, auxTools.ToString(normalisation));
      tensor.AddRowColumn(5, "IsSymmetric");
      tensor.AddRowColumn(5, auxTools.ToString(momentumTensor.IsSymmetric()));
      tensor.AddRowColumn(6, "Determinant");
      tensor.AddRowColumn(6, auxTools.ToString(momentumTensor.Determinant()));
      tensor.AddRowColumn(7, "Trace");
      tensor.AddRowColumn(7, auxTools.ToString(trace));
      tensor.Print(false);
    }

  if ( abs(trace-1.0) > 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the Momentum-Tensor (r = " << r << ") Trace is 1.0. Found that abs(trace-1) = " << abs(trace-1) << ", instead.";
    }

  return momentumTensor;
}



TMatrixDSym ExtrabQuarksInTT::ComputeMomentumTensor2D(std::vector<math::XYZTLorentzVector> jets)
{

  TMatrixDSym momentumTensor(3);
  momentumTensor.Zero();

  // Sanity Check                                     
  if ( jets.size() < 2 )
    {
      return momentumTensor;
    }

  // Declare the Matrix normalisation (sum of momentum magnitutes to power r). That is: sum(|p|^{r})              
  double normalisation = 0.0;
  double trace = 0.0;

  // For-loop: Jets                                                                               
  for (auto& jet: jets){

    // Get the pT                                                                                     
    double pT = jet.Pt();

    // Add pT to the matrix normalisation factor                                                                   
    normalisation += pow(pT,2);

    // Fill the two-dimensional momentum tensor                                      
    momentumTensor(0,0) += jet.px()*jet.px(); // xx                                  
    momentumTensor(0,1) += jet.px()*jet.py(); // xy                                                                
    momentumTensor(1,0) += jet.py()*jet.px(); // yx                                            
    momentumTensor(1,1) += jet.py()*jet.py(); // yy                                               

  }// for (auto& jet: jets){                                                                              

  // Normalise tensor to get the normalised 2-d momentum tensor                                                            
  momentumTensor *= (1.0/normalisation);
  trace = momentumTensor(0,0) + momentumTensor(1,1);

  // Print the tensor                                                                                  
  if (0)
    {
      std::cout << "\nNormalied 2-D Momentum Tensor"  << std::endl;
      Table tensor(" |  | ", "Text"); //LaTeX or Text                                     
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,0) ) );
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,1) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,0) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,1) ) );
      tensor.AddRowColumn(2, "");
      tensor.AddRowColumn(3, "Normalisation");
      tensor.AddRowColumn(3, auxTools.ToString(normalisation));
      tensor.AddRowColumn(4, "IsSymmetric");
      tensor.AddRowColumn(4, auxTools.ToString(momentumTensor.IsSymmetric()));
      tensor.AddRowColumn(5, "Determinant");
      tensor.AddRowColumn(5, auxTools.ToString(momentumTensor.Determinant()));
      tensor.AddRowColumn(6, "Trace");
      tensor.AddRowColumn(6, auxTools.ToString(trace));
      tensor.Print(false);
    }

  if ( abs(trace-1.0) > 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the 2D Momentum-Tensor Trace is 1.0. Found that abs(trace-1) = " << abs(trace-1) << ", instead.";
    }

  return momentumTensor;
}


double ExtrabQuarksInTT::GetAlphaT(std::vector<math::XYZTLorentzVector> jets,
			       float &HT,
			       float &JT,
			       float &MHT,
			       float &Centrality){

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////          
  /// AlphaT:                                                                                                   
  /// Calculates the AlphaT variable, defined as an N-jets. This definition reproduces the kinematics of a  
  // di-jet system by constructing two pseudo-jets, which balance one another in Ht.                           
  // The two pseudo-jets are formed from the combination of the N objects that minimizes the quantity               
  /// DeltaHt = |Ht_pseudoJet1 - Ht_pseudoJet2| of the pseudo-jets.                                                              
  //                                                                                                 
  // Detailed Explanation:                                                                  
  // The method "alphaT()" of this class takes as input all jets in the event and uses them to form      
  // two Pseudo-Jets to describe the event.                                                              
  // If there are NJets in a given event this means there are 2^{NJets-1} combinations to do this.       
  // The methods does exactly that and for the combination which minimises the quantity         
  // DeltaHt = Ht_PseudoJet1 - Ht_PseudoJet2,                                                                                                    
  // it calculates the quantity alphaT.                        
  // The method "alphaT()" employs a double loop to recreate all the possilbe jet combinations            
  // out of NJets, by the use of an NJets-binary system. For example, if NJets=5, the loop          
  // indices ("k" outside, "l" inside) run both from "k"=0 to "k"=2^{4}=16 . The upper limit of          
  // the outside loop is given by the expression:                                                                
  // 1<<(NJets-1)) = shift the number 1 by (NJets-1) positions to the left.                                        
  // So, for NJets=5  (i.e. 1  --> 1 0 0 0 0 )                                                                       
  // This is now the way we will represent grouping into 2 Pseudo-Jets. The 0's represent one group and the 1's the other.           
  // So, for example 1 0 0 0 0 means 1 jet forms Pseudo-Jet1 and 4 jets form Pseudo-Jet2.                     
  // Also, for example, 1 0 0 1 0 means 2 jets form Pseudo-Jet1 and 3 jets form Pseudo-Jet2.           
  // The inside loop performs a bitwise right shift of index "k" by "l" positions and then                
  // compares the resulting bit to 1. So, for "k"=0, all the resulting comparisons in the         
  // inside loop will result to 0, except the one with "l"=4.                                                
  // This gives the first combination: 0 0 0 0 0   ( i.e. 0 jets form Pseudo-Jet1 and 5 jets form Pseudo-Jet2 )          
  // For "k"=1 (00000001 in 8bit representation), the first comparison is 1, since k is shifted by zero positions 
  // and then compared to 1. The rest comparisons yield zero, since by shifting the bit by any position and comparing to 1 gives zero.          
  // Thus, for "k"=1 we have after the second loop: 0 0 0 0 1                                                                                        
  // In the same manner, we get for "k"=2 (00000001 in 8bit representation) we have after the second loop: 0 0 0 1 0     
  //  To summarise, for NJets=5 we have 16 combinations:                                                     
  // For "k"=0  ( 00000000 in 8bit representation) we have after the second loop: 0 0 0 0 0         
  // For "k"=1  ( 00000001 in 8bit representation) we have after the second loop: 0 0 0 0 1           
  // For "k"=2  ( 00000001 in 8bit representation) we have after the second loop: 0 0 0 1 0        
  // For "k"=3  ( 00000011 in 8bit representation) we have after the second loop: 0 0 0 1 1                 
  // For "k"=4  ( 00000100 in 8bit representation) we have after the second loop: 0 0 1 0 0                             
  // For "k"=5  ( 00000101 in 8bit representation) we have after the second loop: 0 0 1 0 1 
  // For "k"=6  ( 00000110 in 8bit representation) we have after the second loop: 0 0 1 1 0                       
  // For "k"=7  ( 00000111 in 8bit representation) we have after the second loop: 0 0 1 1 1        
  // For "k"=8  ( 00001000 in 8bit representation) we have after the second loop: 0 1 0 0 0     
  // For "k"=9  ( 00001001 in 8bit representation) we have after the second loop: 0 1 0 0 1        
  // For "k"=10 ( 00010000 in 8bit representation) we have after the second loop: 0 1 0 0 0               
  // For "k"=11 ( 00010001 in 8bit representation) we have after the second loop: 0 1 0 0 1               
  // For "k"=12 ( 00010010 in 8bit representation) we have after the second loop: 0 1 0 1 0           
  // For "k"=13 ( 00010011 in 8bit representation) we have after the second loop: 0 1 0 1 1                   
  // For "k"=14 ( 00010100 in 8bit representation) we have after the second loop: 0 1 1 0 0          
  // For "k"=15 ( 00010101 in 8bit representation) we have after the second loop: 0 1 1 0 1               
  // For "k"=16 ( 00010110 in 8bit representation) we have after the second loop: 0 1 1 1 0                                                          
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////      

  /// Declaration of variables                                                                                                 
  unsigned nJets = jets.size();

  // Sanity Check                                          
  if ( jets.size() < 2 )
    {
      return -1.0;
    }

  std::vector<float> vE, vEt, vPx, vPy, vPz;
  std::vector<bool> vPseudo_jet1;
  const bool bList = true;

  // For-Loop: All jets                                                                    
  for (auto& jet: jets){
    vE.push_back( jet.E() );
    vEt.push_back( jet.Et() );
    vPx.push_back( jet.Px() );
    vPy.push_back( jet.Py() );
    vPz.push_back( jet.Pz() );
  }

  // Calculate sums                                                                                   
  float fSum_e  = accumulate( vE.begin() , vE.end() , 0.0 );
  float fSum_et = accumulate( vEt.begin(), vEt.end(), 0.0 );
  float fSum_px = accumulate( vPx.begin(), vPx.end(), 0.0 );
  float fSum_py = accumulate( vPy.begin(), vPy.end(), 0.0 );

  // Minimum Delta Et for two pseudo-jets                                                   
  float fMin_delta_sum_et = -1.0;

  // Iterate through different combinations                                                           
  for ( unsigned k=0; k < unsigned(1<<(nJets-1)); k++ ) {
    float fDelta_sum_et = 0.0;
    std::vector<bool> jet;

    // Iterate through jets                                                          
    for ( unsigned l=0; l < vEt.size(); l++ ) {
      /// Bitwise shift of "k" by "l" positions to the right and compare to 1 (&1)                               
      /// i.e.: fDelta_sum_et += vEt[l] * ( 1 - 2*0 );  if comparison is un-successful                        
      ///  or   fDelta_sum_et += vEt[l] * ( 1 - 2*1 );  if comparison is successful                      
      // in this way you add up all Et from PseudoJetsGroupA (belonging to 0's group) and subtract that from PseudoJetsGroupB (1's group)     
      fDelta_sum_et += vEt[l] * ( 1 - 2 * (int(k>>l)&1) );
      if ( bList ) { jet.push_back( (int(k>>l)&1) == 0 ); }
    }

    // Find configuration with minimum value of DeltaHt                              
    if ( ( fabs(fDelta_sum_et) < fMin_delta_sum_et || fMin_delta_sum_et < 0.0 ) ) {
      fMin_delta_sum_et = fabs(fDelta_sum_et);
      if ( bList && jet.size() == vEt.size() ){vPseudo_jet1.resize(jet.size());}
    }

  }

  // Sanity check                                                                                       
  if ( fMin_delta_sum_et < 0.0 )
    {
      throw hplus::Exception("LogicError") << "Minimum Delta(Sum_Et) is less than zero! fMin_delta_sum_et = " << fMin_delta_sum_et;
    }

  // Calculate Event-Shape Variables                                                       
  float dHT = fMin_delta_sum_et;
  HT  = fSum_et;
  JT  = fSum_et - vEt.at(0); // Ht without considering the Ldg Jet of the Event                   
  MHT = sqrt(pow(fSum_px,2) + pow(fSum_py,2));
  Centrality = fSum_et/fSum_e;
  float AlphaT = ( 0.5 * ( HT - dHT ) / sqrt( pow(HT,2) - pow(MHT,2) ) );

  if (0)
    {

      Table vars("Variable | Value | Definition", "Text"); //LaTeX or Text          
      vars.AddRowColumn(0, "HT");
      vars.AddRowColumn(0, auxTools.ToString(HT) );
      vars.AddRowColumn(0, "HT = Sum(Jet_Et)");
      //                                                                                                                     
      vars.AddRowColumn(1, "JT");
      vars.AddRowColumn(1, auxTools.ToString(JT) );
      vars.AddRowColumn(1, "JT = Ht - Jet1_Et");
      //                                                             
      vars.AddRowColumn(2, "dHT");
      vars.AddRowColumn(2, auxTools.ToString(dHT) );
      vars.AddRowColumn(2, "DeltaHT = min[Delta(Pseudojet1_Et, Pseudojet2_Et)]");
      //                                                                              
      vars.AddRowColumn(3, "MHT");
      vars.AddRowColumn(3, auxTools.ToString(MHT) );
      vars.AddRowColumn(3, "MHT = sqrt( pow(Sum(px), 2) + pow(Sum(py), 2))");
      //                                                                                                
      vars.AddRowColumn(4, "AlphaT");
      vars.AddRowColumn(4, auxTools.ToString(AlphaT) );
      vars.AddRowColumn(4, "AlphaT = 0.5 x (HT - dHT) /sqr( pow(HT, 2) - pow(MHT, 2))");
      vars.Print();
    }

  return AlphaT;
}

vector<GenJet> ExtrabQuarksInTT::GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCuts, std::vector<float> etaCuts, vector<genParticle> genParticlesToMatch)
{
  /*                                                                                                                        
    Jet-Flavour Definitions (https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools)                             
                                                                                                 
    Algorithmic definition: (NOTE: Algorithmic definition is used by default for all b-tagging purposes)                     
    - Try to find the parton that most likely determines the properties of the jet and assign that flavour as the true flavour    
    - Here, the ¡°final state¡± partons (after showering, radiation) are analyzed (within ¦¤R < 0.3 of the reconstructed jet axis).       
    Partons selected for the algorithmic definition are those partons that don't have other partons as daughters,                                    
    without any explicit requirement on their status (for Pythia6 these are status=2 partons).                         
    - Jets from radiation are matched with full efficiency                                                  
    -If there is a b/c within the jet cone: label as b/c           
    -Otherwise: assign flavour of the hardest parton                                                           
  */

  // Definitions                                                                         
  std::vector<GenJet> jets;
  unsigned int jet_index   = -1;
  unsigned int ptCut_index  = 0;
  unsigned int etaCut_index = 0;

  // For-loop: All genParticles                                                                         
  for (auto& p: genParticlesToMatch)
    {

      // Comparison variables                                                                                 
      double dR   = 1e6;
      double dPt  = 1e6;
      double dEta = 1e6;
      double dPhi = 1e6;

      // For-loop: All Generated Jets      
      for (auto jet: genJets)
        {

          // Jet index (for pT and eta cuts)                                        
          jet_index++;

          dPt  = jet.pt() - p.pt();
          dEta = jet.eta() - p.eta();
          dPhi = jet.phi() - p.phi();
          dR   = ROOT::Math::VectorUtil::DeltaR(jet.p4(), p.p4());

          // Fail to match                                                                    
          if (dR > 0.3) continue;

          // Apply cuts after the matching                                  
          const float ptCut  = ptCuts.at(ptCut_index);
	  const float etaCut = etaCuts.at(etaCut_index);
	  if (jet.pt() < ptCut) continue;
          if (std::abs(jet.eta()) > etaCut) continue;

          // Save this particle                                           
          jets.push_back(jet);
          if (0) std::cout << "dR = " << dR << ": dPt = " << dPt << ", dEta = " << dEta << ", dPhi = " << dPhi << std::endl;

          // Increment cut index only. Cannot be bigger than the size of the cut list provided                    
	  if (ptCut_index  < ptCuts.size()-1  ) ptCut_index++;
	  if (etaCut_index < etaCuts.size()-1  ) etaCut_index++;
          break;
        }
    }
  if (0) std::cout << "bjets.size() = " << jets.size() << std::endl;
  return jets;
}

vector<GenJet> ExtrabQuarksInTT::GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCuts, std::vector<float> etaCuts)
{
  std::vector<GenJet> jets;

  // Definitions                                                  
  unsigned int jet_index   = -1;
  unsigned int ptCut_index  = 0;
  unsigned int etaCut_index = 0;

  // For-loop: All Generated Jets                                                     
  for (auto jet: genJets)
    {

      // Jet index (for pT and eta cuts)                                   
      jet_index++;

      // Apply cuts                     
      const float ptCut  = ptCuts.at(ptCut_index);
      const float etaCut = etaCuts.at(etaCut_index);
      if (jet.pt() < ptCut) continue;
      if (std::abs(jet.eta()) > etaCut) continue;

      // Save this particle                                                                                            
      jets.push_back(jet);

      // Increment cut index only. Cannot be bigger than the size of the cut list provided                                 
      if (ptCut_index  < ptCuts.size()-1  ) ptCut_index++;
      if (etaCut_index < etaCuts.size()-1  ) etaCut_index++;
    }
  return jets;
}

vector<genParticle> ExtrabQuarksInTT::GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy, const bool hasNoDaughters)
{

  std::vector<genParticle> particles;

  // For-loop: All genParticles                                 
  for (auto& p: genParticles)
    {
      // Find last copy of a given particle                                              
      if (isLastCopy) if (!p.isLastCopy()) continue;

      // Commonly enables for parton-based jet flavour definition                                             
      if (hasNoDaughters) if (p.daughters().size() > 0) continue;

      // Consider only particles                                                                     
      if (std::abs(p.pdgId()) != pdgId) continue;

      // Apply cuts                                                                                                        
      if ( p.pt() < ptCut ) continue;
      if (std::abs(p.eta()) > etaCut) continue;

      // Save this particle                                                                                                                          
      particles.push_back(p);
    }
  return particles;
}
