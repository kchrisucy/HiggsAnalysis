// The analyzer is based on the Kinematics analyzer.  Search for the "kchristo" to find out what I changed and my analysis on BoostTopologies
// The first attempt is under "Study Boosted Topologies" comment with sub-comment "old".  The whole thing is into an if(0) statement so it is
// disabled.  I 've cleaned the code from many unnecessary histograms and I put it under the "Study Boosted Topologies" comment with sub-comment "With fat Jets".
// Because the code is huge, due to of all the histos that i had to make, under the comment:
//                              " // another approach, first the matching with AK8 and after the dR between the objects"
// someone can find the final and most relevant analysis about Boosted Topologies, using AK8 fat-jets

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


class BoostedTopologies: public BaseSelector {
public:
  explicit BoostedTopologies(const ParameterSet& config, const TH1* skimCounters);
  virtual ~BoostedTopologies() {}

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
  //virtual vector<float> GetMomentumTensorEigenValues2D(std::vector<math::XYZTLorentzVector> jets, float &Circularity);
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

  WrappedTH1 *h_GenJets_N;

  
  //kchristo ---------------------------------------------------------------------------------------                                          

  // -----------Study Boosted Topologies----------
  WrappedTH1 *h_bfromH_Higgs_dR;
  WrappedTH2 *h_bfromH_Higgs_dR_Vs_Higgs_Pt;

  WrappedTH1 *h_objectsfromHiggstop_maxdR;
  WrappedTH2 *h_objectsfromHiggstop_maxdR_Vs_Higgstop_Pt;
  WrappedTH2 *h_objectsfromHiggstop_maxdR_Vs_Higgs_Pt;
  WrappedTH1 *h_objectsfromHiggstop_mindR;
  WrappedTH1 *h_objectsfromHiggstop_Prob_mindR_lt_p8;
  WrappedTH1 *h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than;
  WrappedTH2 *h_objectsfromHiggstop_mindR_Vs_Higgstop_Pt;
  WrappedTH2 *h_objectsfromHiggstop_mindR_Vs_Higgs_Pt;
  WrappedTH1 *h_bfromHiggstop_Higgstop_dR;
  WrappedTH2 *h_bfromHiggstop_Higgstop_dR_Vs_Higgstop_Pt;
  WrappedTH1 *h_b_higgs_underdR;
  WrappedTH1 *h_b_top_fromhiggs_underdR;

  WrappedTH1 *h_WfromHiggstop_closestb_dR;
  WrappedTH1 *h_WfromHiggstop_samedecayb_dR;
  WrappedTH1 *h_for_WbfromH_intofatjet_W_Pt;
  WrappedTH1 *h_for_WbfromH_intofatjet_Top_Pt;
  WrappedTH1 *h_Probdecayproductsintofatjet_Hs; // Probability of the decay products from Top from Higgs side to fall within 0.8 from the barycenter
  WrappedTH1 *h_Hs_QuarksintofatjetMultiplicity;
  WrappedTH1 *h_Hs_isbQuarkintofatjet;
  WrappedTH1 *h_Pairsinbarycenter_enoughdeltaR_Hs; // Probability of the decay products from Top from Higgs side that are in 0.8 from the barycenter to have between them more than 0.4 deltaR                                            
  WrappedTH2 *h_WfromHiggstop_closestb_dR_Vs_W_Pt;
  WrappedTH2 *h_WfromHiggstop_samedecayb_dR_Vs_W_Pt;
  WrappedTH2 *h_WfromHiggstop_closestb_dR_Vs_samedecaytop_Pt;
  WrappedTH2 *h_WfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt;
  
  // try with gen-jets ++++++++++++++++++++++
  WrappedTH1 *h_Hs_deltaR_bQuark_closestbJet;
  WrappedTH1 *h_Hs_deltaR_bQuark_fromH_closestbJet;
  WrappedTH1 *h_Hs_deltaR_QfromW_closestJet;
  WrappedTH1 *h_Hs_numof_obj_fromtop_matchedwith_genjet;
  WrappedTH1 *h_Hs_ProbdecayJetsintofatjet;
  WrappedTH1 *h_Hs_JetsintofatjetMultiplicity;
  WrappedTH1 *h_Hs_isbJetintofatjet;
  WrappedTH1Triplet *h_Hs_BjetInsideFatJet_Top_Pt;
  WrappedTH1Triplet *h_Hs_BjetInsideFatJet_W_Pt;
  WrappedTH1Triplet *h_Hs_BjetInsideFatJet_H_Pt;

  WrappedTH1 *h_NotHs_deltaR_bQuark_closestbJet;
  WrappedTH1 *h_NotHs_deltaR_QfromW_closestJet;
  WrappedTH1 *h_NotHs_numof_obj_fromtop_matchedwith_genjet;
  WrappedTH1 *h_NotHs_ProbdecayJetsintofatjet;
  WrappedTH1 *h_NotHs_JetsintofatjetMultiplicity;
  WrappedTH1 *h_NotHs_isbJetintofatjet;
  WrappedTH1Triplet *h_NotHs_BjetInsideFatJet_Top_Pt;
  WrappedTH1Triplet *h_NotHs_BjetInsideFatJet_W_Pt;
  // ++++++++++++++++++++++++++++++++++++++++
  
  // notfrom Higgs side
  WrappedTH1 *h_bNOTfromHiggstop_top_dR;
  WrappedTH1 *h_objectsNOTfromHiggstop_maxdR;
  WrappedTH2 *h_objectsNOTfromHiggstop_maxdR_Vs_top_Pt;
  WrappedTH1 *h_objectsNOTfromHiggstop_mindR;
  WrappedTH1 *h_objectsNOTfromHiggstop_Prob_mindR_lt_p8;
  WrappedTH1 *h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than;
  WrappedTH2 *h_objectsNOTfromHiggstop_mindR_Vs_top_Pt;
  WrappedTH2 *h_bNOTfromHiggstop_top_dR_Vs_top_Pt;
  WrappedTH1 *h_b_topNOTfromhiggs_underdR;
  
  WrappedTH1 *h_WNOTfromHiggstop_closestb_dR;
  WrappedTH1 *h_WNOTfromHiggstop_samedecayb_dR;
  WrappedTH1 *h_for_WbNOTfromH_intofatjet_W_Pt;
  WrappedTH1 *h_for_WbNOTfromH_intofatjet_Top_Pt;
  WrappedTH1 *h_Probdecayproductsintofatjet_NOTHs;// Probability of the decay products from Top NOT from Higgs side to fall within 0.8 from the barycenter
  WrappedTH1 *h_NotHs_QuarksintofatjetMultiplicity;
  WrappedTH1 *h_NotHs_isbQuarkintofatjet;         
  WrappedTH1 *h_Pairsinbarycenter_enoughdeltaR_NOTHs; // Probability of the decay products from Top NOT from Higgs side that are in 0.8 from the barycenter to have between them more than 0.4 deltaR
  WrappedTH2 *h_WNOTfromHiggstop_closestb_dR_Vs_W_Pt;
  WrappedTH2 *h_WNOTfromHiggstop_samedecayb_dR_Vs_W_Pt;
  WrappedTH2 *h_WNOTfromHiggstop_closestb_dR_Vs_samedecaytop_Pt;
  WrappedTH2 *h_WNOTfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt;
  

  WrappedTH1 *h_Hs_objectsfromtop_top_maxdR;
  WrappedTH1 *h_Hs_objectsfromtop_top_mindR;
  WrappedTH1 *h_Hs_which_objectfromtop_maxdR;
  WrappedTH1 *h_NotHs_objectsfromtop_top_maxdR;
  WrappedTH1 *h_NotHs_objectsfromtop_top_mindR;
  WrappedTH1 *h_NotHs_which_objectfromtop_maxdR;

  // All the plots with "Iffatjet" corresponds to the case where all 3 quarks are in 0.8 from barycenter(the top direction) 
  WrappedTH1 *h_Iffatjet_Hs_Top_pT;                //Hs stands for Higgs side
  WrappedTH1 *h_Iffatjet_Hs_W_pT;
  WrappedTH2 *h_Iffatjet_Hs_Top_pT_Vs_W_pT;
  WrappedTH1 *h_Iffatjet_Hs_bfromH_pT;
  WrappedTH1 *h_Iffatjet_Hs_EventsWithHighTop_pT;
  WrappedTH1 *h_Iffatjet_Hs_EventsWithHighW_pT;
  WrappedTH1 *h_Iffatjet_Hs_HighWandTop_pT_dRmin;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_Top_dEta;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_Top_dPhi;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_Top_dR;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_W_dEta;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_W_dPhi;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_W_dR;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_bfromH_dEta;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_bfromH_dPhi;
  WrappedTH1 *h_Iffatjet_Hs_bfromT_bfromH_dR;
  
  WrappedTH1 *h_Iffatjet_NotHs_Top_pT;                //NotHs stands for Not from Higgs side       
  WrappedTH1 *h_Iffatjet_NotHs_W_pT;
  WrappedTH2 *h_Iffatjet_NotHs_Top_pT_Vs_W_pT;
  WrappedTH1 *h_Iffatjet_NotHs_bfromH_pT;
  WrappedTH1 *h_Iffatjet_NotHs_EventsWithHighTop_pT;
  WrappedTH1 *h_Iffatjet_NotHs_EventsWithHighW_pT;
  WrappedTH1 *h_Iffatjet_NotHs_HighWandTop_pT_dRmin;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_Top_dEta;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_Top_dPhi;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_Top_dR;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_W_dEta;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_W_dPhi;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_W_dR;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_bfromH_dEta;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_bfromH_dPhi;
  WrappedTH1 *h_Iffatjet_NotHs_bfromT_bfromH_dR;
  
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__top_pT;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__W_pT;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__b_pT;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRmax_b_top;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb_dRmin_b_objfromW;
  WrappedTH2 *h_Hs_mostdistantfromtop_isb__dRqq_Vs_W_pT;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than;
  WrappedTH2 *h_Hs_mostdistantfromtop_isb__dRqq_Vs_top_pT;
  
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__top_pT;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__W_pT;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__b_pT;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRmax_b_top;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb_dRmin_b_objfromW;
  WrappedTH2 *h_NotHs_mostdistantfromtop_isb__dRqq_Vs_W_pT;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than;
  WrappedTH2 *h_NotHs_mostdistantfromtop_isb__dRqq_Vs_top_pT;
  
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7_tf;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf;

  //Study Boosted Topologies With fat Jets ----------------
  WrappedTH1 *h_Hs_QuarksFromW_deltaR;
  WrappedTH1 *h_Hs_QuarksFromW_Prob_deltaR;
  WrappedTH1 *h_Hs_QuarksintoBaryCenterMultiplicity;
  WrappedTH1 *h_Hs_isbQuarkintoBaryCenter;
  WrappedTH1 *h_Hs_OnlyQQ_dR_less_p7;
  WrappedTH1 *h_Hs_Prob_Diquark_match_with_fj;
  WrappedTH1 *h_Hs_MatchedWithDiquark_fj_pT;
  
  WrappedTH1 *h_Hs_objectsfromtop_dRmax_pTcuts;
  WrappedTH1 *h_Hs_objectsfromtop_dRmin_pTcuts;
  WrappedTH1 *h_Hs_QuarksFromW_deltaR_pTcuts;
  WrappedTH1 *h_Hs_objectsfromtop_Prob_dRmax_pTcuts;
  WrappedTH1 *h_Hs_objectsfromtop_Prob_dRmin_pTcuts;
  WrappedTH1 *h_Hs_QuarksFromW_Prob_deltaR_pTcuts;
  WrappedTH1 *h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts;
  WrappedTH1 *h_Hs_isbQuarkintoBaryCenter_pTcuts;
  WrappedTH1 *h_Hs_OnlyQQ_dR_less_p7_pTcuts;
  WrappedTH1 *h_Hs_Prob_Diquark_match_with_fj_pTcuts;
  WrappedTH1 *h_Hs_MatchedWithDiquark_fj_pT_pTcuts;
  WrappedTH1 *h_Hs_MatchedWithDiquark_Prob_fj_pT_pTcuts;
  WrappedTH1 *h_Hs_MatchedWithDiquark_Prob_fj_pT;
  WrappedTH1 *h_Hs_QuarksFromTop_Passed_pTcuts;
  
  WrappedTH1 *h_NotHs_QuarksFromW_deltaR;
  WrappedTH1 *h_NotHs_QuarksFromW_Prob_deltaR;
  WrappedTH1 *h_NotHs_QuarksintoBaryCenterMultiplicity;
  WrappedTH1 *h_NotHs_isbQuarkintoBaryCenter;
  WrappedTH1 *h_NotHs_OnlyQQ_dR_less_p7;
  WrappedTH1 *h_NotHs_Prob_Diquark_match_with_fj;
  WrappedTH1 *h_NotHs_MatchedWithDiquark_fj_pT;

  WrappedTH1 *h_NotHs_objectsfromtop_dRmax_pTcuts;
  WrappedTH1 *h_NotHs_objectsfromtop_dRmin_pTcuts;
  WrappedTH1 *h_NotHs_objectsfromtop_Prob_dRmax_pTcuts;
  WrappedTH1 *h_NotHs_objectsfromtop_Prob_dRmin_pTcuts;
  WrappedTH1 *h_NotHs_QuarksFromW_deltaR_pTcuts;
  WrappedTH1 *h_NotHs_QuarksFromW_Prob_deltaR_pTcuts;
  WrappedTH1 *h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts;
  WrappedTH1 *h_NotHs_isbQuarkintoBaryCenter_pTcuts;
  WrappedTH1 *h_NotHs_OnlyQQ_dR_less_p7_pTcuts;
  WrappedTH1 *h_NotHs_Prob_Diquark_match_with_fj_pTcuts;
  WrappedTH1 *h_NotHs_MatchedWithDiquark_fj_pT_pTcuts;
  WrappedTH1 *h_NotHs_MatchedWithDiquark_Prob_fj_pT_pTcuts;
  WrappedTH1 *h_NotHs_MatchedWithDiquark_Prob_fj_pT;
  WrappedTH1 *h_NotHs_QuarksFromTop_Passed_pTcuts;
  
  WrappedTH1 *h_Hs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj;
  WrappedTH1 *h_Hs_objectsfromtop_mindR_ltp8_matchedfj_pt;
  WrappedTH1 *h_Hs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt;
  WrappedTH1 *h_Hs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt;
  
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt;
  WrappedTH1 *h_Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt;


  WrappedTH1 *h_NotHs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj;
  WrappedTH1 *h_NotHs_objectsfromtop_mindR_ltp8_matchedfj_pt;
  WrappedTH1 *h_NotHs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt;
  WrappedTH1 *h_NotHs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt;

  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt;
  WrappedTH1 *h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt;
  
  // update 24.01.2018 after ptcuts on quarks
  WrappedTH1 *h_Hs_Bc_MassOf_bq;
  WrappedTH1 *h_Hs_Bc_bqq_Top_pT;
  WrappedTH1 *h_Hs_Bc_qq_Top_pT;
  WrappedTH1 *h_Hs_Bc_qq_bq_Top_pT;
  WrappedTH1 *h_Hs_Bc_bq_Top_pT;
  WrappedTH1 *h_Hs_Bc_bqq_W_pT;
  WrappedTH1 *h_Hs_Bc_qq_W_pT;
  WrappedTH1 *h_Hs_Bc_qq_bq_W_pT;
  WrappedTH1 *h_Hs_Bc_bq_W_pT;
  
  WrappedTH1 *h_Hs_bqcase_deltaR_bq;
  WrappedTH1 *h_Hs_bqcase_Prob_deltaR_lt;
  WrappedTH1 *h_Hs_BoostedW_deltaR_qq;
  WrappedTH1 *h_Hs_BoostedW_Prob_deltaR_lt;

  WrappedTH1 *h_Hs_otherBcloseToTopProd;
  WrappedTH1 *h_Hs_otherBcloseToTopProd_whichProd;
  WrappedTH1 *h_Hs_otherBclose_BoostedTop;
  WrappedTH1 *h_Hs_otherBclose_BoostedW;
   
  WrappedTH1 *h_Hs_MatchedWithDiquark_fj_NumOf_Subjets;
  WrappedTH1 *h_Hs_MatchedWithDiquark_fj_hasBsubjet;
  WrappedTH1 *h_Hs_MatchedWithDiquark_fj_CSV;
  WrappedTH2 *h_Hs_BoostedW_W_pT_Vs_Fatjet_pT;
  
  // another approach, first the matching with AK8 and after the dR between the objects
  WrappedTH1 *h_Hs_QuarksintoFatJetMultiplicity;
  WrappedTH1 *h_Hs_BoostedWinFatJet_dR_qq;
  WrappedTH1 *h_Hs_BoostedWinFatJet_Prob_deltaR_lt;
  
  WrappedTH1 *h_Hs_TopProdInFatJet_fatjet_pT;
  WrappedTH2 *h_Hs_TopProdInFatJet_Top_pT_Vs_fatjet_pT;
  WrappedTH1 *h_Hs_TopProdInFatJet_Higgs_pT;
  WrappedTH1 *h_Hs_TopProdInFatJet_hasBsubjet;
  WrappedTH1 *h_Hs_TopProdInFatJet_Njettinesstau1;
  WrappedTH1 *h_Hs_TopProdInFatJet_Njettinesstau2;
  WrappedTH1 *h_Hs_TopProdInFatJet_Njettinesstau3;
  WrappedTH1 *h_Hs_TopProdInFatJet_Njettinesstau4;
  WrappedTH1 *h_Hs_TopProdInFatJet_tau2DIVtau1;
  WrappedTH1 *h_Hs_TopProdInFatJet_tau3DIVtau2;
  WrappedTH2 *h_Hs_TopProdInFatJet_tau2DIVtau1_Vs_fatjet_pT;
  WrappedTH2 *h_Hs_TopProdInFatJet_tau3DIVtau2_Vs_fatjet_pT;
  WrappedTH1 *h_Hs_TopProdInFatJet_ldgORsubldg;
  WrappedTH1Triplet *h_Hs_TopProdInFatJet_TFtau32cut_fatjet_pT;
  WrappedTH1Triplet *h_Hs_TopProdInFatJet_TFtau21cut_fatjet_pT;
  WrappedTH2 *h_Hs_TopProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT;
  WrappedTH2 *h_Hs_TopProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT;

  WrappedTH1 *h_Hs_WProdInFatJet_fatjet_pT;
  WrappedTH2 *h_Hs_WProdInFatJet_W_pT_Vs_fatjet_pT;
  WrappedTH1 *h_Hs_WProdInFatJet_Higgs_pT;
  WrappedTH1 *h_Hs_WProdInFatJet_hasBsubjet;
  WrappedTH1 *h_Hs_WProdInFatJet_Njettinesstau1;
  WrappedTH1 *h_Hs_WProdInFatJet_Njettinesstau2;
  WrappedTH1 *h_Hs_WProdInFatJet_Njettinesstau3;
  WrappedTH1 *h_Hs_WProdInFatJet_Njettinesstau4;
  WrappedTH1 *h_Hs_WProdInFatJet_tau2DIVtau1;
  WrappedTH1 *h_Hs_WProdInFatJet_tau3DIVtau2;
  WrappedTH2 *h_Hs_WProdInFatJet_tau2DIVtau1_Vs_fatjet_pT;
  WrappedTH2 *h_Hs_WProdInFatJet_tau3DIVtau2_Vs_fatjet_pT;
  WrappedTH1 *h_Hs_WProdInFatJet_ldgORsubldg;
  WrappedTH1Triplet *h_Hs_WProdInFatJet_TFtau32cut_fatjet_pT;
  WrappedTH1Triplet *h_Hs_WProdInFatJet_TFtau21cut_fatjet_pT;
  WrappedTH2 *h_Hs_WProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT;
  WrappedTH2 *h_Hs_WProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT;

  WrappedTH1 *h_Hs_bqInFatJet_fatjet_pT;
  WrappedTH2 *h_Hs_bqInFatJet_Top_pT_Vs_fatjet_pT;
  WrappedTH2 *h_Hs_bqInFatJet_W_pT_Vs_fatjet_pT;
  WrappedTH1 *h_Hs_bqInFatJet_Higgs_pT;
  WrappedTH1 *h_Hs_bqInFatJet_hasBsubjet;
  WrappedTH1 *h_Hs_bqInFatJet_Njettinesstau1;
  WrappedTH1 *h_Hs_bqInFatJet_Njettinesstau2;
  WrappedTH1 *h_Hs_bqInFatJet_Njettinesstau3;
  WrappedTH1 *h_Hs_bqInFatJet_Njettinesstau4;
  WrappedTH1 *h_Hs_bqInFatJet_tau2DIVtau1;
  WrappedTH1 *h_Hs_bqInFatJet_tau3DIVtau2;
  WrappedTH2 *h_Hs_bqInFatJet_tau2DIVtau1_Vs_fatjet_pT;
  WrappedTH2 *h_Hs_bqInFatJet_tau3DIVtau2_Vs_fatjet_pT;
  WrappedTH1 *h_Hs_bqInFatJet_ldgORsubldg;
  WrappedTH1Triplet *h_Hs_bqInFatJet_TFtau32cut_fatjet_pT;
  WrappedTH1Triplet *h_Hs_bqInFatJet_TFtau21cut_fatjet_pT;
  WrappedTH2 *h_Hs_bqInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT;
  WrappedTH2 *h_Hs_bqInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT;

  // ---------------------------------------------------------------------------------------------------
                                                  
};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(BoostedTopologies);

BoostedTopologies::BoostedTopologies(const ParameterSet& config, const TH1* skimCounters)
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

void BoostedTopologies::book(TDirectory *dir) {

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
  TDirectory* myInclusiveDir   = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myInclusiveLabel);
  TDirectory* myBoolIsFalseDir = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myFalseLabel);
  TDirectory* myBoolIsTrueDir  = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myTrueLabel);
  std::vector<TDirectory*> myTripletDirs = {myInclusiveDir, myBoolIsFalseDir, myBoolIsTrueDir};
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Event Variables                                                                                      
  h_genMET_Et         =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital      , th1, "genMET_Et"    , ";Gen E_{T}^{miss} (GeV)"       , 100,  0.0,   +500.0);
  h_genMET_Phi        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "genMET_Phi"   , ";Gen E_{T}^{miss} #phi (rads)" , nBinsPhi, minPhi, maxPhi);
  h_genHT_GenJets     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital      , th1, "genHT_GenJets", ";GenJ H_{T} (GeV)"             , nBinsM  , minM  , maxM   );

  // GenJets                                                                                                                                  
  h_GenJets_N   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJets_N" , ";genJet multiplicity" , 30, 0.0, 30.0);


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


  // kchristo ---------Boosted Topologies ------------------------------------------------------------------------

  h_bfromH_Higgs_dR           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromH_Higgs_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_objectsfromHiggstop_maxdR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "objectsfromHiggstop_maxdR", ";#DeltaR_{max}", 50 , 0.0 , 5.0 );
  h_objectsfromHiggstop_mindR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "objectsfromHiggstop_mindR", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );
  h_objectsfromHiggstop_Prob_mindR_lt_p8 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "objectsfromHiggstop_Prob_mindR_lt_p8",
								      ";#DeltaR_{min}", 2 , 0.0, 2.0);
  h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than", ";p_{T} (GeV)", 3 , 0.0, 3.0);
  h_bfromHiggstop_Higgstop_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bfromHiggstop_Higgstop_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_b_top_fromhiggs_underdR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "b_top_fromhiggs_underdR"," ", 2 , 0.0, 2.0);
  h_b_higgs_underdR           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "b_higgs_underdR"," ", 2 , 0.0, 2.0);

  h_bfromH_Higgs_dR_Vs_Higgs_Pt              = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bfromH_Higgs_dR_Vs_Higgs_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_objectsfromHiggstop_maxdR_Vs_Higgstop_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "objectsfromHiggstop_maxdR_Vs_Higgstop_Pt", ";#DeltaR_{max};p_{T} (GeV)", nBinsdR,0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_objectsfromHiggstop_maxdR_Vs_Higgs_Pt    = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "objectsfromHiggstop_maxdR_Vs_Higgs_Pt",";#DeltaR_{max};p_{T} (GeV)", nBinsdR,0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_objectsfromHiggstop_mindR_Vs_Higgstop_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "objectsfromHiggstop_mindR_Vs_Higgstop_Pt", ";#DeltaR_{min};p_{T} (GeV)", nBinsdR,0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_objectsfromHiggstop_mindR_Vs_Higgs_Pt    = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "objectsfromHiggstop_mindR_Vs_Higgs_Pt",";#DeltaR_{min};p_{T} (GeV)", nBinsdR,0.0, 5.0, nBinsPt, 0.0, 1000.0); 
  h_bfromHiggstop_Higgstop_dR_Vs_Higgstop_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bfromHiggstop_Higgstop_dR_Vs_Higgstop_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  
  h_WfromHiggstop_closestb_dR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "WfromHiggstop_closestb_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_WfromHiggstop_samedecayb_dR= fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "WfromHiggstop_samedecayb_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_for_WbfromH_intofatjet_W_Pt   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "for_WbfromH_intofatjet_W_Pt", ";p_{t} (GeV)", 200 , 0.0, 1000.0 );
  h_for_WbfromH_intofatjet_Top_Pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "for_WbfromH_intofatjet_Top_Pt", ";p_{t} (GeV)", 200 , 0.0, 1000.0 );
  
  h_WfromHiggstop_closestb_dR_Vs_W_Pt   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "WfromHiggstop_closestb_dR_Vs_W_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_WfromHiggstop_samedecayb_dR_Vs_W_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "WfromHiggstop_samedecayb_dR_Vs_W_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_WfromHiggstop_closestb_dR_Vs_samedecaytop_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "WfromHiggstop_closestb_dR_Vs_samedecaytop_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_WfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "WfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_Probdecayproductsintofatjet_Hs           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Probdecayproductsintofatjet_Hs"," ", 4 , 0.0, 4.0);
  h_Hs_QuarksintofatjetMultiplicity          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksintofatjetMultiplicity", " ", 6 , 0.0 , 6.0);
  h_Hs_isbQuarkintofatjet                    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_isbQuarkintofatjet", " ", 2, 0.0 , 2.0);
  h_Pairsinbarycenter_enoughdeltaR_Hs        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Pairsinbarycenter_enoughdeltaR_Hs"," ", 4 , 0.0, 4.0);


  h_bNOTfromHiggstop_top_dR      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "bNOTfromHiggstop_top_dR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_objectsNOTfromHiggstop_maxdR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "objectsNOTfromHiggstop_maxdR", ";#DeltaR_{max}", 50 , 0.0 , 5.0 );
  h_objectsNOTfromHiggstop_maxdR_Vs_top_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "objectsNOTfromHiggstop_maxdR_Vs_top_Pt",";#DeltaR_{max};p_{T} (GeV)", nBinsdR,0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_objectsNOTfromHiggstop_mindR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "objectsNOTfromHiggstop_mindR", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );
  h_objectsNOTfromHiggstop_Prob_mindR_lt_p8 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, 
									 "objectsNOTfromHiggstop_Prob_mindR_lt_p8",";#DeltaR_{min}", 
									 2 , 0.0, 2.0);
  h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "objectsfromNOTHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than", ";p_{T} (GeV)", 3 , 0.0, 3.0);
  h_objectsNOTfromHiggstop_mindR_Vs_top_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "objectsNOTfromHiggstop_mindR_Vs_top_Pt",";#DeltaR_{min};p_{T} (GeV)", nBinsdR,0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_bNOTfromHiggstop_top_dR_Vs_top_Pt      = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bNOTfromHiggstop_top_dR_Vs_top_Pt",";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_b_topNOTfromhiggs_underdR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "b_topNOTfromhiggs_underdR"," ", 2 , 0.0, 2.0);  

  h_WNOTfromHiggstop_closestb_dR        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "WNOTfromHiggstop_closestb_dR", ";#DeltaR", 50 , 0.0 ,5.0 );
  h_WNOTfromHiggstop_samedecayb_dR      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "WNOTfromHiggstop_samedecayb_dR", ";#DeltaR", 50 , 0.0, 5.0 );
  h_WNOTfromHiggstop_closestb_dR_Vs_W_Pt   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "WNOTfromHiggstop_closestb_dR_Vs_W_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_WNOTfromHiggstop_samedecayb_dR_Vs_W_Pt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "WNOTfromHiggstop_samedecayb_dR_Vs_W_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_WNOTfromHiggstop_closestb_dR_Vs_samedecaytop_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "WNOTfromHiggstop_closestb_dR_Vs_samedecaytop_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_WNOTfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "WNOTfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, nBinsPt, 0.0, 1000.0);
  h_for_WbNOTfromH_intofatjet_W_Pt   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "for_WbNOTfromH_intofatjet_W_Pt", ";p_{t} (GeV)", 200 , 0.0, 1000.0 );
  h_for_WbNOTfromH_intofatjet_Top_Pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "for_WbNOTfromH_intofatjet_Top_Pt", ";p_{t} (GeV)", 200 , 0.0, 1000.0 );
  h_Probdecayproductsintofatjet_NOTHs           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Probdecayproductsintofatjet_NOTHs"," ", 4 , 0.0, 4.0);
  h_NotHs_QuarksintofatjetMultiplicity          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_QuarksintofatjetMultiplicity", " ", 6 , 0.0 , 6.0);
  h_NotHs_isbQuarkintofatjet                    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_isbQuarkintofatjet", " ", 2 , 0.0 , 2.0);
  h_Pairsinbarycenter_enoughdeltaR_NOTHs        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Pairsinbarycenter_enoughdeltaR_NOTHs"," ", 4 , 0.0, 4.0);
 
  // try with gen-jets +++++++++++++++++++++++++++++++++++
  h_Hs_deltaR_bQuark_closestbJet            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_deltaR_bQuark_closestbJet", ";#DeltaR_{min}", 30 , 0.0 , 3.0 );
  h_Hs_deltaR_bQuark_fromH_closestbJet      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_deltaR_bQuark_fromH_closestbJet", ";#DeltaR_{min}", 30 , 0.0 , 3.0 );
  h_Hs_deltaR_QfromW_closestJet             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_deltaR_QfromW_closestJet", ";#DeltaR_{min}", 30 , 0.0 , 3.0 );
  h_Hs_numof_obj_fromtop_matchedwith_genjet = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_numof_obj_fromtop_matchedwith_genjet", "Objects Matched with genJets", 4 , -0.5 , 3.5);
  h_Hs_ProbdecayJetsintofatjet              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_ProbdecayJetsintofatjet", " ", 4 , 0.0 , 4.0);
  h_Hs_JetsintofatjetMultiplicity           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_JetsintofatjetMultiplicity", " ", 6 , 0.0 , 6.0);
  h_Hs_isbJetintofatjet                     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_isbJetintofatjet", " ", 2 , 0.0 , 2.0);
  h_Hs_BjetInsideFatJet_Top_Pt              = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
										"Hs_BjetInsideFatJet_Top_Pt",
										";p_{T} (GeV)",
										200, 0.0, 1000.0);
  h_Hs_BjetInsideFatJet_W_Pt                = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
                                                                                "Hs_BjetInsideFatJet_W_Pt",
                                                                                ";p_{T} (GeV)",
                                                                                200, 0.0, 1000.0);

  h_Hs_BjetInsideFatJet_H_Pt                = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
									        "Hs_BjetInsideFatJet_H_Pt",
									        ";p_{T} (GeV)",
									        200, 0.0, 1000.0);


  h_NotHs_deltaR_bQuark_closestbJet            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_deltaR_bQuark_closestbJet", ";#DeltaR_{min}", 30 , 0.0 , 3.0 );
  h_NotHs_deltaR_QfromW_closestJet             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_deltaR_QfromW_closestJet", ";#DeltaR_{min}", 30 , 0.0 , 3.0 );
  h_NotHs_numof_obj_fromtop_matchedwith_genjet = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_numof_obj_fromtop_matchedwith_genjet", "Objects Matched with genJets", 4 , -0.5 , 3.5);
  h_NotHs_ProbdecayJetsintofatjet              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_ProbdecayJetsintofatjet", " ", 4 , 0.0 , 4.0);
  h_NotHs_JetsintofatjetMultiplicity           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_JetsintofatjetMultiplicity", " ", 6 , 0.0 , 6.0);
  h_NotHs_isbJetintofatjet                     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_isbJetintofatjet", " ", 2 , 0.0 , 2.0);

  h_NotHs_BjetInsideFatJet_Top_Pt              = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
                                                                                "NotHs_BjetInsideFatJet_Top_Pt",
                                                                                ";p_{T} (GeV)",
                                                                                200, 0.0, 1000.0);
  h_NotHs_BjetInsideFatJet_W_Pt                = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
                                                                                "NotHs_BjetInsideFatJet_W_Pt",
                                                                                ";p_{T} (GeV)",
                                                                                200, 0.0, 1000.0);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++

  h_Hs_objectsfromtop_top_mindR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_objectsfromtop_top_mindR", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );
  h_Hs_objectsfromtop_top_maxdR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_objectsfromtop_top_maxdR", ";#DeltaR_{max}", 50 , 0.0 , 5.0 );
  h_Hs_which_objectfromtop_maxdR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_which_objectfromtop_maxdR"," ", 2 , 0.0, 2.0);
  h_NotHs_objectsfromtop_top_maxdR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_objectsfromtop_top_maxdR", ";#DeltaR_{max}",  50 , 0.0 , 5.0 );
  h_NotHs_objectsfromtop_top_mindR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_objectsfromtop_top_mindR", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );
  h_NotHs_which_objectfromtop_maxdR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_which_objectfromtop_maxdR"," ", 2 , 0.0, 2.0);

  h_Iffatjet_Hs_Top_pT         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_Top_pT", ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Iffatjet_Hs_W_pT           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_W_pT", ";p_{t} (GeV)", 100 , 0.0, 1000.0);  
  h_Iffatjet_Hs_Top_pT_Vs_W_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Iffatjet_Hs_Top_pT_Vs_W_pT", ";p_{T} (GeV);p_{T} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  h_Iffatjet_Hs_bfromH_pT      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromH_pT", ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Iffatjet_Hs_EventsWithHighTop_pT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_EventsWithHighTop_pT"," ", 2 , 0.0, 2.0);
  h_Iffatjet_Hs_EventsWithHighW_pT   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_EventsWithHighW_pT"," ", 2 , 0.0, 2.0);
  h_Iffatjet_Hs_HighWandTop_pT_dRmin = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_HighWandTop_pT_dRmin", ";#DeltaR_{min}", 20 , 0.0 , 2.0 );
  h_Iffatjet_Hs_bfromT_Top_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_Top_dEta", ";#Delta#eta", 50 , 0.0 ,5.0 );
  h_Iffatjet_Hs_bfromT_Top_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_Top_dPhi", ";#Delta#phi", 50 , 0.0 ,5.0 );
  h_Iffatjet_Hs_bfromT_Top_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_Top_dR", ";#DeltaR", 50 , 0.0 ,5.0 );
  h_Iffatjet_Hs_bfromT_W_dEta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_W_dEta", ";#Delta#eta", 50 ,0.0 ,5.0 );
  h_Iffatjet_Hs_bfromT_W_dPhi   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_W_dPhi", ";#Delta#phi", 50 ,0.0 ,5.0 );
  h_Iffatjet_Hs_bfromT_W_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_W_dR", ";#DeltaR", 50 , 0.0 ,5.0 );
  h_Iffatjet_Hs_bfromT_bfromH_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_bfromH_dEta", ";#Delta#eta", 50 ,0.0 ,5.0 );
  h_Iffatjet_Hs_bfromT_bfromH_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_bfromH_dPhi", ";#Delta#phi", 50 ,0.0 ,5.0 );
  h_Iffatjet_Hs_bfromT_bfromH_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_Hs_bfromT_bfromH_dR", ";#DeltaR", 50 , 0.0 ,5.0 );
  
  h_Iffatjet_NotHs_Top_pT         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_Top_pT", ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Iffatjet_NotHs_W_pT           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_W_pT", ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Iffatjet_NotHs_bfromH_pT      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromH_pT", ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Iffatjet_NotHs_Top_pT_Vs_W_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Iffatjet_NotHs_Top_pT_Vs_W_pT", ";p_{T} (GeV);p_{T} (GeV)", 200, 0.0, 1000.0, 200, 0.0, 1000.0);
  h_Iffatjet_NotHs_EventsWithHighTop_pT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_EventsWithHighTop_pT"," ", 2 , 0.0, 2.0);
  h_Iffatjet_NotHs_EventsWithHighW_pT   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_EventsWithHighW_pT"," ", 2 , 0.0, 2.0);
  h_Iffatjet_NotHs_HighWandTop_pT_dRmin = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_HighWandTop_pT_dRmin", ";#DeltaR_{min}", 20 , 0.0 , 2.0 );
  h_Iffatjet_NotHs_bfromT_Top_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_Top_dEta", ";#Delta#eta", 50 , 0.0 ,5.0 );
  h_Iffatjet_NotHs_bfromT_Top_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_Top_dPhi", ";#Delta#phi", 50 , 0.0 ,5.0 );
  h_Iffatjet_NotHs_bfromT_Top_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_Top_dR", ";#DeltaR", 50 , 0.0 ,5.0 );
  h_Iffatjet_NotHs_bfromT_W_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_W_dEta", ";#Delta#eta", 50 ,0.0 ,5.0 );
  h_Iffatjet_NotHs_bfromT_W_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_W_dPhi", ";#Delta#phi", 50 ,0.0 ,5.0 );
  h_Iffatjet_NotHs_bfromT_W_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_W_dR", ";#DeltaR", 50 , 0.0 ,5.0 );
  h_Iffatjet_NotHs_bfromT_bfromH_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_bfromH_dEta", ";#Delta#eta", 50 ,0.0 ,5.0 );
  h_Iffatjet_NotHs_bfromT_bfromH_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_bfromH_dPhi", ";#Delta#phi", 50 ,0.0 ,5.0 );
  h_Iffatjet_NotHs_bfromT_bfromH_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Iffatjet_NotHs_bfromT_bfromH_dR", ";#DeltaR", 50 , 0.0 ,5.0 );

  h_Hs_mostdistantfromtop_isb__top_pT          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__top_pT" , ";p_{t} (GeV)", 100 , 0.0, 1000.0); 
  h_Hs_mostdistantfromtop_isb__W_pT            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__W_pT" , ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_mostdistantfromtop_isb__b_pT            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__b_pT" , ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_mostdistantfromtop_isb__dRmax_b_top     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRmax_b_top", ";#DeltaR",  50 , 0.0 , 5.0 );
  h_Hs_mostdistantfromtop_isb_dRmin_b_objfromW = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb_dRmin_b_objfromW", ";#DeltaR_{min}",  50 , 0.0 , 5.0 );
  h_Hs_mostdistantfromtop_isb__dRqq_Vs_W_pT    = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_mostdistantfromtop_isb__dRqq_Vs_W_pT", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, 100, 0.0, 1000.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than", ";p_{T} (GeV)", 3 , 0.0, 3.0);

  h_Hs_mostdistantfromtop_isb__dRqq_Vs_top_pT  = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_mostdistantfromtop_isb__dRqq_Vs_top_pT", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, 100, 0.0, 1000.0);
  
  h_NotHs_mostdistantfromtop_isb__top_pT          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__top_pT" , ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__W_pT            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__W_pT" , ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__b_pT            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__b_pT" , ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__dRmax_b_top     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRmax_b_top", ";#DeltaR",  50 , 0.0 , 5.0 );
  h_NotHs_mostdistantfromtop_isb_dRmin_b_objfromW = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb_dRmin_b_objfromW", ";#DeltaR_{min}",  50 , 0.0 , 5.0 );
  h_NotHs_mostdistantfromtop_isb__dRqq_Vs_W_pT   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "NotHs_mostdistantfromtop_isb__dRqq_Vs_W_pT", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, 100, 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than", ";p_{T} (GeV)", 3 , 0.0, 3.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than", ";p_{T} (GeV)", 3 , 0.0, 3.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7_tf = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7_tf", " ", 2 , 0.0, 2.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than", ";p_{T} (GeV)", 3 , 0.0, 3.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf", " ", 2 , 0.0, 2.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_Vs_top_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "NotHs_mostdistantfromtop_isb__dRqq_Vs_top_pT", ";#DeltaR;p_{t} (GeV)", nBinsdR, 0.0, 5.0, 100, 0.0, 1000.0);

  //Study Boosted Topologies With fat Jets --------------------------------------------------------------------------------------- 
  
  h_Hs_QuarksFromW_deltaR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksFromW_deltaR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_Hs_QuarksFromW_Prob_deltaR          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksFromW_Prob_deltaR", ";#DeltaR" , 2 , 0.0 , 2.0 );
  h_Hs_QuarksintoBaryCenterMultiplicity = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksintoBaryCenterMultiplicity", " ", 5 , 0.0 , 5.0 );
  h_Hs_isbQuarkintoBaryCenter = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_isbQuarkintoBaryCenter", " " , 3 , 0.0 , 3.0 );
  h_Hs_OnlyQQ_dR_less_p7      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_OnlyQQ_dR_less_p7", " " , 2 , 0.0 , 2.0 );
  h_Hs_Prob_Diquark_match_with_fj = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Prob_Diquark_match_with_fj", " " , 2 , 0.0 , 2.0 );
  h_Hs_MatchedWithDiquark_fj_pT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_MatchedWithDiquark_fj_pT",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  
  h_Hs_objectsfromtop_dRmax_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_objectsfromtop_dRmax_pTcuts", ";#DeltaR_{max}", 50 , 0.0 , 5.0 );
  h_Hs_objectsfromtop_dRmin_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_objectsfromtop_dRmin_pTcuts", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );
  h_Hs_objectsfromtop_Prob_dRmax_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_objectsfromtop_Prob_dRmax_pTcuts", ";#DeltaR_{max}", 2 , 0.0 , 2.0 );
  h_Hs_objectsfromtop_Prob_dRmin_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_objectsfromtop_Prob_dRmin_pTcuts", ";#DeltaR_{min}", 2 , 0.0 , 2.0 );
  h_Hs_QuarksFromW_deltaR_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksFromW_deltaR_pTcuts", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_Hs_QuarksFromW_Prob_deltaR_pTcuts  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksFromW_Prob_deltaR_pTcuts", ";#DeltaR" ,2 , 0.0 , 2.0 );
  h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksintoBaryCenterMultiplicity_pTcuts", " ", 5 , 0.0 , 5.0 );
  h_Hs_isbQuarkintoBaryCenter_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_isbQuarkintoBaryCenter_pTcuts", " " , 3 , 0.0 , 3.0 );
  h_Hs_OnlyQQ_dR_less_p7_pTcuts      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_OnlyQQ_dR_less_p7_pTcuts", " " , 2 , 0.0 , 2.0 );
  h_Hs_Prob_Diquark_match_with_fj_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Prob_Diquark_match_with_fj_pTcuts", " " , 2 , 0.0 ,2.0 );
  h_Hs_MatchedWithDiquark_fj_pT_pTcuts   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_MatchedWithDiquark_fj_pT_pTcuts",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_MatchedWithDiquark_Prob_fj_pT_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_MatchedWithDiquark_Prob_fj_pT_pTcuts", " " , 2, 0.0 ,2.0 );
  h_Hs_MatchedWithDiquark_Prob_fj_pT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_MatchedWithDiquark_Prob_fj_pT", " ", 2, 0.0 ,2.0 );

  h_Hs_QuarksFromTop_Passed_pTcuts       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksFromTop_Passed_pTcuts", " " , 2 , 0.0 ,2.0 );

  h_NotHs_QuarksFromW_deltaR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_QuarksFromW_deltaR", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_NotHs_QuarksFromW_Prob_deltaR          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_QuarksFromW_Prob_deltaR", ";#DeltaR" ,2 , 0.0 , 2.0 );
  h_NotHs_QuarksintoBaryCenterMultiplicity = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_QuarksintoBaryCenterMultiplicity", " ", 5 , 0.0 , 5.0 );
  h_NotHs_isbQuarkintoBaryCenter = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_isbQuarkintoBaryCenter", " " , 3 , 0.0 , 3.0 );
  h_NotHs_OnlyQQ_dR_less_p7      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_OnlyQQ_dR_less_p7", " " , 2 , 0.0 , 2.0 );
  h_NotHs_Prob_Diquark_match_with_fj = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_Prob_Diquark_match_with_fj", " " , 2 , 0.0 ,2.0 );
  h_NotHs_MatchedWithDiquark_fj_pT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_MatchedWithDiquark_fj_pT",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  
  h_NotHs_objectsfromtop_dRmax_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_objectsfromtop_dRmax_pTcuts", ";#DeltaR_{max}", 50 , 0.0 , 5.0 );
  h_NotHs_objectsfromtop_dRmin_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_objectsfromtop_dRmin_pTcuts", ";#DeltaR_{min}", 50 , 0.0 , 5.0 );
  h_NotHs_objectsfromtop_Prob_dRmax_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_objectsfromtop_Prob_dRmax_pTcuts", ";#DeltaR_{max}", 2 , 0.0 , 2.0 );
  h_NotHs_objectsfromtop_Prob_dRmin_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_objectsfromtop_Prob_dRmin_pTcuts", ";#DeltaR_{min}", 2 , 0.0 , 2.0 );
  h_NotHs_QuarksFromW_deltaR_pTcuts   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_QuarksFromW_deltaR_pTcuts", ";#DeltaR", 50 , 0.0 , 5.0 );
  h_NotHs_QuarksFromW_Prob_deltaR_pTcuts  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_QuarksFromW_Prob_deltaR_pTcuts", ";#DeltaR" ,2 , 0.0 , 2.0 );
  h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_QuarksintoBaryCenterMultiplicity_pTcuts", " ", 5 , 0.0 , 5.0 );
  h_NotHs_isbQuarkintoBaryCenter_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_isbQuarkintoBaryCenter_pTcuts", " " , 3 ,0.0 , 3.0 );
  h_NotHs_OnlyQQ_dR_less_p7_pTcuts      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_OnlyQQ_dR_less_p7_pTcuts", " " , 2 , 0.0 , 2.0 );
  h_NotHs_Prob_Diquark_match_with_fj_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_Prob_Diquark_match_with_fj_pTcuts", " " , 2 , 0.0 ,2.0 );
  h_NotHs_MatchedWithDiquark_fj_pT_pTcuts   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_MatchedWithDiquark_fj_pT_pTcuts",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_MatchedWithDiquark_Prob_fj_pT_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_MatchedWithDiquark_Prob_fj_pT_pTcuts", " " , 2, 0.0 ,2.0 );
  h_NotHs_MatchedWithDiquark_Prob_fj_pT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_MatchedWithDiquark_Prob_fj_pT", " " , 2, 0.0 ,2.0 );
  h_NotHs_QuarksFromTop_Passed_pTcuts       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_QuarksFromTop_Passed_pTcuts", " " , 2, 0.0 ,2.0 );



  h_Hs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1,
										     "Hs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj",
										     " ", 2 , 0.0, 2.0);
  h_Hs_objectsfromtop_mindR_ltp8_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, 
									   "Hs_objectsfromtop_mindR_ltp8_matchedfj_pt" , 
									   ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj"," ", 2 , 0.0, 2.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj"," ", 2 , 0.0, 2.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);


  h_NotHs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1,
                                                                                     "NotHs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj",
                                                                                     " ", 2 , 0.0, 2.0);
  h_NotHs_objectsfromtop_mindR_ltp8_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1,
                                                                           "NotHs_objectsfromtop_mindR_ltp8_matchedfj_pt" ,
                                                                           ";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);

  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj"," ", 2 , 0.0, 2.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj"," ", 2 , 0.0, 2.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);
  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotHs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt",";p_{t} (GeV)", 100 , 0.0, 1000.0);

  // update 24.01.2018 after ptcuts
  h_Hs_Bc_MassOf_bq    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_MassOf_bq",";m_{bq} (GeV)", 70, 0.0, 350.0);
  h_Hs_Bc_bqq_Top_pT   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_bqq_Top_pT",";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_Bc_qq_Top_pT    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_qq_Top_pT",";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_Bc_qq_bq_Top_pT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_qq_bq_Top_pT",";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_Bc_bq_Top_pT    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_bq_Top_pT",";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_Bc_bqq_W_pT   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_bqq_W_pT",";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_Bc_qq_W_pT    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_qq_W_pT",";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_Bc_qq_bq_W_pT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_qq_bq_W_pT",";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_Bc_bq_W_pT    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_Bc_bq_W_pT",";p_{T} (GeV)", 100 , 0.0, 1000.0);
  
  h_Hs_bqcase_deltaR_bq  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqcase_deltaR_bq", ";#DeltaR_{qq}", 10 , 0.0 , 1.0 );
  h_Hs_bqcase_Prob_deltaR_lt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqcase_Prob_deltaR_lt", ";#DeltaR_{qq}", 2 , 0.0 , 2.0 );
  h_Hs_BoostedW_deltaR_qq  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_BoostedW_deltaR_qq", ";#DeltaR_{bq}", 10 , 0.0 , 1.0 );
  h_Hs_BoostedW_Prob_deltaR_lt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_BoostedW_Prob_deltaR_lt", ";#DeltaR_{bq}", 2 , 0.0, 2.0 );

  h_Hs_otherBcloseToTopProd = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_otherBcloseToTopProd",";Other b close", 4, -0.5, 3.5);
  h_Hs_otherBcloseToTopProd_whichProd = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_otherBcloseToTopProd_whichProd"," ", 2, 0.0, 2.0);
  h_Hs_otherBclose_BoostedW = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_otherBclose_BoostedW"," ", 2, 0.0, 2.0);
  h_Hs_otherBclose_BoostedTop = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_otherBclose_BoostedTop"," ", 2, 0.0, 2.0);

  h_Hs_MatchedWithDiquark_fj_NumOf_Subjets = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_MatchedWithDiquark_fj_NumOf_Subjets",";Subjets", 4, -0.5, 3.5);
  h_Hs_MatchedWithDiquark_fj_hasBsubjet    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_MatchedWithDiquark_fj_hasBsubjet",";Has b-subjet", 2, -0.5, 1.5);
  h_Hs_MatchedWithDiquark_fj_CSV           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_MatchedWithDiquark_fj_CSV",";CSV", 100, 0.0, 1.0);
  h_Hs_BoostedW_W_pT_Vs_Fatjet_pT          = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_BoostedW_W_pT_Vs_Fatjet_pT"\
									, ";p_{t} (GeV);p_{t} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  // another approach, first the matching with AK8 and after the dR between the objects
  h_Hs_QuarksintoFatJetMultiplicity = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksintoFatJetMultiplicity"," ", 4, 0.0, 4.0);
  h_Hs_BoostedWinFatJet_dR_qq       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_BoostedWinFatJet_dR_qq", ";#DeltaR_{qq}", 20 , 0.0 , 2.0 );
  h_Hs_BoostedWinFatJet_Prob_deltaR_lt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_BoostedWinFatJet_Prob_deltaR_lt", ";#DeltaR_{qq}", 2 , 0.0, 2.0 );
  
  h_Hs_TopProdInFatJet_fatjet_pT           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_fatjet_pT",
									";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_TopProdInFatJet_Top_pT_Vs_fatjet_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_TopProdInFatJet_Top_pT_Vs_fatjet_pT" \
                                                                        , ";p_{T,jet} (GeV);p_{T,top} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  h_Hs_TopProdInFatJet_Higgs_pT            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_Higgs_pT",
									";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_TopProdInFatJet_hasBsubjet          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_hasBsubjet",
									";Has b-subjet", 2 , 0.0, 2.0 );
  h_Hs_TopProdInFatJet_Njettinesstau1      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_Njettinesstau1",
                                                                        ";#tau_{1}", 50 , 0.0, 1.0);
  h_Hs_TopProdInFatJet_Njettinesstau2      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_Njettinesstau2",
								        ";#tau_{2}", 50 , 0.0, 1.0);
  h_Hs_TopProdInFatJet_Njettinesstau3      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_Njettinesstau3",
								        ";#tau_{3}", 50 , 0.0, 1.0);
  h_Hs_TopProdInFatJet_Njettinesstau4      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_Njettinesstau4",
									";#tau_{4}", 50 , 0.0, 1.0);
  h_Hs_TopProdInFatJet_tau2DIVtau1         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_tau2DIVtau1",
                                                                        ";#tau_{2} / #tau_{1}", 50 , 0.0, 1.0);
  h_Hs_TopProdInFatJet_tau3DIVtau2         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_tau3DIVtau2",
								      ";#tau_{3} / #tau_{2}", 50 , 0.0, 1.0);
  h_Hs_TopProdInFatJet_tau2DIVtau1_Vs_fatjet_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_TopProdInFatJet_tau2DIVtau1_Vs_fatjet_pT", ";p_{t,jet} (GeV);#tau_{2} / #tau_{1}", 100, 0.0, 1000.0, 100, 0.0, 1.0);
  h_Hs_TopProdInFatJet_tau3DIVtau2_Vs_fatjet_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_TopProdInFatJet_tau3DIVtau2_Vs_fatjet_pT", ";p_{t,jet} (GeV);#tau_{3} / #tau_{2}", 100, 0.0, 1000.0, 100, 0.0, 1.0); 
  h_Hs_TopProdInFatJet_ldgORsubldg = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_TopProdInFatJet_ldgORsubldg"," ", 3, 0.0, 3.0);
  h_Hs_TopProdInFatJet_TFtau32cut_fatjet_pT      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
										     "Hs_TopProdInFatJet_TFtau32cut_fatjet_pT",
										     ";p_{T,jet} (GeV)",
										     100, 0.0, 1000.0);
  h_Hs_TopProdInFatJet_TFtau21cut_fatjet_pT      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
                                                                                     "Hs_TopProdInFatJet_TFtau21cut_fatjet_pT",
                                                                                     ";p_{T,jet} (GeV)",
                                                                                     100, 0.0, 1000.0);
  h_Hs_TopProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT= fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_TopProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT",
										      ";p_{T,subldg} (GeV);p_{T,BoostTop} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  h_Hs_TopProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT= fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_TopProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT",
										      ";p_{T,ldg} (GeV);p_{T,BoostTop} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);


  h_Hs_WProdInFatJet_fatjet_pT           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_fatjet_pT",
								      ";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_WProdInFatJet_W_pT_Vs_fatjet_pT   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_WProdInFatJet_W_pT_Vs_fatjet_pT" \
								      , ";p_{T,jet} (GeV);p_{T,W} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  h_Hs_WProdInFatJet_Higgs_pT            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_Higgs_pT",
								      ";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_WProdInFatJet_hasBsubjet          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_hasBsubjet",
								      ";Has b-subjet", 2 , 0.0, 2.0 );
  h_Hs_WProdInFatJet_Njettinesstau1      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_Njettinesstau1",
								      ";#tau_{1}", 50 , 0.0, 1.0);
  h_Hs_WProdInFatJet_Njettinesstau2      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_Njettinesstau2",
								      ";#tau_{2}", 50 , 0.0, 1.0);
  h_Hs_WProdInFatJet_Njettinesstau3      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_Njettinesstau3",
								      ";#tau_{3}", 50 , 0.0, 1.0);
  h_Hs_WProdInFatJet_Njettinesstau4      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_Njettinesstau4",
								      ";#tau_{4}", 50 , 0.0, 1.0);
  h_Hs_WProdInFatJet_tau2DIVtau1         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_tau2DIVtau1",
								      ";#tau_{2} / #tau_{1}", 50 , 0.0, 1.0);
  h_Hs_WProdInFatJet_tau3DIVtau2         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_tau3DIVtau2",
								      ";#tau_{3} / #tau_{2}", 50 , 0.0, 1.0);
  h_Hs_WProdInFatJet_tau2DIVtau1_Vs_fatjet_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_WProdInFatJet_tau2DIVtau1_Vs_fatjet_pT", ";p_{T,jet} (GeV);#tau_{2} / #tau_{1}", 100, 0.0, 1000.0, 100, 0.0, 1.0);
  h_Hs_WProdInFatJet_tau3DIVtau2_Vs_fatjet_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_WProdInFatJet_tau3DIVtau2_Vs_fatjet_pT", ";p_{T,jet} (GeV);#tau_{3} / #tau_{2}", 100, 0.0, 1000.0, 100, 0.0, 1.0);
  h_Hs_WProdInFatJet_ldgORsubldg = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_WProdInFatJet_ldgORsubldg"," ", 3, 0.0, 3.0);
  h_Hs_WProdInFatJet_TFtau32cut_fatjet_pT      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
                                                                                   "Hs_WProdInFatJet_TFtau32cut_fatjet_pT",
                                                                                   ";p_{T,jet} (GeV)",
                                                                                   100, 0.0, 1000.0);
  h_Hs_WProdInFatJet_TFtau21cut_fatjet_pT      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
                                                                                   "Hs_WProdInFatJet_TFtau21cut_fatjet_pT",
                                                                                   ";p_{T,jet} (GeV)",
                                                                                   100, 0.0, 1000.0);
  h_Hs_WProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT= fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_WProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT",
										    ";p_{T,subldg} (GeV);p_{T,BoostW} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  h_Hs_WProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT= fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_WProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT",
										    ";p_{T,ldg} (GeV);p_{T,BoostW} (GeV)",  100, 0.0, 1000.0, 100, 0.0, 1000.0);



  h_Hs_bqInFatJet_fatjet_pT           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_fatjet_pT",
                                                                      ";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_bqInFatJet_Top_pT_Vs_fatjet_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_bqInFatJet_Top_pT_Vs_fatjet_pT" \
                                                                      , ";p_{T,jet} (GeV);p_{T,top} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  h_Hs_bqInFatJet_W_pT_Vs_fatjet_pT   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_bqInFatJet_W_pT_Vs_fatjet_pT" \
								   , ";p_{T,jet} (GeV);p_{T,W} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  h_Hs_bqInFatJet_Higgs_pT            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_Higgs_pT",
                                                                      ";p_{T} (GeV)", 100 , 0.0, 1000.0);
  h_Hs_bqInFatJet_hasBsubjet          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_hasBsubjet",
                                                                      ";Has b-subjet", 2 , 0.0, 2.0 );
  h_Hs_bqInFatJet_Njettinesstau1      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_Njettinesstau1",
                                                                      ";#tau_{1}", 50 , 0.0, 1.0);
  h_Hs_bqInFatJet_Njettinesstau2      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_Njettinesstau2",
                                                                      ";#tau_{2}", 50 , 0.0, 1.0);
  h_Hs_bqInFatJet_Njettinesstau3      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_Njettinesstau3",
                                                                      ";#tau_{3}", 50 , 0.0, 1.0);
  h_Hs_bqInFatJet_Njettinesstau4      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_Njettinesstau4",
                                                                      ";#tau_{4}", 50 , 0.0, 1.0);
  h_Hs_bqInFatJet_tau2DIVtau1         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_tau2DIVtau1",
                                                                      ";#tau_{2} / #tau_{1}", 50 , 0.0, 1.0);
  h_Hs_bqInFatJet_tau3DIVtau2         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_tau3DIVtau2",
                                                                      ";#tau_{3} / #tau_{2}", 50 , 0.0, 1.0);
  h_Hs_bqInFatJet_tau2DIVtau1_Vs_fatjet_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_bqInFatJet_tau2DIVtau1_Vs_fatjet_pT", ";p_{T,jet} (GeV);#tau_{2} / #tau_{1}", 100, 0.0, 1000.0, 100, 0.0, 1.0);
  h_Hs_bqInFatJet_tau3DIVtau2_Vs_fatjet_pT = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_bqInFatJet_tau3DIVtau2_Vs_fatjet_pT", ";p_{T,,jet} (GeV);#tau_{3} / #tau_{2}", 100, 0.0, 1000.0, 100, 0.0, 1.0);
  h_Hs_bqInFatJet_ldgORsubldg = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_bqInFatJet_ldgORsubldg"," ", 3, 0.0, 3.0);
  h_Hs_bqInFatJet_TFtau32cut_fatjet_pT      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
										"Hs_bqInFatJet_TFtau32cut_fatjet_pT",
										";p_{T,jet} (GeV)",
										100, 0.0, 1000.0);
  h_Hs_bqInFatJet_TFtau21cut_fatjet_pT      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myTripletDirs,
                                                                                "Hs_bqInFatJet_TFtau21cut_fatjet_pT",
                                                                                ";p_{T,jet} (GeV)",
                                                                                100, 0.0, 1000.0);
  h_Hs_bqInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT= fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_bqInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT",
                                                                                    ";p_{T,subldg} (GeV);p_{T,bq} (GeV)", 100, 0.0, 1000.0, 100, 0.0, 1000.0);
  h_Hs_bqInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT= fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Hs_bqInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT",
                                                                                    ";p_{T,ldg} (GeV);p_{T,bq} (GeV)",  100, 0.0, 1000.0, 100, 0.0, 1000.0);

  // -----------------------------------------------------------------------------------------------------------------------------------
  // end of plots

  return;
}

void BoostedTopologies::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
}


void BoostedTopologies::process(Long64_t entry) {

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
  // useful vector with the p4 of the B jets
  
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
  // Event-Shape Variables
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

   // -----------------------------------------------------------------------------                                                              // Assign mother/daughers for Lxy & d0 calculation       
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
   // if (!p.isHardProcess()) continue;
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

 //////kchristo////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 ////////////////////////////////////////////////////// Study Boosted Topologies///////////////////////////////////////////////////////////// 
 //////////////old ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 if(0)
   {
     math::XYZTLorentzVector b_fromH_p4(0,0,0,0);
     for (auto& p: fEvent.genparticles().getGenParticles())
       {
	 if(std::abs(p.pdgId()) == 37 && p.isFirstCopy()) // find the Higgs 
	   {
	     genParticle lastHiggs;
	     unsigned int lastHiggscopy = GetTheLastCopy(p.index()); // we need to find the last copy in order to find the real daughters
	     lastHiggs =  fEvent.genparticles().getGenParticles()[lastHiggscopy];
	     std::vector<short> Higgs_d = lastHiggs.daughters();
	     for (unsigned int i = 0; i < Higgs_d.size() ; i++)
	       {
		 genParticle d = fEvent.genparticles().getGenParticles()[Higgs_d.at(i)];
		 if(std::abs(d.pdgId()) == 5) // find the b from Higgs
		   {
		     b_fromH_p4 = d.p4(); 
		   }
	       } // for Higgs daus
	   }     // if Higgs
       }         // for gen Particles
     
     for (auto& p: fEvent.genparticles().getGenParticles())
       {
	 
	 if(std::abs(p.pdgId()) == 6 && p.isFirstCopy()) // find the top
	   {
	     std::vector<short> top_mothers = p.mothers();
             //for (unsigned int j = 0; j < top_mothers.size() ; j++)
	     genParticle m = fEvent.genparticles().getGenParticles()[top_mothers.at(0)]; //create the particle object        

	     // try with gen-jets ++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                                                                                       
	     math::XYZTLorentzVector matched_genbJet_p4(0,0,0,0), matched_genbJet_fromH_p4(0,0,0,0);
	     std::vector<math::XYZTLorentzVector> genJets_from_W_p4(2);
	     math::XYZTLorentzVector p4_initializer(0,0,0,0);
	     for(size_t gj=0; gj < genJets_from_W_p4.size(); gj++)
	       {
		 genJets_from_W_p4.at(gj) = p4_initializer;
	       }
	     unsigned int numof_obj_fromtop_matchedwith_genjet = 0;
	     bool bfromHiggsMatched = false;
	     bool bJetInsideFatJet  = false;
	     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	     if(std::abs(m.pdgId()) == 37)//if the top comes from Higgs///////////////////////////////////////////////////////////////////////////
	       {
		 math::XYZTLorentzVector b_fromTop_fromH_p4(0,0,0,0), obj1_fromW_fromTop_p4(0,0,0,0), obj2_fromW_fromTop_p4(0,0,0,0);
		 math::XYZTLorentzVector W_fromTop_fromH_p4(0,0,0,0), Top_fromH_p4(0,0,0,0);
		 math::XYZTLorentzVector W_genReco_fromTop_fromH_p4(0,0,0,0), Top_genReco_fromH_p4(0,0,0,0), Higgs_genReco_p4(0,0,0,0);
		 bool isTheRightTop    = true;
		 bool obj1_exist       = false;
		 bool bQuarkIntoFatJet = false;

		 std::vector<short> Higgs_daughters = m.daughters();
		 for (unsigned int i = 0; i < Higgs_daughters.size() ; i++) 
		   {
		     genParticle d = fEvent.genparticles().getGenParticles()[Higgs_daughters.at(i)]; //create the particle object
		     if(std::abs(d.pdgId()) == 5) // find the b from Higgs
		       {
			 double deltaR_bfromH_Higgs = ROOT::Math::VectorUtil::DeltaR(d.p4(), m.p4()); //find dR(bfromH,H)
			 h_bfromH_Higgs_dR             -> Fill(deltaR_bfromH_Higgs);

			 h_b_higgs_underdR             -> Fill("< 0.8",0);  //just to determinate the label of the first bin
			 if(deltaR_bfromH_Higgs < 0.8)  h_b_higgs_underdR   -> Fill("< 0.8",1);
			 else                           h_b_higgs_underdR   -> Fill("> 0.8",1);

			 h_bfromH_Higgs_dR_Vs_Higgs_Pt -> Fill(deltaR_bfromH_Higgs,m.pt());
			 
		       }// find the b from Higgs
		   }//for Higgs daus 

		 Top_fromH_p4 = p.p4();
		 
		 unsigned int lastcopy = GetTheLastCopy(p.index()); // we need to find the last copy in order to find the real daughters     
		 genParticle lastT; //create the particle object of the last copy of top                        
		 lastT =  fEvent.genparticles().getGenParticles()[lastcopy];
		 std::vector<short> lastT_daughters = lastT.daughters();
		 
		 for (unsigned int i = 0; i < lastT_daughters.size() ; i++)
		   {
		     genParticle d; //create the particle object                                                                       
		     d =  fEvent.genparticles().getGenParticles()[lastT_daughters.at(i)];
		     
		     if (std::abs(d.pdgId()) == 24) // if W from top
		       {          
			 W_fromTop_fromH_p4 = d.p4();
			 genParticle lastW;
			 unsigned int lastWcopy = GetTheLastCopy(d.index()); // we need to find the last copy in order to find the real daughters 
			 lastW = fEvent.genparticles().getGenParticles()[lastWcopy];
			 std::vector<short> W_daughters = lastW.daughters();
			 for (size_t k = 0; k < W_daughters.size() ; k++)
			   {
			     genParticle top_grand_d; //create the particle object, i define this as grand_d, the daus of the W, so the grand-daughters of Top                                                                                              
			     top_grand_d =  fEvent.genparticles().getGenParticles()[W_daughters.at(k)];
			     
			     if ( !mcTools.IsQuark(top_grand_d.pdgId()) ) //if the daus of W are not quarks (f.e. leptons), we reject the top 
			       {
				 isTheRightTop = false;
				 continue;
			       }

			     if (obj1_exist == false)
			       {
				 obj1_exist = true;
				 obj1_fromW_fromTop_p4 = top_grand_d.p4();
			       }
			     else
			       {
				 obj2_fromW_fromTop_p4 = top_grand_d.p4();
			       }

			   }      //for loop, W's daus                                                                  
			 if (!isTheRightTop) continue; // move out from the checking loop, reject this top  
			 
		       } // if W from top 
		     
		     if (std::abs(d.pdgId()) == 5) // if b from top
		       {
			 b_fromTop_fromH_p4 = d.p4();
		       }
		   }//for top from Higgs daus
		 
		 if (!isTheRightTop) continue; // move out from the checking loop, reject this top
		 
		 double deltaR_b_fromTop_W    = ROOT::Math::VectorUtil::DeltaR(W_fromTop_fromH_p4,b_fromTop_fromH_p4);
		 h_WfromHiggstop_samedecayb_dR         -> Fill(deltaR_b_fromTop_W);
		 h_WfromHiggstop_samedecayb_dR_Vs_W_Pt -> Fill(deltaR_b_fromTop_W, W_fromTop_fromH_p4.pt());
		 h_WfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt -> Fill(deltaR_b_fromTop_W, Top_fromH_p4.pt());

		 if(deltaR_b_fromTop_W < 0.8) 
		   {
		     h_for_WbfromH_intofatjet_W_Pt     -> Fill(W_fromTop_fromH_p4.pt());
		     h_for_WbfromH_intofatjet_Top_Pt   -> Fill(Top_fromH_p4.pt());
		   }

		 double deltaR_b_fromTop_obj1 = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, obj1_fromW_fromTop_p4);
		 double deltaR_b_fromTop_obj2 = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, obj2_fromW_fromTop_p4);
		 double deltaR_obj1_obj2      = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, obj2_fromW_fromTop_p4);
		 double deltaRmax = -1;  // give an initial, non sense value
		 double deltaRmin = 1e6; //give an initial, non sense value

		 // find the dRmax
		 if(deltaR_b_fromTop_obj1 > deltaR_b_fromTop_obj2 && deltaR_b_fromTop_obj1 > deltaR_obj1_obj2)      deltaRmax = deltaR_b_fromTop_obj1;
		 else if(deltaR_b_fromTop_obj2 > deltaR_b_fromTop_obj1 && deltaR_b_fromTop_obj2 > deltaR_obj1_obj2) deltaRmax = deltaR_b_fromTop_obj2;
		 else                                                                                               deltaRmax = deltaR_obj1_obj2;
	 
		 h_objectsfromHiggstop_maxdR                  -> Fill(deltaRmax);
		 h_objectsfromHiggstop_maxdR_Vs_Higgstop_Pt   -> Fill(deltaRmax,p.pt());
		 h_objectsfromHiggstop_maxdR_Vs_Higgs_Pt      -> Fill(deltaRmax,m.pt());

		 //find the dRmin
		 if(deltaR_b_fromTop_obj1 < deltaR_b_fromTop_obj2 && deltaR_b_fromTop_obj1 < deltaR_obj1_obj2)      deltaRmin = deltaR_b_fromTop_obj1;
                 else if(deltaR_b_fromTop_obj2 < deltaR_b_fromTop_obj1 && deltaR_b_fromTop_obj2 < deltaR_obj1_obj2) deltaRmin = deltaR_b_fromTop_obj2;
                 else                                                                                               deltaRmin = deltaR_obj1_obj2;
		 
		 h_objectsfromHiggstop_mindR                                -> Fill(deltaRmin);
		 //
		 h_objectsfromHiggstop_Prob_mindR_lt_p8                     -> Fill("<0.8",0);
		 if(deltaRmin < 0.8) 
		   {
		     h_objectsfromHiggstop_Prob_mindR_lt_p8 -> Fill("<0.8",1);
		     h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("< 400",0);
		     h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("400<pT<500",0);
		     if     (Top_fromH_p4.pt() > 500.0) h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("> 500",1);
		     else if(Top_fromH_p4.pt() > 400.0) h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("400<pT<500",1);
		     else                               h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("< 400",1);
		   }
			    
		 else  h_objectsfromHiggstop_Prob_mindR_lt_p8               -> Fill(">0.8",1);
		 //
                 h_objectsfromHiggstop_mindR_Vs_Higgstop_Pt                 -> Fill(deltaRmin,p.pt());
                 h_objectsfromHiggstop_mindR_Vs_Higgs_Pt                    -> Fill(deltaRmin,m.pt());

		 // try with gen-jets ++++Matching++++++++++++++++++++++++++++++++++++++++++
		 matched_genbJet_p4 = p4_initializer;
                 for(size_t gj=0; gj < genJets_from_W_p4.size(); gj++)
                   {
                     genJets_from_W_p4.at(gj) = p4_initializer;
                   }
                 numof_obj_fromtop_matchedwith_genjet = 0;
		 bJetInsideFatJet = false;

		 double _deltaRMin = 1e6;
		 int matchedbjet = -1;
		 // match the bQuark
		 for (size_t bjet=0; bjet < selectedBJets_p4.size(); bjet++)
                   {
		     math::XYZTLorentzVector genbJet_p4 = selectedBJets_p4.at(bjet);
                     // find the b-Jet,  that has the min dR with the b from top
                     double deltaR = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, genbJet_p4);

                     if (deltaR < _deltaRMin)
                       {
                         _deltaRMin = deltaR;
			 matched_genbJet_p4 = genbJet_p4;
			 matchedbjet = bjet;
                       }
                   } //For-loop: All selected b-jets 
		 if (_deltaRMin < 0.4)
		   {
		     numof_obj_fromtop_matchedwith_genjet++ ;
		   }
		 h_Hs_deltaR_bQuark_closestbJet   -> Fill(_deltaRMin);
		 

		 //match the bquark from Higgs
		 _deltaRMin = 1e6;
		 for (size_t bjet=0; bjet < selectedBJets_p4.size(); bjet++)
                   {
		     if (int(bjet) == matchedbjet) continue; //make sure that you will not match the same b-jet to both b-quarks  
		     math::XYZTLorentzVector genbJet_p4 = selectedBJets_p4.at(bjet);
                     // find the b-Jet,  that has the min dR with the b from top       
                     double deltaR = ROOT::Math::VectorUtil::DeltaR(b_fromH_p4, genbJet_p4);

                     if (deltaR < _deltaRMin)
                       {
                         _deltaRMin = deltaR;
                         matched_genbJet_fromH_p4 = genbJet_p4;
		       }
                   } //For-loop: All selected b-jets
		 h_Hs_deltaR_bQuark_fromH_closestbJet -> Fill(_deltaRMin);
		 if(_deltaRMin < 0.4) bfromHiggsMatched = true;
		 
		 // match the quarks from W
		 _deltaRMin = 1e6;
		 int matchedjet = -1; //the jet that will be matched, give an initial non sense value  

		 // For-loop: All selected jets, No b-jets (i)                                                                                                                                              
		 for (size_t nobjet=0; nobjet < selJets_NoBJets_p4.size(); nobjet++)
		   {
		     math::XYZTLorentzVector gen_nobJet_p4 = selJets_NoBJets_p4.at(nobjet);

		     // find the Genjet, no b Jets,  that has the min dR with the q from W                                                                                                                     
		     double deltaR = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, gen_nobJet_p4);

		     if (deltaR < _deltaRMin)
		       {
			 _deltaRMin = deltaR;
			 genJets_from_W_p4.at(0) = gen_nobJet_p4;
			 matchedjet = int(nobjet);
		       }
		   } //For-loop: All selected jets, No b-jets (i)
		 if (_deltaRMin < 0.4)
                   {
                     numof_obj_fromtop_matchedwith_genjet++ ;
                   }
		 h_Hs_deltaR_QfromW_closestJet   -> Fill(_deltaRMin);

		 _deltaRMin = 1e6;
		 // For-loop: All selected jets, No b-jets (ii)		 
		 for (size_t nobjet=0; nobjet < selJets_NoBJets_p4.size(); nobjet++)
                   {
                     if (int(nobjet) == matchedjet) continue; //make sure that you will not match the same jet to both quarks from W                         
		     math::XYZTLorentzVector gen_nobJet_p4 = selJets_NoBJets_p4.at(nobjet);

                     // find the Genjet, no b Jets,  that has the min dR with the q from W                                                                                                   
                     double deltaR = ROOT::Math::VectorUtil::DeltaR(obj2_fromW_fromTop_p4, gen_nobJet_p4);

                     if (deltaR < _deltaRMin)
                       {
                         _deltaRMin = deltaR;
                         genJets_from_W_p4.at(1) = gen_nobJet_p4;
                       }
                   } //For-loop: All selected jets, No b-jets (ii)
		 if (_deltaRMin < 0.4)
                   {
                     numof_obj_fromtop_matchedwith_genjet++ ;
                   }
		 h_Hs_deltaR_QfromW_closestJet   -> Fill(_deltaRMin);

		 h_Hs_numof_obj_fromtop_matchedwith_genjet -> Fill(numof_obj_fromtop_matchedwith_genjet);

		 //Gen-Reconstruction of top, W and Higgs from jets
		 if(numof_obj_fromtop_matchedwith_genjet == 3)
		   {
		     W_genReco_fromTop_fromH_p4             = genJets_from_W_p4.at(0)    + genJets_from_W_p4.at(1);
		     Top_genReco_fromH_p4                   = W_genReco_fromTop_fromH_p4 + matched_genbJet_p4;
		     if(bfromHiggsMatched) Higgs_genReco_p4 = Top_genReco_fromH_p4       + matched_genbJet_fromH_p4;
		   }
		   
		 // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 /////////////////////////////////////////////////BaryCenter////////////////////////////////////////////////////////////////
		 double deltaR_b_fromTop_Top    = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, Top_fromH_p4);
                 double deltaR_obj1_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, Top_fromH_p4);
                 double deltaR_obj2_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(obj2_fromW_fromTop_p4, Top_fromH_p4);
		 unsigned int decayproductsintofatjet = 0;
		 double objtop_deltaRmax = -1;  // give an initial, non sense value
		 double objtop_deltaRmin = 999;  // give an initial, non sense value
		 
		 h_Hs_which_objectfromtop_maxdR -> Fill ("W'obj",0) ; //just to determinate the first label

		 if(deltaR_b_fromTop_Top > deltaR_obj1_fromTop_Top && deltaR_b_fromTop_Top > deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmax = deltaR_b_fromTop_Top;
		     h_Hs_which_objectfromtop_maxdR -> Fill("b",1) ;

		     if(deltaR_b_fromTop_Top > 0.8)
		       {
			 h_Hs_mostdistantfromtop_isb__top_pT      -> Fill(Top_fromH_p4.pt());
			 h_Hs_mostdistantfromtop_isb__W_pT        -> Fill(W_fromTop_fromH_p4.pt());
			 h_Hs_mostdistantfromtop_isb__b_pT        -> Fill(b_fromTop_fromH_p4.pt());
			 h_Hs_mostdistantfromtop_isb__dRmax_b_top -> Fill(deltaR_b_fromTop_Top);
			 			 
			 if(deltaR_b_fromTop_obj1 < deltaR_b_fromTop_obj2) h_Hs_mostdistantfromtop_isb_dRmin_b_objfromW -> Fill(deltaR_b_fromTop_obj1);
			 else                                              h_Hs_mostdistantfromtop_isb_dRmin_b_objfromW -> Fill(deltaR_b_fromTop_obj2);
		       }
		     
		     if(deltaR_b_fromTop_Top > 0.8 && deltaR_obj1_obj2 < 0.7)
		       {
				 h_Hs_mostdistantfromtop_isb__dRqq_Vs_W_pT   -> Fill(deltaR_obj1_obj2,W_fromTop_fromH_p4.pt());
			 h_Hs_mostdistantfromtop_isb__dRqq_Vs_top_pT -> Fill(deltaR_obj1_obj2,Top_fromH_p4.pt());
			 
			 h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("< 400",0);
			 h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("400<pT<500",0);
			 h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("< 200",0);
                         h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("200<pT<300",0);
			 if     (Top_fromH_p4.pt() > 500.0) h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("> 500",1);
			 else if(Top_fromH_p4.pt() > 400.0) h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("400<pT<500",1);
			 else                               h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("< 400",1);
			 
			 if     (W_fromTop_fromH_p4.pt() > 300.0) h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("> 300",1);
                         else if(W_fromTop_fromH_p4.pt() > 200.0) h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("200<pT<300",1);
                         else                               h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("< 200",1);
		       }
		     
		     h_Hs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("True",0);
		     if(deltaR_b_fromTop_Top > 0.8 && deltaR_obj1_obj2 < 0.7) h_Hs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("True",1);
		     else                                                     h_Hs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("False",1);
                     		     		     
                   }
                 else if(deltaR_obj1_fromTop_Top > deltaR_b_fromTop_Top && deltaR_obj1_fromTop_Top > deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmax = deltaR_obj1_fromTop_Top;
		     h_Hs_which_objectfromtop_maxdR-> Fill("W'obj",1);
                   }
                 else
		   {
		     objtop_deltaRmax = deltaR_obj2_fromTop_Top;
		     h_Hs_which_objectfromtop_maxdR-> Fill("W'obj",1);
                   }
		 h_Hs_objectsfromtop_top_maxdR -> Fill(objtop_deltaRmax);
		 
		 if(deltaR_b_fromTop_Top < deltaR_obj1_fromTop_Top && deltaR_b_fromTop_Top < deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmin = deltaR_b_fromTop_Top;
                   }
                 else if(deltaR_obj1_fromTop_Top < deltaR_b_fromTop_Top && deltaR_obj1_fromTop_Top < deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmin = deltaR_obj1_fromTop_Top;
                   }
                 else objtop_deltaRmin = deltaR_obj2_fromTop_Top;
		 h_Hs_objectsfromtop_top_mindR -> Fill(objtop_deltaRmin);   



		 if(deltaR_b_fromTop_Top < 0.8)    // b into fatjet 
		   {
		     decayproductsintofatjet++;
		     bQuarkIntoFatJet = true;
		   }
		 else bQuarkIntoFatJet = false;
		 if(deltaR_obj1_fromTop_Top < 0.8) decayproductsintofatjet++;// obj1 into fatjet
		 if(deltaR_obj2_fromTop_Top < 0.8) decayproductsintofatjet++;// obj2 into fatjet
	     
		 h_Probdecayproductsintofatjet_Hs         -> Fill("0",0);  //just to determinate the label of the first bin
		 h_Probdecayproductsintofatjet_Hs         -> Fill("1",0);  //just to determinate the label of the second bin      
		 h_Probdecayproductsintofatjet_Hs         -> Fill("2",0);  //just to determinate the label of the third bin 

		 h_Hs_QuarksintofatjetMultiplicity -> Fill("0",0);  //just to determinate the label of the first bin                          
		 h_Hs_QuarksintofatjetMultiplicity -> Fill("j",0);
		 h_Hs_QuarksintofatjetMultiplicity -> Fill("b",0);
		 h_Hs_QuarksintofatjetMultiplicity -> Fill("jj",0);
		 h_Hs_QuarksintofatjetMultiplicity -> Fill("bj",0);

		 h_Hs_isbQuarkintofatjet           -> Fill("No b",0);  //just to determinate the label of the first bin

                 if      (decayproductsintofatjet == 0)   
		   {
		     h_Probdecayproductsintofatjet_Hs  -> Fill("0",1);
		     h_Hs_QuarksintofatjetMultiplicity -> Fill("0",1);
		     h_Hs_isbQuarkintofatjet           -> Fill("No b",1);
		   }
		 else if (decayproductsintofatjet == 1) 
		   {
		     h_Probdecayproductsintofatjet_Hs -> Fill("1",1); 
		     if(bQuarkIntoFatJet)
		       {
			 h_Hs_QuarksintofatjetMultiplicity -> Fill("b",1);
			 h_Hs_isbQuarkintofatjet           -> Fill("b",1);
		       }
		     else
		       {
			 h_Hs_QuarksintofatjetMultiplicity -> Fill("j",1);
                         h_Hs_isbQuarkintofatjet           -> Fill("No b",1); 
		       }
		   }
                 else if (decayproductsintofatjet == 2)  
		   {
		     h_Probdecayproductsintofatjet_Hs -> Fill("2",1);
		     if(bQuarkIntoFatJet)
                       {
                         h_Hs_QuarksintofatjetMultiplicity -> Fill("bj",1);
                         h_Hs_isbQuarkintofatjet           -> Fill("b",1);
                       }
                     else
                       {
                         h_Hs_QuarksintofatjetMultiplicity -> Fill("jj",1);
                         h_Hs_isbQuarkintofatjet           -> Fill("No b",1);
                       }
		   }
		 else                                    
		   {
		     h_Probdecayproductsintofatjet_Hs  -> Fill("All 3",1);
		     h_Hs_QuarksintofatjetMultiplicity -> Fill("bjj",1);
		     h_Hs_isbQuarkintofatjet           -> Fill("b",1);
		   }
		 
		 if (decayproductsintofatjet == 3)
		   {
		     unsigned int Num_pairs_withenough_deltaR_amongthemselves = 0;
		     if      (deltaR_b_fromTop_obj1 < 0.4 && deltaR_b_fromTop_obj2 < 0.4 && deltaR_obj1_obj2 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 0;

		     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 < 0.4 && deltaR_b_fromTop_obj2 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 1;
		     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 < 0.4 && deltaR_b_fromTop_obj1 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 1;
		     else if (deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj1 < 0.4 && deltaR_b_fromTop_obj2 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 1;
		     
		     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj2 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 2;
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj1 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 2;
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 <0.4 ) Num_pairs_withenough_deltaR_amongthemselves = 2;
		     
		     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 > 0.4) Num_pairs_withenough_deltaR_amongthemselves = 3;		     
		     
		     
		     h_Pairsinbarycenter_enoughdeltaR_Hs         -> Fill("0",0);  //just to determinate the label of the first bin     
		     h_Pairsinbarycenter_enoughdeltaR_Hs         -> Fill("1 pair",0);  //just to determinate the label of the second bin      
		     h_Pairsinbarycenter_enoughdeltaR_Hs         -> Fill("2 pairs",0);  //just to determinate the label of the third bin  

		     if      (Num_pairs_withenough_deltaR_amongthemselves == 0)   h_Pairsinbarycenter_enoughdeltaR_Hs -> Fill("0",1);
		     else if (Num_pairs_withenough_deltaR_amongthemselves == 1)   h_Pairsinbarycenter_enoughdeltaR_Hs -> Fill("1 pair",1);
		     else if (Num_pairs_withenough_deltaR_amongthemselves == 2)   h_Pairsinbarycenter_enoughdeltaR_Hs -> Fill("2 pairs",1);
		     else                                                         h_Pairsinbarycenter_enoughdeltaR_Hs -> Fill("All 3",1);
		     
		     h_Iffatjet_Hs_Top_pT         -> Fill (Top_fromH_p4.pt());
		     h_Iffatjet_Hs_W_pT           -> Fill (W_fromTop_fromH_p4.pt());
		     h_Iffatjet_Hs_Top_pT_Vs_W_pT -> Fill (Top_fromH_p4.pt(),W_fromTop_fromH_p4.pt());
		     h_Iffatjet_Hs_bfromH_pT      -> Fill (b_fromH_p4.pt());
		     
		     h_Iffatjet_Hs_EventsWithHighTop_pT -> Fill("top.pT < 400 GeV", 0);
		     h_Iffatjet_Hs_EventsWithHighW_pT   -> Fill("W.pT < 200 GeV", 0);
		     if(Top_fromH_p4.pt() > 400.0)       h_Iffatjet_Hs_EventsWithHighTop_pT -> Fill("top.pT > 400 GeV", 1);
		     else                                h_Iffatjet_Hs_EventsWithHighTop_pT -> Fill("top.pT < 400 GeV", 1);
		     if(W_fromTop_fromH_p4.pt() > 200.0) h_Iffatjet_Hs_EventsWithHighW_pT   -> Fill("W.pT > 200 GeV", 1);
                     else                                h_Iffatjet_Hs_EventsWithHighW_pT   -> Fill("W.pT < 200 GeV", 1);
		     if(Top_fromH_p4.pt() > 400.0 && W_fromTop_fromH_p4.pt() > 200.0) 
		       {
			 h_Iffatjet_Hs_HighWandTop_pT_dRmin -> Fill (deltaRmin);
		       }
		     
		     double dR_bfromHiggstop_Higgstop   = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, Top_fromH_p4);
		     h_Iffatjet_Hs_bfromT_Top_dR        -> Fill(dR_bfromHiggstop_Higgstop);
		     double dEta_bfromHiggstop_Higgstop = std::abs(b_fromTop_fromH_p4.eta() - Top_fromH_p4.eta());
		     h_Iffatjet_Hs_bfromT_Top_dEta      -> Fill(dEta_bfromHiggstop_Higgstop);
		     double dPhi_bfromHiggstop_Higgstop = std::abs(b_fromTop_fromH_p4.phi() - Top_fromH_p4.phi());
		     h_Iffatjet_Hs_bfromT_Top_dPhi      -> Fill(dPhi_bfromHiggstop_Higgstop);
		     
		     double dR_bfromHiggstop_W    = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, W_fromTop_fromH_p4);
                     h_Iffatjet_Hs_bfromT_W_dR    -> Fill(dR_bfromHiggstop_W);
                     double dEta_bfromHiggstop_W  = std::abs(b_fromTop_fromH_p4.eta() - W_fromTop_fromH_p4.eta());
                     h_Iffatjet_Hs_bfromT_W_dEta  -> Fill(dEta_bfromHiggstop_W);
                     double dPhi_bfromHiggstop_W  = std::abs(b_fromTop_fromH_p4.phi() - W_fromTop_fromH_p4.phi());
                     h_Iffatjet_Hs_bfromT_W_dPhi  -> Fill(dPhi_bfromHiggstop_W);
		     
		     double dR_bfromHiggstop_bfromH    = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, b_fromH_p4);
                     h_Iffatjet_Hs_bfromT_bfromH_dR    -> Fill(dR_bfromHiggstop_bfromH);
                     double dEta_bfromHiggstop_bfromH  = std::abs(b_fromTop_fromH_p4.eta() - b_fromH_p4.eta());
                     h_Iffatjet_Hs_bfromT_bfromH_dEta  -> Fill(dEta_bfromHiggstop_bfromH);
                     double dPhi_bfromHiggstop_bfromH  = std::abs(b_fromTop_fromH_p4.phi() - b_fromH_p4.phi());
                     h_Iffatjet_Hs_bfromT_bfromH_dPhi  -> Fill(dPhi_bfromHiggstop_bfromH);
		     
		     
		   } // if all 3 decay products are into 0.8

		 // try with gen-jets ++++Barycenter+++++++++++++++++++++++++++++++++++++++++
		 if(numof_obj_fromtop_matchedwith_genjet == 3)
		   {
		     double deltaR_bjet_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(matched_genbJet_p4, Top_fromH_p4);
		     double deltaR_jet1_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(genJets_from_W_p4.at(0), Top_fromH_p4);
		     double deltaR_jet2_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(genJets_from_W_p4.at(1), Top_fromH_p4);
		     unsigned int jets_intofatjet = 0;

		     if(deltaR_bjet_fromTop_Top < 0.8)
		       {
			 jets_intofatjet++;// bjet into fatjet        
			 bJetInsideFatJet = true;
		       }
		     else bJetInsideFatJet = false;
		     if(deltaR_jet1_fromTop_Top < 0.8) jets_intofatjet++;// jet1 into fatjet
		     if(deltaR_jet2_fromTop_Top < 0.8) jets_intofatjet++;// jet2 into fatjet         

		     // labels
		     h_Hs_ProbdecayJetsintofatjet   -> Fill("0",0);  //just to determinate the label of the first bin                             
		     h_Hs_ProbdecayJetsintofatjet   -> Fill("1",0);  //just to determinate the label of the second bin                                                           
		     h_Hs_ProbdecayJetsintofatjet   -> Fill("2",0);  //just to determinate the label of the third bin                                                           
		     
		     h_Hs_JetsintofatjetMultiplicity-> Fill("0",0);  //just to determinate the label of the first bin
		     h_Hs_JetsintofatjetMultiplicity-> Fill("j",0);
		     h_Hs_JetsintofatjetMultiplicity-> Fill("b",0);
		     h_Hs_JetsintofatjetMultiplicity-> Fill("jj",0);
                     h_Hs_JetsintofatjetMultiplicity-> Fill("bj",0);

		     h_Hs_isbJetintofatjet          -> Fill("No b-jet",0);  //just to determinate the label of the first bin 
		     //

		     if      (jets_intofatjet == 0)
		       {
			 h_Hs_ProbdecayJetsintofatjet     -> Fill("0",1);
			 h_Hs_JetsintofatjetMultiplicity  -> Fill("0",1);
			 h_Hs_isbJetintofatjet            -> Fill("No b-jet",1);
		       }

		     else if (jets_intofatjet == 1) 
		       {
			 h_Hs_ProbdecayJetsintofatjet -> Fill("1",1);
			 if(bJetInsideFatJet)
			   {
			     h_Hs_JetsintofatjetMultiplicity-> Fill("b",1);
			     h_Hs_isbJetintofatjet          -> Fill("b-jet",1);
			   }
			 else
			   {
			     h_Hs_JetsintofatjetMultiplicity-> Fill("j",1);
			     h_Hs_isbJetintofatjet          -> Fill("No b-jet",1);
			   }
		       }

		     else if (jets_intofatjet == 2)
		       {
			 h_Hs_ProbdecayJetsintofatjet -> Fill("2",1);
			 if(bJetInsideFatJet)
                           {
                             h_Hs_JetsintofatjetMultiplicity-> Fill("bj",1);
                             h_Hs_isbJetintofatjet          -> Fill("b-jet",1);
                           }
                         else
                           {
                             h_Hs_JetsintofatjetMultiplicity-> Fill("jj",1);
                             h_Hs_isbJetintofatjet          -> Fill("No b-jet",1);
                           }
                       }

		     else
		       {
			 h_Hs_ProbdecayJetsintofatjet -> Fill("All 3",1);
			 h_Hs_JetsintofatjetMultiplicity-> Fill("bjj",1);
			 h_Hs_isbJetintofatjet          -> Fill("b-jet",1);
		       }
		     
		     h_Hs_BjetInsideFatJet_Top_Pt -> Fill(bJetInsideFatJet, Top_genReco_fromH_p4.pt());
		     h_Hs_BjetInsideFatJet_W_Pt   -> Fill(bJetInsideFatJet, W_genReco_fromTop_fromH_p4.pt());
                     if(bfromHiggsMatched) h_Hs_BjetInsideFatJet_H_Pt   -> Fill(bJetInsideFatJet, Higgs_genReco_p4.pt());
		   }
		 
		 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 
		 
		 double deltaR_bfromHiggstop_Higgstop = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, p.p4()); //find dR(bfromtfromH,tfromH)       
		 h_bfromHiggstop_Higgstop_dR                -> Fill(deltaR_bfromHiggstop_Higgstop);
		 
		 h_b_top_fromhiggs_underdR                  -> Fill("< 0.8",0);  //just to determinate the label of the first bin             
		 if(deltaR_bfromHiggstop_Higgstop < 0.8)  h_b_top_fromhiggs_underdR   -> Fill("< 0.8",1);
		 else                                     h_b_top_fromhiggs_underdR   -> Fill("> 0.8",1);

		 h_bfromHiggstop_Higgstop_dR_Vs_Higgstop_Pt -> Fill(deltaR_bfromHiggstop_Higgstop,p.pt());
		 ///////////////////////////////////////////////////////////////////////////////
		 double deltaR_WfromHiggs_closest_b = 9999.9; //non sense initial value
		 for (auto& gp: fEvent.genparticles().getGenParticles()) // for genpar %2
		   {
		     
		     if(std::abs(gp.pdgId()) == 5 && gp.isFirstCopy()) // find the b
		       {
			 double deltaR_WfromHiggs_any_b = ROOT::Math::VectorUtil::DeltaR(W_fromTop_fromH_p4,gp.p4());
			 if ( deltaR_WfromHiggs_any_b < deltaR_WfromHiggs_closest_b) deltaR_WfromHiggs_closest_b = deltaR_WfromHiggs_any_b;
			 else continue;			   
		       }
		   }// for genpar %2
		 
		 h_WfromHiggstop_closestb_dR                                 -> Fill(deltaR_WfromHiggs_closest_b);
		 h_WfromHiggstop_closestb_dR_Vs_W_Pt                         -> Fill(deltaR_WfromHiggs_closest_b,W_fromTop_fromH_p4.pt());
		 h_WfromHiggstop_closestb_dR_Vs_samedecaytop_Pt  -> Fill(deltaR_WfromHiggs_closest_b,Top_fromH_p4.pt());
		 //////////////////////////////////////////////////////////////////////////////////
		 
	       }   // if top from Higgs

	     else //if top is not from Higgs//////////////////////////////////////////////////////////////////////////////////////////////
	       {
		 math::XYZTLorentzVector b_fromTop_NOTfromH_p4(0,0,0,0), obj1_fromW_fromTop_p4(0,0,0,0), obj2_fromW_fromTop_p4(0,0,0,0);
		 math::XYZTLorentzVector W_fromTop_NOTfromH_p4(0,0,0,0), Top_NOTfromH_p4(0,0,0,0);
		 math::XYZTLorentzVector W_genReco_fromTop_NOTfromH_p4(0,0,0,0), Top_genReco_NOTfromH_p4(0,0,0,0);
		 bool isTheRightTop = true;
		 bool obj1_exist = false;
		 bool bQuarkIntoFatJet = false;

		 Top_NOTfromH_p4 = p.p4();
		 
		 unsigned int lastcopy = GetTheLastCopy(p.index()); // we need to find the last copy in order to find the real daughters      
                 genParticle lastT; //create the particle object of the last copy of top                                                      
                 lastT =  fEvent.genparticles().getGenParticles()[lastcopy];
		 std::vector<short> lastT_daughters = lastT.daughters();
		 
                 for (unsigned int i = 0; i < lastT_daughters.size() ; i++)
                   {
                     genParticle d; //create the particle object 
                     d =  fEvent.genparticles().getGenParticles()[lastT_daughters.at(i)];
		     if (std::abs(d.pdgId()) == 24) // if W from top                                                                         
		       {
			 W_fromTop_NOTfromH_p4 = d.p4();
                         genParticle lastW;
                         unsigned int lastWcopy = GetTheLastCopy(d.index()); // we need to find the last copy in order to find the real daughters                                                                                                                                           
			 lastW = fEvent.genparticles().getGenParticles()[lastWcopy];
			 std::vector<short> W_daughters = lastW.daughters();
			 for (size_t k = 0; k < W_daughters.size() ; k++)
			   {
			     genParticle top_grand_d; //create the particle object, i define this as grand_d, the daus of the W, so the grand-daughters of Top                                                                                                                              
			     top_grand_d =  fEvent.genparticles().getGenParticles()[W_daughters.at(k)];
			     if ( !mcTools.IsQuark(top_grand_d.pdgId()) ) //if the daus of W are not quarks (f.e. leptons), we reject the top 
			       {
				 isTheRightTop = false;
				 continue;
			       }
			     
			     if (obj1_exist == false)
			       {
				 obj1_exist = true;
				 obj1_fromW_fromTop_p4 = top_grand_d.p4();
                               }
                             else
			       {
                                 obj2_fromW_fromTop_p4 = top_grand_d.p4();
                               }
                           }      //for loop, W's daus                                                                   
                         if (!isTheRightTop) continue; // move out from the checking loop, reject this top          
                       } // if W from top                                  
		     
                     if (std::abs(d.pdgId()) == 5) // if b from top                                                                           
                       {
                         b_fromTop_NOTfromH_p4 = d.p4();
                       }
                   }//for top NOT from Higgs daus                                                             
		 
                 if (!isTheRightTop) continue; // move out from the checking loop, reject this top 
		 
		 double deltaR_b_fromTop_W    = ROOT::Math::VectorUtil::DeltaR(W_fromTop_NOTfromH_p4,b_fromTop_NOTfromH_p4);
                 h_WNOTfromHiggstop_samedecayb_dR         -> Fill(deltaR_b_fromTop_W);
                 h_WNOTfromHiggstop_samedecayb_dR_Vs_W_Pt -> Fill(deltaR_b_fromTop_W, W_fromTop_NOTfromH_p4.pt());
                 h_WNOTfromHiggstop_samedecayb_dR_Vs_samedecaytop_Pt -> Fill(deltaR_b_fromTop_W, Top_NOTfromH_p4.pt());

		 if(deltaR_b_fromTop_W < 0.8)
                   {
                     h_for_WbNOTfromH_intofatjet_W_Pt     -> Fill(W_fromTop_NOTfromH_p4.pt());
                     h_for_WbNOTfromH_intofatjet_Top_Pt   -> Fill(Top_NOTfromH_p4.pt());
                   }

		 double deltaR_b_fromTop_obj1 = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, obj1_fromW_fromTop_p4);
                 double deltaR_b_fromTop_obj2 = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, obj2_fromW_fromTop_p4);
                 double deltaR_obj1_obj2      = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, obj2_fromW_fromTop_p4);
                 double deltaRmax = -1; // give an initial, non sense value                                        
		 double deltaRmin = 1e6; // give an initial, non sense value

		 // find the dRmin
                 if(deltaR_b_fromTop_obj1 > deltaR_b_fromTop_obj2 && deltaR_b_fromTop_obj1 > deltaR_obj1_obj2)
                   {
                     deltaRmax = deltaR_b_fromTop_obj1;
                   }
                 else if(deltaR_b_fromTop_obj2 > deltaR_b_fromTop_obj1 && deltaR_b_fromTop_obj2 > deltaR_obj1_obj2)
                   {
                     deltaRmax = deltaR_b_fromTop_obj2;
                   }
                 else deltaRmax = deltaR_obj1_obj2;

                 h_objectsNOTfromHiggstop_maxdR             -> Fill(deltaRmax);
                 h_objectsNOTfromHiggstop_maxdR_Vs_top_Pt   -> Fill(deltaRmax,p.pt());
		 
		 //find the dRmin among the top decay products 
                 if(deltaR_b_fromTop_obj1 < deltaR_b_fromTop_obj2 && deltaR_b_fromTop_obj1 < deltaR_obj1_obj2)      deltaRmin = deltaR_b_fromTop_obj1;
                 else if(deltaR_b_fromTop_obj2 < deltaR_b_fromTop_obj1 && deltaR_b_fromTop_obj2 < deltaR_obj1_obj2) deltaRmin = deltaR_b_fromTop_obj2;
                 else                                                                                               deltaRmin = deltaR_obj1_obj2;

                 h_objectsNOTfromHiggstop_mindR                                -> Fill(deltaRmin);

                 h_objectsNOTfromHiggstop_Prob_mindR_lt_p8                     -> Fill("<0.8",0);
                 if(deltaRmin < 0.8) 
		   {
		     h_objectsNOTfromHiggstop_Prob_mindR_lt_p8 -> Fill("<0.8",1);
		     h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("< 400",0);
                     h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("400<pT<500",0);
		     if     (Top_NOTfromH_p4.pt() > 500.0) h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("> 500",1);
                     else if(Top_NOTfromH_p4.pt() > 400.0) h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("400<pT<500",1);
                     else                                  h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("< 400",1);
		   }
                 else  h_objectsNOTfromHiggstop_Prob_mindR_lt_p8               -> Fill(">0.8",1);
                 h_objectsNOTfromHiggstop_mindR_Vs_top_Pt                      -> Fill(deltaRmin,p.pt());
                 
		 // try with gen-jets+++Matching++++++++++++++++++++++++++++++++++++++++++++
                 matched_genbJet_p4 = p4_initializer;
                 for(size_t gj=0; gj < genJets_from_W_p4.size(); gj++)
                   {
                     genJets_from_W_p4.at(gj) = p4_initializer;
                   }
                 numof_obj_fromtop_matchedwith_genjet = 0;
		 bJetInsideFatJet = false;

		 double _deltaRMin = 1e6;
                 // match the bQuark         
                 for (size_t bjet=0; bjet < selectedBJets_p4.size(); bjet++)
                   {
		     math::XYZTLorentzVector genbJet_p4 = selectedBJets_p4.at(bjet);
                     // find the b-Jet,  that has the min dR with the b from top        
                     double deltaR = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, genbJet_p4);

                     if (deltaR < _deltaRMin)
		       {
                         _deltaRMin = deltaR;
                         matched_genbJet_p4 = genbJet_p4;
                       }
                   } //For-loop: All selected b-jets         
                 if (_deltaRMin < 0.4)
                   {
                     numof_obj_fromtop_matchedwith_genjet++ ;
                   }
                 h_NotHs_deltaR_bQuark_closestbJet   -> Fill(_deltaRMin);

                 // match the quarks from W          
                 _deltaRMin = 1e6;
		 int matchedjet = -1; //the jet that will be matched, give an initial non sense value                                                                                                             

                 // For-loop: All selected jets, No b-jets (i)        
                 for (size_t nobjet=0; nobjet < selJets_NoBJets_p4.size(); nobjet++)
                   {
		     math::XYZTLorentzVector gen_nobJet_p4 = selJets_NoBJets_p4.at(nobjet);

                     // find the Genjet, no b Jets,  that has the min dR with the q from W         
                     double deltaR = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, gen_nobJet_p4);

                     if (deltaR < _deltaRMin)
                       {
                         _deltaRMin = deltaR;
                         genJets_from_W_p4.at(0) = gen_nobJet_p4;
                         matchedjet = int(nobjet);
                       }
                   } //For-loop: All selected jets, No b-jets (i)           

                 if (_deltaRMin < 0.4)
                   {
                     numof_obj_fromtop_matchedwith_genjet++ ;
                   }
                 h_NotHs_deltaR_QfromW_closestJet   -> Fill(_deltaRMin);
		 
                 _deltaRMin = 1e6;
                 // For-loop: All selected jets, No b-jets (ii) 
		 for (size_t nobjet=0; nobjet < selJets_NoBJets_p4.size(); nobjet++)
                   {
                     if (int(nobjet) == matchedjet) continue; //make sure that you will not match the same jet to both quarks from W          
		     math::XYZTLorentzVector gen_nobJet_p4 = selJets_NoBJets_p4.at(nobjet);
		     
                     // find the Genjet, no b Jets,  that has the min dR with the q from W      
                     double deltaR = ROOT::Math::VectorUtil::DeltaR(obj2_fromW_fromTop_p4, gen_nobJet_p4);
		     
                     if (deltaR < _deltaRMin)
                       {
                         _deltaRMin = deltaR;
                         genJets_from_W_p4.at(1) = gen_nobJet_p4;
                       }
                   } //For-loop: All selected jets, No b-jets (ii)    

                 if (_deltaRMin < 0.4)
                   {
                     numof_obj_fromtop_matchedwith_genjet++ ;
                   }
                 h_NotHs_deltaR_QfromW_closestJet   -> Fill(_deltaRMin);

                 h_NotHs_numof_obj_fromtop_matchedwith_genjet -> Fill(numof_obj_fromtop_matchedwith_genjet);
		 
		 if(numof_obj_fromtop_matchedwith_genjet == 3)
                   {
                     W_genReco_fromTop_NOTfromH_p4   = genJets_from_W_p4.at(0)       + genJets_from_W_p4.at(1);
                     Top_genReco_NOTfromH_p4         = W_genReco_fromTop_NOTfromH_p4 + matched_genbJet_p4;
		   }
                 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                 /////////////////////////////////////////////////BaryCenter//////////////////////////////////////////////////////////////////
		 /////////////assuming that is the top-direction /////////////////////////////////////////////////////////////////////////////
                 double deltaR_b_fromTop_Top    = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, Top_NOTfromH_p4);
                 double deltaR_obj1_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, Top_NOTfromH_p4);
                 double deltaR_obj2_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(obj2_fromW_fromTop_p4, Top_NOTfromH_p4);
                 unsigned int decayproductsintofatjet = 0;
		 double objtop_deltaRmax = -1;  // give an initial, non sense value                 
                 double objtop_deltaRmin = 999;  // give an initial, non sense value                 

                 h_NotHs_which_objectfromtop_maxdR -> Fill ("W'obj",0) ; //just to determinate the first label              

                 if(deltaR_b_fromTop_Top > deltaR_obj1_fromTop_Top && deltaR_b_fromTop_Top > deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmax = deltaR_b_fromTop_Top;
                     h_NotHs_which_objectfromtop_maxdR -> Fill("b",1) ;
		     
                     if(deltaR_b_fromTop_Top > 0.8)
                       {
                         h_NotHs_mostdistantfromtop_isb__top_pT      -> Fill(Top_NOTfromH_p4.pt());
                         h_NotHs_mostdistantfromtop_isb__W_pT        -> Fill(W_fromTop_NOTfromH_p4.pt());
                         h_NotHs_mostdistantfromtop_isb__b_pT        -> Fill(b_fromTop_NOTfromH_p4.pt());
			 h_NotHs_mostdistantfromtop_isb__dRmax_b_top -> Fill(deltaR_b_fromTop_Top);
			 
                         if(deltaR_b_fromTop_obj1 < deltaR_b_fromTop_obj2) h_NotHs_mostdistantfromtop_isb_dRmin_b_objfromW -> Fill(deltaR_b_fromTop_obj1);
                         else                                              h_NotHs_mostdistantfromtop_isb_dRmin_b_objfromW -> Fill(deltaR_b_fromTop_obj2);
                       }
		     
		     if(deltaR_b_fromTop_Top > 0.8 && deltaR_obj1_obj2 < 0.7)
		       {
                         h_NotHs_mostdistantfromtop_isb__dRqq_Vs_W_pT   -> Fill(deltaR_obj1_obj2,W_fromTop_NOTfromH_p4.pt());
                         h_NotHs_mostdistantfromtop_isb__dRqq_Vs_top_pT -> Fill(deltaR_obj1_obj2,Top_NOTfromH_p4.pt());
			 
			 h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("< 400",0);
                         h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("400<pT<500",0);
			 h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("< 200",0);
                         h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("200<pT<300",0);
                         if     (Top_NOTfromH_p4.pt() > 500.0) h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("> 500",1);
                         else if(Top_NOTfromH_p4.pt() > 400.0) h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("400<pT<500",1);
                         else                                  h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("< 400",1);

                         if     (W_fromTop_NOTfromH_p4.pt() > 300.0) h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill(">300",1);
                         else if(W_fromTop_NOTfromH_p4.pt() > 200.0) h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("200<pT<300",1);
                         else                                        h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("< 200",1);
                       }

                     h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("True",0);
                     if(deltaR_b_fromTop_Top > 0.8 && deltaR_obj1_obj2 < 0.7) h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("True",1);
                     else                                                     h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("False",1);
                   }
                 else if(deltaR_obj1_fromTop_Top > deltaR_b_fromTop_Top && deltaR_obj1_fromTop_Top > deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmax = deltaR_obj1_fromTop_Top;
                     h_NotHs_which_objectfromtop_maxdR-> Fill("W'obj",1);
                   }
                 else
                   {
                     objtop_deltaRmax = deltaR_obj2_fromTop_Top;
                     h_NotHs_which_objectfromtop_maxdR-> Fill("W'obj",1);
                   }
                 h_NotHs_objectsfromtop_top_maxdR -> Fill(objtop_deltaRmax);

                 if(deltaR_b_fromTop_Top < deltaR_obj1_fromTop_Top && deltaR_b_fromTop_Top < deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmin = deltaR_b_fromTop_Top;
                   }
                 else if(deltaR_obj1_fromTop_Top < deltaR_b_fromTop_Top && deltaR_obj1_fromTop_Top < deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmin = deltaR_obj1_fromTop_Top;
                   }
                 else objtop_deltaRmin = deltaR_obj2_fromTop_Top;
                 h_NotHs_objectsfromtop_top_mindR -> Fill(objtop_deltaRmin);
		 

                 if(deltaR_b_fromTop_Top < 0.8)  // b into fatjet
		   {
		     decayproductsintofatjet++;
		     bQuarkIntoFatJet = true;
		   }
		 else bQuarkIntoFatJet = false;
                 if(deltaR_obj1_fromTop_Top < 0.8) decayproductsintofatjet++;  // obj1 into fatjet                
                 if(deltaR_obj2_fromTop_Top < 0.8) decayproductsintofatjet++;  // obj2 into fatjet                               

                 h_Probdecayproductsintofatjet_NOTHs          -> Fill("0",0);  //just to determinate the label of the first bin  
                 h_Probdecayproductsintofatjet_NOTHs          -> Fill("1",0);  //just to determinate the label of the second bin       
                 h_Probdecayproductsintofatjet_NOTHs          -> Fill("2",0);  //just to determinate the label of the third bin             

		 h_NotHs_QuarksintofatjetMultiplicity -> Fill("0",0);  //just to determinate the label of the first bin 
                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("j",0);
                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("b",0);
                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("jj",0);
                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("bj",0);

                 h_NotHs_isbQuarkintofatjet           -> Fill("No b",0);  //just to determinate the label of the first bin                        

                 if      (decayproductsintofatjet == 0)
                   {
                     h_Probdecayproductsintofatjet_NOTHs  -> Fill("0",1);
                     h_NotHs_QuarksintofatjetMultiplicity -> Fill("0",1);
                     h_NotHs_isbQuarkintofatjet           -> Fill("No b",1);
                   }
                 else if (decayproductsintofatjet == 1)
                   {
                     h_Probdecayproductsintofatjet_NOTHs -> Fill("1",1);
                     if(bQuarkIntoFatJet)
                       {
                         h_NotHs_QuarksintofatjetMultiplicity -> Fill("b",1);
                         h_NotHs_isbQuarkintofatjet           -> Fill("b",1);
                       }
                     else
                       {
                         h_NotHs_QuarksintofatjetMultiplicity -> Fill("j",1);
                         h_NotHs_isbQuarkintofatjet           -> Fill("No b",1);
                       }
                   }
		 else if (decayproductsintofatjet == 2)
                   {
                     h_Probdecayproductsintofatjet_NOTHs -> Fill("2",1);
                     if(bQuarkIntoFatJet)
                       {
                         h_NotHs_QuarksintofatjetMultiplicity -> Fill("bj",1);
                         h_NotHs_isbQuarkintofatjet           -> Fill("b",1);
                       }
                     else
                       {
                         h_NotHs_QuarksintofatjetMultiplicity -> Fill("jj",1);
                         h_NotHs_isbQuarkintofatjet           -> Fill("No b",1);
		       }
                   }
                 else
                   {
                     h_Probdecayproductsintofatjet_NOTHs  -> Fill("All 3",1);
                     h_NotHs_QuarksintofatjetMultiplicity -> Fill("bjj",1);
                     h_NotHs_isbQuarkintofatjet           -> Fill("b",1);
                   }

		 if (decayproductsintofatjet == 3)
                   {
                     unsigned int Num_pairs_withenough_deltaR_amongthemselves = 0;
                     if      (deltaR_b_fromTop_obj1 < 0.4 && deltaR_b_fromTop_obj2 < 0.4 && deltaR_obj1_obj2 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 0;

                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 < 0.4 && deltaR_b_fromTop_obj2 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 1;
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 < 0.4 && deltaR_b_fromTop_obj1 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 1;
                     else if (deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj1 < 0.4 && deltaR_b_fromTop_obj2 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 1;

                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj2 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 2;
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj1 < 0.4) Num_pairs_withenough_deltaR_amongthemselves = 2;
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 <0.4 ) Num_pairs_withenough_deltaR_amongthemselves = 2;

                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 > 0.4) Num_pairs_withenough_deltaR_amongthemselves = 3;


                     h_Pairsinbarycenter_enoughdeltaR_NOTHs         -> Fill("0",0);  //just to determinate the label of the first bin     
                     h_Pairsinbarycenter_enoughdeltaR_NOTHs         -> Fill("1 pair",0);  //just to determinate the label of the second bin
                     h_Pairsinbarycenter_enoughdeltaR_NOTHs         -> Fill("2 pairs",0);  //just to determinate the label of the third bin 
		 
		     if      (Num_pairs_withenough_deltaR_amongthemselves == 0)   h_Pairsinbarycenter_enoughdeltaR_NOTHs -> Fill("0",1);
                     else if (Num_pairs_withenough_deltaR_amongthemselves == 1)   h_Pairsinbarycenter_enoughdeltaR_NOTHs -> Fill("1 pair",1);
                     else if (Num_pairs_withenough_deltaR_amongthemselves == 2)   h_Pairsinbarycenter_enoughdeltaR_NOTHs -> Fill("2 pairs",1);
                     else                                                         h_Pairsinbarycenter_enoughdeltaR_NOTHs -> Fill("All 3",1);
		     
		     h_Iffatjet_NotHs_Top_pT         -> Fill (Top_NOTfromH_p4.pt());
                     h_Iffatjet_NotHs_W_pT           -> Fill (W_fromTop_NOTfromH_p4.pt());
		     h_Iffatjet_NotHs_Top_pT_Vs_W_pT -> Fill(Top_NOTfromH_p4.pt(),W_fromTop_NOTfromH_p4.pt());
                     h_Iffatjet_NotHs_bfromH_pT      -> Fill (b_fromH_p4.pt());

		     h_Iffatjet_NotHs_EventsWithHighTop_pT -> Fill("top.pT < 400 GeV", 0);
		     h_Iffatjet_NotHs_EventsWithHighW_pT   -> Fill("W.pT < 200 GeV", 0);
                     if(Top_NOTfromH_p4.pt() > 400.0)       h_Iffatjet_NotHs_EventsWithHighTop_pT -> Fill("top.pT > 400 GeV", 1);
                     else                                   h_Iffatjet_NotHs_EventsWithHighTop_pT -> Fill("top.pT < 400 GeV", 1);
                     if(W_fromTop_NOTfromH_p4.pt() > 200.0) h_Iffatjet_NotHs_EventsWithHighW_pT   -> Fill("W.pT > 200 GeV", 1);
                     else                                   h_Iffatjet_NotHs_EventsWithHighW_pT   -> Fill("W.pT < 200 GeV", 1);
		     
		     if(Top_NOTfromH_p4.pt() > 400.0 && W_fromTop_NOTfromH_p4.pt() > 200.0)
                       {
                         h_Iffatjet_NotHs_HighWandTop_pT_dRmin -> Fill (deltaRmin);
                       }

                     double dR_bfromtop_top              = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, Top_NOTfromH_p4);
                     h_Iffatjet_NotHs_bfromT_Top_dR      -> Fill(dR_bfromtop_top);
                     double dEta_bfromtop_top            = std::abs(b_fromTop_NOTfromH_p4.eta() - Top_NOTfromH_p4.eta());
                     h_Iffatjet_NotHs_bfromT_Top_dEta    -> Fill(dEta_bfromtop_top);
                     double dPhi_bfromtop_top            = std::abs(b_fromTop_NOTfromH_p4.phi() - Top_NOTfromH_p4.phi());
                     h_Iffatjet_NotHs_bfromT_Top_dPhi    -> Fill(dPhi_bfromtop_top);

                     double dR_bfromtop_W                = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, W_fromTop_NOTfromH_p4);
                     h_Iffatjet_NotHs_bfromT_W_dR        -> Fill(dR_bfromtop_W);
                     double dEta_bfromtop_W              = std::abs(b_fromTop_NOTfromH_p4.eta() - W_fromTop_NOTfromH_p4.eta());
                     h_Iffatjet_NotHs_bfromT_W_dEta      -> Fill(dEta_bfromtop_W);
                     double dPhi_bfromtop_W              = std::abs(b_fromTop_NOTfromH_p4.phi() - W_fromTop_NOTfromH_p4.phi());
                     h_Iffatjet_NotHs_bfromT_W_dPhi      -> Fill(dPhi_bfromtop_W);

                     double dR_bfromtop_bfromH           = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, b_fromH_p4);
                     h_Iffatjet_NotHs_bfromT_bfromH_dR   -> Fill(dR_bfromtop_bfromH);
                     double dEta_bfromtop_bfromH         = std::abs(b_fromTop_NOTfromH_p4.eta() - b_fromH_p4.eta());
                     h_Iffatjet_NotHs_bfromT_bfromH_dEta -> Fill(dEta_bfromtop_bfromH);
                     double dPhi_bfromtop_bfromH         = std::abs(b_fromTop_NOTfromH_p4.phi() - b_fromH_p4.phi());
                     h_Iffatjet_NotHs_bfromT_bfromH_dPhi -> Fill(dPhi_bfromtop_bfromH);
		     
		     
                   } // All 3 decay products into 0.8

		 // try with gen-jets ++++Barycenter+++++++++++++++++++++++++++++++++++++++++
                 if(numof_obj_fromtop_matchedwith_genjet == 3)
                   {
                     double deltaR_bjet_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(matched_genbJet_p4, Top_NOTfromH_p4);
                     double deltaR_jet1_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(genJets_from_W_p4.at(0), Top_NOTfromH_p4);
                     double deltaR_jet2_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(genJets_from_W_p4.at(1), Top_NOTfromH_p4);
                     unsigned int jets_intofatjet = 0;

                     if(deltaR_bjet_fromTop_Top < 0.8) 
		       {
			 jets_intofatjet++;// bjet into fatjet  
			 bJetInsideFatJet = true;
		       }
		     else bJetInsideFatJet = false;
                     if(deltaR_jet1_fromTop_Top < 0.8) jets_intofatjet++;// jet1 into fatjet  
                     if(deltaR_jet2_fromTop_Top < 0.8) jets_intofatjet++;// jet2 into fatjet    

		     //labels
                     h_NotHs_ProbdecayJetsintofatjet   -> Fill("0",0);  //just to determinate the label of the first bin    
                     h_NotHs_ProbdecayJetsintofatjet   -> Fill("1",0);  //just to determinate the label of the second bin          
                     h_NotHs_ProbdecayJetsintofatjet   -> Fill("2",0);  //just to determinate the label of the third bin 

		     h_NotHs_JetsintofatjetMultiplicity-> Fill("0",0);  //just to determinate the label of the first bin   
                     h_NotHs_JetsintofatjetMultiplicity-> Fill("j",0);
                     h_NotHs_JetsintofatjetMultiplicity-> Fill("b",0);
                     h_NotHs_JetsintofatjetMultiplicity-> Fill("jj",0);
                     h_NotHs_JetsintofatjetMultiplicity-> Fill("bj",0);

                     h_NotHs_isbJetintofatjet          -> Fill("No b-jet",0);  //just to determinate the label of the first bin  
		     //

                     if      (jets_intofatjet == 0)
                       {
                         h_NotHs_ProbdecayJetsintofatjet     -> Fill("0",1);
                         h_NotHs_JetsintofatjetMultiplicity  -> Fill("0",1);
                         h_NotHs_isbJetintofatjet            -> Fill("No b-jet",1);
                       }

                     else if (jets_intofatjet == 1)
                       {
                         h_NotHs_ProbdecayJetsintofatjet -> Fill("1",1);
                         if(bJetInsideFatJet)
                           {
                             h_NotHs_JetsintofatjetMultiplicity-> Fill("b",1);
                             h_NotHs_isbJetintofatjet          -> Fill("b-jet",1);
                           }
                         else
                           {
                             h_NotHs_JetsintofatjetMultiplicity-> Fill("j",1);
                             h_NotHs_isbJetintofatjet          -> Fill("No b-jet",1);
                           }
                       }

                     else if (jets_intofatjet == 2)
                       {
                         h_NotHs_ProbdecayJetsintofatjet -> Fill("2",1);
                         if(bJetInsideFatJet)
                           {
                             h_NotHs_JetsintofatjetMultiplicity-> Fill("bj",1);
                             h_NotHs_isbJetintofatjet          -> Fill("b-jet",1);
                           }
                         else
                           {
                             h_NotHs_JetsintofatjetMultiplicity-> Fill("jj",1);
                             h_NotHs_isbJetintofatjet          -> Fill("No b-jet",1);
                           }
                       }

                     else
                       {
                         h_NotHs_ProbdecayJetsintofatjet -> Fill("All 3",1);
                         h_NotHs_JetsintofatjetMultiplicity-> Fill("bjj",1);
                         h_NotHs_isbJetintofatjet          -> Fill("b-jet",1);
                       }
		     
		     h_NotHs_BjetInsideFatJet_Top_Pt -> Fill(bJetInsideFatJet, Top_genReco_NOTfromH_p4.pt());
                     h_NotHs_BjetInsideFatJet_W_Pt   -> Fill(bJetInsideFatJet, W_genReco_fromTop_NOTfromH_p4.pt());
                   }
                 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                 double deltaR_bNOTfromHiggstop_top = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, p.p4()); //find dR(bNOTfromtfromH,tNOTfromH)                      
		 h_bNOTfromHiggstop_top_dR                  -> Fill(deltaR_bNOTfromHiggstop_top);
		 
		 h_b_topNOTfromhiggs_underdR                -> Fill("< 0.8",0);  //just to determinate the label of the first bin             
                 if(deltaR_bNOTfromHiggstop_top < 0.8)    h_b_topNOTfromhiggs_underdR -> Fill("< 0.8",1);
                 else                                     h_b_topNOTfromhiggs_underdR -> Fill("> 0.8",1);
		 
		 h_bNOTfromHiggstop_top_dR_Vs_top_Pt        -> Fill(deltaR_bNOTfromHiggstop_top,p.pt());
		 
		 ////////////////////////////////////////////////////////////////////////////////////////         
                 double deltaR_WNOTfromHiggs_closest_b = 9999.9; //non sense initial value                            
                 for (auto& gp: fEvent.genparticles().getGenParticles()) // for genpar %3                    
                   {

                     if(std::abs(gp.pdgId()) == 5 && gp.isFirstCopy()) // find the b                                      
                       {
                         double deltaR_WNOTfromHiggs_any_b = ROOT::Math::VectorUtil::DeltaR(W_fromTop_NOTfromH_p4,gp.p4());
                         if ( deltaR_WNOTfromHiggs_any_b < deltaR_WNOTfromHiggs_closest_b) deltaR_WNOTfromHiggs_closest_b = deltaR_WNOTfromHiggs_any_b;
                         else continue;
                       }
                   }// for genpar %3                                                  
                 h_WNOTfromHiggstop_closestb_dR                            -> Fill(deltaR_WNOTfromHiggs_closest_b);
                 h_WNOTfromHiggstop_closestb_dR_Vs_W_Pt                    -> Fill(deltaR_WNOTfromHiggs_closest_b,W_fromTop_NOTfromH_p4.pt());
                 h_WNOTfromHiggstop_closestb_dR_Vs_samedecaytop_Pt  -> Fill(deltaR_WNOTfromHiggs_closest_b,Top_NOTfromH_p4.pt());
		 ///////////////////////////////////////////////////////////////////////////////////////////
		 
	       } //else the top is not from Higgs
	     
	   }     //if top
       }         //for gen particles
   } 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 //////kchristo////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 ////////////////////////////////////////////////////// Study Boosted Topologies///////////////////////////////////////////////////////////// 
 /////////////////////////////////////////////////////////With fat Jets//////////////////////////////////////////////////////////////////////
 if(1)
   {
     for (auto& p: fEvent.genparticles().getGenParticles())
       {

         if(std::abs(p.pdgId()) == 6 && p.isFirstCopy()) // find the top   
           {
	     std::vector<short> top_mothers = p.mothers();
             genParticle m = fEvent.genparticles().getGenParticles()[top_mothers.at(0)]; //create the particle object

	     if(std::abs(m.pdgId()) == 37)//if the top comes from Higgs//////////////////////////////////////////////////////////////
               {
		 math::XYZTLorentzVector b_fromTop_fromH_p4(0,0,0,0), obj1_fromW_fromTop_p4(0,0,0,0), obj2_fromW_fromTop_p4(0,0,0,0);
		 math::XYZTLorentzVector W_fromTop_fromH_p4(0,0,0,0), Top_fromH_p4(0,0,0,0);
                 bool isTheRightTop    = true;
                 bool obj1_exist       = false;
                 bool bQuarkIntoFatJet = false;
		 unsigned int b_fromTop_fromH_index = 0;

                 Top_fromH_p4 = p.p4();
		 
		 unsigned int lastcopy = GetTheLastCopy(p.index()); // we need to find the last copy in order to find the real daughters 
                 genParticle lastT =  fEvent.genparticles().getGenParticles()[lastcopy]; //create the particle object of the last copy of top
		 std::vector<short> lastT_daughters = lastT.daughters();

                 for (unsigned int i = 0; i < lastT_daughters.size() ; i++)
                   {
                     genParticle d; //create the particle object
                     d =  fEvent.genparticles().getGenParticles()[lastT_daughters.at(i)];

                     if (std::abs(d.pdgId()) == 24) // if W from top  
                       {
                         W_fromTop_fromH_p4 = d.p4();
                         genParticle lastW;
                         unsigned int lastWcopy = GetTheLastCopy(d.index()); // we need to find the last copy in order to find the real daughters
			 lastW = fEvent.genparticles().getGenParticles()[lastWcopy];
			 std::vector<short> W_daughters = lastW.daughters();
                         for (size_t k = 0; k < W_daughters.size() ; k++)
                           {
                             genParticle top_grand_d; //create the particle object, i define this as grand_d, the daus of the W, so the grand-daughters of Top     
			     top_grand_d =  fEvent.genparticles().getGenParticles()[W_daughters.at(k)];
			     
			     if ( !mcTools.IsQuark(top_grand_d.pdgId()) ) //if the daus of W are not quarks (f.e. leptons), we reject the top 
			       {
				 isTheRightTop = false;
				 continue;
			       }
			     
			     if (obj1_exist == false)
			       {
				 obj1_exist = true;
				 obj1_fromW_fromTop_p4 = top_grand_d.p4();
			       }
			     else
                               {
                                 obj2_fromW_fromTop_p4 = top_grand_d.p4();
                               }

                           }      //for loop, W's daus                            
                         if (!isTheRightTop) continue; // move out from the checking loop, reject this top    

                       } // if W from top   

                     if (std::abs(d.pdgId()) == 5) // if b from top 
                       {
                         b_fromTop_fromH_p4 = d.p4();
			 b_fromTop_fromH_index = d.index();
                       }
                   }//for top from Higgs daus

		 if (!isTheRightTop) continue; // move out from the checking loop, reject this top 
		 
		 // Matching top and W with fat-jet++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                 bool TheTopQuarkIsMatched_WithFatJet = false;
		 bool TheWQuarkIsMatched_WithFatJet = false;

                 double deltaRmin_fatJet_topQuark = 1e6; //give an initial, non sense value  
		 double deltaRmin_fatJet_WQuark = 1e6; //give an initial, non sense value
		 math::XYZTLorentzVector closest_toTop_fatJet_p4(0,0,0,0), closest_toW_fatJet_p4(0,0,0,0);
                 int AK8jetSD_index = -1;
                 for(AK8JetsSoftDrop fatjet: fEvent.ak8jetsSoftDrop())
                   {
                     AK8jetSD_index++;
		     math::XYZTLorentzVector fatJet_p4(0,0,0,0);
                     fatJet_p4 = fatjet.p4();
                     double deltaR_fatJet_topQuark  = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,Top_fromH_p4);
		     double deltaR_fatJet_WQuark    = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,W_fromTop_fromH_p4);
                     if(deltaR_fatJet_topQuark < deltaRmin_fatJet_topQuark)
                       {
                         deltaRmin_fatJet_topQuark = deltaR_fatJet_topQuark;
                         closest_toTop_fatJet_p4 = fatjet.p4();
                       }
		     if(deltaR_fatJet_WQuark < deltaRmin_fatJet_WQuark)
                       {
                         deltaRmin_fatJet_WQuark = deltaR_fatJet_WQuark;
                         closest_toW_fatJet_p4 = fatjet.p4();
                       }
                   }   //for fat-jets loop
		 
                 if(deltaRmin_fatJet_topQuark <0.8) TheTopQuarkIsMatched_WithFatJet = true;
		 if(deltaRmin_fatJet_WQuark <0.8)   TheWQuarkIsMatched_WithFatJet = true;
                 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		 double deltaR_b_fromTop_obj1 = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, obj1_fromW_fromTop_p4);
                 double deltaR_b_fromTop_obj2 = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, obj2_fromW_fromTop_p4);
                 double deltaR_obj1_obj2      = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, obj2_fromW_fromTop_p4);
                 double deltaRmax = -1;  // give an initial, non sense value           
                 double deltaRmin = 1e6; //give an initial, non sense value         

		 // find the dRmax 
                 if(deltaR_b_fromTop_obj1 > deltaR_b_fromTop_obj2 && deltaR_b_fromTop_obj1 > deltaR_obj1_obj2)   
		   {
		     deltaRmax = deltaR_b_fromTop_obj1;
		   }
                 else if(deltaR_b_fromTop_obj2 > deltaR_b_fromTop_obj1 && deltaR_b_fromTop_obj2 > deltaR_obj1_obj2) 
		   {
		     deltaRmax = deltaR_b_fromTop_obj2;
		   }
                 else deltaRmax = deltaR_obj1_obj2;

                 h_objectsfromHiggstop_maxdR                  -> Fill(deltaRmax);

                 //find the dRmin 
                 if(deltaR_b_fromTop_obj1 < deltaR_b_fromTop_obj2 && deltaR_b_fromTop_obj1 < deltaR_obj1_obj2)      
		   {
		     deltaRmin = deltaR_b_fromTop_obj1;
		   }
                 else if(deltaR_b_fromTop_obj2 < deltaR_b_fromTop_obj1 && deltaR_b_fromTop_obj2 < deltaR_obj1_obj2) 
		   {
		     deltaRmin = deltaR_b_fromTop_obj2;
		   }
                 else deltaRmin = deltaR_obj1_obj2;

                 h_objectsfromHiggstop_mindR                                -> Fill(deltaRmin);
		 
		 h_objectsfromHiggstop_Prob_mindR_lt_p8                     -> Fill("<0.8",0);
                 if(deltaRmin < 0.8)
                   {
		     h_Hs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj -> Fill("Matched",0);
		     if(TheTopQuarkIsMatched_WithFatJet) 
		       {
			 h_Hs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj -> Fill("Matched",1);
			 h_Hs_objectsfromtop_mindR_ltp8_matchedfj_pt -> Fill(closest_toTop_fatJet_p4.pt());
		       }
		     else                                h_Hs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj -> Fill("Not Matched",1);
		       
                     h_objectsfromHiggstop_Prob_mindR_lt_p8 -> Fill("<0.8",1);
                     h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("< 400",0);
                     h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("400<pT<500",0);
                     if     (Top_fromH_p4.pt() > 500.0) 
		       {
			 h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("> 500",1);
			 if(TheTopQuarkIsMatched_WithFatJet)
			   {
			     h_Hs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt -> Fill(closest_toTop_fatJet_p4.pt());
			   }
		       }
                     else if(Top_fromH_p4.pt() > 400.0) h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("400<pT<500",1);
                     else                               h_objectsfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("< 400",1);
		     
		     if     (Top_fromH_p4.pt() > 400.0)
                       {
			 if(TheTopQuarkIsMatched_WithFatJet)
                           {
                             h_Hs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt -> Fill(closest_toTop_fatJet_p4.pt());
			   }
                       }
                   }

                 else  h_objectsfromHiggstop_Prob_mindR_lt_p8               -> Fill(">0.8",1);

                 h_objectsfromHiggstop_mindR_Vs_Higgstop_Pt                 -> Fill(deltaRmin,p.pt());
                 h_objectsfromHiggstop_mindR_Vs_Higgs_Pt                    -> Fill(deltaRmin,m.pt());
		 
		 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                 /////////////////////////////////////////////////BaryCenter//////////////////////////////////////////////////////////////// 
		 /////////////assuming that is the top-direction ///////////////////////////////////////////////////////////////////////////

                 double deltaR_b_fromTop_Top    = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4, Top_fromH_p4);
                 double deltaR_obj1_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, Top_fromH_p4);
                 double deltaR_obj2_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(obj2_fromW_fromTop_p4, Top_fromH_p4);
                 unsigned int decayproductsintofatjet = 0; // into fat jet, means the case where the products are in 0.8 from top-direction
                 double objtop_deltaRmax = -1;  // give an initial, non sense value  
                 double objtop_deltaRmin = 999;  // give an initial, non sense value        

                 h_Hs_which_objectfromtop_maxdR -> Fill ("W'obj",0) ; //just to determinate the first label      

                 if(deltaR_b_fromTop_Top > deltaR_obj1_fromTop_Top && deltaR_b_fromTop_Top > deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmax = deltaR_b_fromTop_Top;
                     h_Hs_which_objectfromtop_maxdR -> Fill("b",1) ;
                     if(deltaR_b_fromTop_Top > 0.8)
                       {
                         h_Hs_mostdistantfromtop_isb__top_pT      -> Fill(Top_fromH_p4.pt());
                         h_Hs_mostdistantfromtop_isb__W_pT        -> Fill(W_fromTop_fromH_p4.pt());
                         h_Hs_mostdistantfromtop_isb__b_pT        -> Fill(b_fromTop_fromH_p4.pt());
                         h_Hs_mostdistantfromtop_isb__dRmax_b_top -> Fill(deltaR_b_fromTop_Top);
			 
			 if(deltaR_b_fromTop_obj1 < deltaR_b_fromTop_obj2) 
			   {
			     h_Hs_mostdistantfromtop_isb_dRmin_b_objfromW -> Fill(deltaR_b_fromTop_obj1);
			   }
                         else   
			   {
			     h_Hs_mostdistantfromtop_isb_dRmin_b_objfromW -> Fill(deltaR_b_fromTop_obj2);
			   }
                       }

                     if(deltaR_b_fromTop_Top > 0.8 && deltaR_obj1_obj2 < 0.7)
                       {
			 h_Hs_mostdistantfromtop_isb__dRqq_Vs_W_pT   -> Fill(deltaR_obj1_obj2,W_fromTop_fromH_p4.pt());
                         h_Hs_mostdistantfromtop_isb__dRqq_Vs_top_pT -> Fill(deltaR_obj1_obj2,Top_fromH_p4.pt());

			 h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj -> Fill("Matched",0);
			 if(TheTopQuarkIsMatched_WithFatJet)
			   {
			     h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj -> Fill("Matched",1);
			     h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt -> Fill(closest_toTop_fatJet_p4.pt());
			   }
			 else h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj -> Fill("Not Matched",1);
			 
			 h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj -> Fill("Matched",0);
                         if(TheWQuarkIsMatched_WithFatJet)
                           {
                             h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj -> Fill("Matched",1);
                             h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt -> Fill(closest_toW_fatJet_p4.pt());
                           }
                         else h_Hs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj -> Fill("Not Matched",1);
			 

                         h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("< 400",0);
                         h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("400<pT<500",0);
                         h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("< 200",0);
                         h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("200<pT<300",0);
                        
			 if     (Top_fromH_p4.pt() > 500.0) 
			   {
			     h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("> 500",1);
			     if(TheTopQuarkIsMatched_WithFatJet)
			       {
				 h_Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt ->Fill(closest_toTop_fatJet_p4.pt());
			       }
			   }
                         else if(Top_fromH_p4.pt() > 400.0)
			   {
			     h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("400<pT<500",1);
			   }
                         else h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("< 400",1);
			 
			 if     (Top_fromH_p4.pt() > 400.0)
			   {
			     if(TheTopQuarkIsMatched_WithFatJet)
			       {
                                 h_Hs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt ->Fill(closest_toTop_fatJet_p4.pt());
                               }
                           }
			 

                         if     (W_fromTop_fromH_p4.pt() > 300.0)
			   {
			     h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("> 300",1);
			     if(TheWQuarkIsMatched_WithFatJet)
			       {
                                 h_Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt ->Fill(closest_toW_fatJet_p4.pt());
                               }
			   }
                         else if(W_fromTop_fromH_p4.pt() > 200.0) 
			   {
			     h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("200<pT<300",1);
			   }
                         else h_Hs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("< 200",1);
			 
			 if     (W_fromTop_fromH_p4.pt() > 200.0)
                           {
                             if(TheWQuarkIsMatched_WithFatJet)
                               {
                                 h_Hs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt ->Fill(closest_toW_fatJet_p4.pt());
                               }
                           }
                       }

                     h_Hs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("True",0);
                     if(deltaR_b_fromTop_Top > 0.8 && deltaR_obj1_obj2 < 0.7) h_Hs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("True",1);
                     else                                                     h_Hs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("False",1);

                   }
                 else if(deltaR_obj1_fromTop_Top > deltaR_b_fromTop_Top && deltaR_obj1_fromTop_Top > deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmax = deltaR_obj1_fromTop_Top;
                     h_Hs_which_objectfromtop_maxdR-> Fill("W'obj",1);
                   }
                 else
                   {
                     objtop_deltaRmax = deltaR_obj2_fromTop_Top;
                     h_Hs_which_objectfromtop_maxdR-> Fill("W'obj",1);
                   }
                 h_Hs_objectsfromtop_top_maxdR -> Fill(objtop_deltaRmax);

                 if(deltaR_b_fromTop_Top < deltaR_obj1_fromTop_Top && deltaR_b_fromTop_Top < deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmin = deltaR_b_fromTop_Top;
                   }
		 else if(deltaR_obj1_fromTop_Top < deltaR_b_fromTop_Top && deltaR_obj1_fromTop_Top < deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmin = deltaR_obj1_fromTop_Top;
                   }
                 else objtop_deltaRmin = deltaR_obj2_fromTop_Top;
                 h_Hs_objectsfromtop_top_mindR -> Fill(objtop_deltaRmin);

		 if(deltaR_b_fromTop_Top < 0.8)    // b into fatjet  
                   {
                     decayproductsintofatjet++;
                     bQuarkIntoFatJet = true;
                   }
                 else bQuarkIntoFatJet = false;
                 if(deltaR_obj1_fromTop_Top < 0.8) decayproductsintofatjet++;// obj1 into fatjet                
                 if(deltaR_obj2_fromTop_Top < 0.8) decayproductsintofatjet++;// obj2 into fatjet      
		 // into fat jet, means the case where the products are in 0.8 from top-direction

                 h_Probdecayproductsintofatjet_Hs         -> Fill("0",0);  //just to determinate the label of the first bin    
                 h_Probdecayproductsintofatjet_Hs         -> Fill("1",0);  //just to determinate the label of the second bin       
                 h_Probdecayproductsintofatjet_Hs         -> Fill("2",0);  //just to determinate the label of the third bin          

                 h_Hs_QuarksintofatjetMultiplicity -> Fill("0",0);  //just to determinate the label of the first bin                          
                 h_Hs_QuarksintofatjetMultiplicity -> Fill("q",0);
                 h_Hs_QuarksintofatjetMultiplicity -> Fill("b",0);
                 h_Hs_QuarksintofatjetMultiplicity -> Fill("qq",0);
                 h_Hs_QuarksintofatjetMultiplicity -> Fill("bq",0);

                 h_Hs_isbQuarkintofatjet           -> Fill("No b",0);  //just to determinate the label of the first bin      

                 if      (decayproductsintofatjet == 0)
                   {
                     h_Probdecayproductsintofatjet_Hs  -> Fill("0",1);
                     h_Hs_QuarksintofatjetMultiplicity -> Fill("0",1);
                     h_Hs_isbQuarkintofatjet           -> Fill("No b",1);
                   }
		 else if (decayproductsintofatjet == 1)
                   {
                     h_Probdecayproductsintofatjet_Hs -> Fill("1",1);
                     if(bQuarkIntoFatJet)
                       {
                         h_Hs_QuarksintofatjetMultiplicity -> Fill("b",1);
                         h_Hs_isbQuarkintofatjet           -> Fill("b",1);
                       }
                     else
                       {
                         h_Hs_QuarksintofatjetMultiplicity -> Fill("q",1);
                         h_Hs_isbQuarkintofatjet           -> Fill("No b",1);
                       }
                   }
                 else if (decayproductsintofatjet == 2)
                   {
                     h_Probdecayproductsintofatjet_Hs -> Fill("2",1);
                     if(bQuarkIntoFatJet)
                       {
                         h_Hs_QuarksintofatjetMultiplicity -> Fill("bq",1);
                         h_Hs_isbQuarkintofatjet           -> Fill("b",1);
                       }
		     else
                       {
                         h_Hs_QuarksintofatjetMultiplicity -> Fill("qq",1);
                         h_Hs_isbQuarkintofatjet           -> Fill("No b",1);
                       }
                   }
                 else
                   {
                     h_Probdecayproductsintofatjet_Hs  -> Fill("All 3",1);
                     h_Hs_QuarksintofatjetMultiplicity -> Fill("bq",1);
                     h_Hs_isbQuarkintofatjet           -> Fill("b",1);
                   }
		 
		 if (decayproductsintofatjet == 3)
                   {
                     unsigned int Num_pairs_withenough_deltaR_amongthemselves = 0;
                     if      (deltaR_b_fromTop_obj1 < 0.4 && deltaR_b_fromTop_obj2 < 0.4 && deltaR_obj1_obj2 < 0.4) 
		       {
			 Num_pairs_withenough_deltaR_amongthemselves = 0;
		       }

                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 < 0.4 && deltaR_b_fromTop_obj2 < 0.4) 
		       {
			 Num_pairs_withenough_deltaR_amongthemselves = 1;
		       }
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 < 0.4 && deltaR_b_fromTop_obj1 < 0.4) 
		       {
			 Num_pairs_withenough_deltaR_amongthemselves = 1;
		       }
                     else if (deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj1 < 0.4 && deltaR_b_fromTop_obj2 < 0.4) 
		       {
			 Num_pairs_withenough_deltaR_amongthemselves = 1;
		       }
		     
                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj2 < 0.4) 
		       {
			 Num_pairs_withenough_deltaR_amongthemselves = 2;
		       }
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj1 < 0.4) 
		       {
			 Num_pairs_withenough_deltaR_amongthemselves = 2;
		       }
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 <0.4 )
		       {
			 Num_pairs_withenough_deltaR_amongthemselves = 2;
		       }

                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 > 0.4) 
		       {
			 Num_pairs_withenough_deltaR_amongthemselves = 3;
		       }
		     
		     h_Pairsinbarycenter_enoughdeltaR_Hs         -> Fill("0",0);  //just to determinate the label of the first bin 
                     h_Pairsinbarycenter_enoughdeltaR_Hs         -> Fill("1 pair",0);  //just to determinate the label of the second bin      
                     h_Pairsinbarycenter_enoughdeltaR_Hs         -> Fill("2 pairs",0);  //just to determinate the label of the third bin 

                     if      (Num_pairs_withenough_deltaR_amongthemselves == 0)   h_Pairsinbarycenter_enoughdeltaR_Hs -> Fill("0",1);
                     else if (Num_pairs_withenough_deltaR_amongthemselves == 1)   h_Pairsinbarycenter_enoughdeltaR_Hs -> Fill("1 pair",1);
                     else if (Num_pairs_withenough_deltaR_amongthemselves == 2)   h_Pairsinbarycenter_enoughdeltaR_Hs -> Fill("2 pairs",1);
                     else                                                         h_Pairsinbarycenter_enoughdeltaR_Hs -> Fill("All 3",1);

                     h_Iffatjet_Hs_Top_pT               -> Fill (Top_fromH_p4.pt());
		     h_Iffatjet_Hs_EventsWithHighTop_pT -> Fill("top.pT < 400 GeV", 0);
                     if(Top_fromH_p4.pt() > 400.0)       h_Iffatjet_Hs_EventsWithHighTop_pT -> Fill("top.pT > 400 GeV", 1);
                     else                                h_Iffatjet_Hs_EventsWithHighTop_pT -> Fill("top.pT < 400 GeV", 1);
		 
		   } // if all 3 decay products are into 0.8
		 
		 /////////////////////////////////////////////////BaryCenter////////////////////////////////////////////////////////////////
		 // barycenter not the top-direction.  Actually here take place a different approach.  If at least 2 objects have dr<0.8,  
		 // had been assumed that they are into the barycenter.  After that, comes the matching with AK8 jets///////////////////////
		 
		 h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("0",0);  //just to determinate the label of the first bin
		 h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("qq",0);
                 h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("bq",0);
		 h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("qq-bq",0);

                 h_Hs_isbQuarkintoBaryCenter           -> Fill("No b",0);  //just to determinate the label of the first bin
		 h_Hs_isbQuarkintoBaryCenter           -> Fill("?",0);
		 
		 h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("0",0);  //just to determinate the label of the first bin     
                 h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq",0);
                 h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bq",0);
                 h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq-bq",0);

                 h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",0);  //just to determinate the label of the first bin 
                 h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("?",0);
		 
		 int otherBcloseToTopProd = 0;
		 std::vector<bool> otherBintoBoostTop(2);
		 std::vector<bool> otherBintoBoostW(2);
		 std::vector<math::XYZTLorentzVector> otherBcloseToTopProd_p4(2);
		 math::XYZTLorentzVector p4_initializer(0,0,0,0);
		 for(size_t otherb=0; otherb < otherBcloseToTopProd_p4.size(); otherb++)
		   {
		     otherBcloseToTopProd_p4.at(otherb) = p4_initializer;
		     otherBintoBoostTop.at(otherb) = false;
		     otherBintoBoostW.at(otherb) = false;
		   }
		 
		 for (unsigned int applyCut = 0; applyCut < 2 ; applyCut++) // just to find the fraction of events that are passing the quarksPtCut = 30.0 and the differences on the shapes of various histograms		                                                            
		   {
		     double quarksPtCut = 30.0;
		     if(applyCut)
		       {
			 h_Hs_QuarksFromTop_Passed_pTcuts -> Fill("Passed",0);
			 if((obj1_fromW_fromTop_p4.pt() < quarksPtCut) || (obj2_fromW_fromTop_p4.pt() < quarksPtCut) ||  (b_fromTop_fromH_p4.pt() < quarksPtCut)) 
			   {
			     h_Hs_QuarksFromTop_Passed_pTcuts -> Fill("Not Passed",1);
			     continue;
			   }
			 else h_Hs_QuarksFromTop_Passed_pTcuts -> Fill("Passed",1);
			 
			 for (auto& gp: fEvent.genparticles().getGenParticles()) // for genpar %4  
			   {
			     if(std::abs(gp.pdgId()) == 5 && gp.isFirstCopy()) // find the b
			       {
				 double deltaR_bfromTop_any_b         = ROOT::Math::VectorUtil::DeltaR(b_fromTop_fromH_p4,gp.p4());
				 double deltaR_obj1fromWfromTop_any_b = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4,gp.p4());
				 double deltaR_obj2fromWfromTop_any_b = ROOT::Math::VectorUtil::DeltaR(obj2_fromW_fromTop_p4,gp.p4());
				 if(b_fromTop_fromH_index == gp.index() || gp.pt() < quarksPtCut) continue; //dont take the b you have or a soft one
				 else if(deltaR_bfromTop_any_b > 0.8 && deltaR_obj1fromWfromTop_any_b > 0.8 && deltaR_obj2fromWfromTop_any_b > 0.8) continue;
				 else 
				   {
				     
				     if (deltaRmax < 0.8) //boosted top
				       {
					 if(deltaR_bfromTop_any_b < 0.8 && deltaR_obj1fromWfromTop_any_b < 0.8 && deltaR_obj2fromWfromTop_any_b) otherBintoBoostTop.at(otherBcloseToTopProd) = true; 
					 //h_Hs_otherBclose_BoostedTop -> Fill("bqq->bbqq",0);
					 //if(deltaR_bfromTop_any_b < 0.8 && deltaR_obj1fromWfromTop_any_b < 0.8 && deltaR_obj2fromWfromTop_any_b) h_Hs_otherBclose_BoostedTop -> Fill("bqq->bbqq",1);
					 //else h_Hs_otherBclose_BoostedTop -> Fill("bqq",1);
				       }
				     else if(deltaR_obj1_obj2 < 0.8 && deltaR_b_fromTop_obj1 > 0.8 && deltaR_b_fromTop_obj2 > 0.8) //BoostedW
				       {
					 if(deltaR_obj1fromWfromTop_any_b < 0.8 && deltaR_obj2fromWfromTop_any_b) otherBintoBoostW.at(otherBcloseToTopProd) = true;
					 //h_Hs_otherBclose_BoostedW -> Fill("qq -> bqq",0);
                                         //if(deltaR_obj1fromWfromTop_any_b < 0.8 && deltaR_obj2fromWfromTop_any_b) h_Hs_otherBclose_BoostedW -> Fill("qq -> bqq",1);
                                         //else h_Hs_otherBclose_BoostedW -> Fill("qq",1);
				       }
				      
				     otherBcloseToTopProd++;
				     if(otherBcloseToTopProd == 2)  break;
				   }
			       }
			   }// for genpar %4 
			 h_Hs_otherBcloseToTopProd -> Fill(otherBcloseToTopProd);
		       }
		     
		     if(applyCut)
		       {
			 h_Hs_objectsfromtop_dRmax_pTcuts        -> Fill(deltaRmax);
			 h_Hs_objectsfromtop_Prob_dRmax_pTcuts        -> Fill("< 0.8",0);
                         if(deltaRmax < 0.8) h_Hs_objectsfromtop_Prob_dRmax_pTcuts -> Fill("< 0.8",1);
                         else                h_Hs_objectsfromtop_Prob_dRmax_pTcuts -> Fill("> 0.8",1);

			 h_Hs_objectsfromtop_dRmin_pTcuts        -> Fill(deltaRmin);
			 h_Hs_objectsfromtop_Prob_dRmin_pTcuts        -> Fill("< 0.8",0);
                         if(deltaRmin < 0.8) h_Hs_objectsfromtop_Prob_dRmin_pTcuts -> Fill("< 0.8",1);
                         else                h_Hs_objectsfromtop_Prob_dRmin_pTcuts -> Fill("> 0.8",1);
			 

			 h_Hs_QuarksFromW_deltaR_pTcuts          -> Fill(deltaR_obj1_obj2);
			 h_Hs_QuarksFromW_Prob_deltaR_pTcuts     -> Fill("< 0.8",0);
			 if(deltaR_obj1_obj2 < 0.8) h_Hs_QuarksFromW_Prob_deltaR_pTcuts -> Fill("< 0.8",1);
			 else                       h_Hs_QuarksFromW_Prob_deltaR_pTcuts -> Fill("> 0.8",1);
		       }
		     else // if the ptCut on the quarks hasn't been applied
		       {
			 h_Hs_QuarksFromW_deltaR                                -> Fill(deltaR_obj1_obj2);
			 h_Hs_QuarksFromW_Prob_deltaR-> Fill("< 0.8",0);
			 if(deltaR_obj1_obj2 < 0.8) h_Hs_QuarksFromW_Prob_deltaR-> Fill("< 0.8",1);
			 else                       h_Hs_QuarksFromW_Prob_deltaR-> Fill("> 0.8",1);
		       }

		     if (deltaRmax < 0.8)             
		       {
			 if(applyCut)
			   {
			     h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bqq",1);
			     h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("b",1);
			     h_Hs_Bc_bqq_Top_pT                           -> Fill(Top_fromH_p4.pt());
			     h_Hs_Bc_bqq_W_pT                             -> Fill(W_fromTop_fromH_p4.pt());
			     
			     h_Hs_otherBclose_BoostedTop -> Fill("bqq->bbqq",0);                                                 
			     if(otherBintoBoostTop.at(0) || otherBintoBoostTop.at(1)) h_Hs_otherBclose_BoostedTop -> Fill("bqq->bbqq",1);    
			     else h_Hs_otherBclose_BoostedTop -> Fill("bqq",1);
			   }
			 else// if the ptCut on the quarks hasn't been applied
			   {
			     h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("bqq",1);
			     h_Hs_isbQuarkintoBaryCenter           -> Fill("b",1);
			   }
		       }
		     else if(deltaR_obj1_obj2 < 0.8 && deltaR_b_fromTop_obj1 > 0.8 && deltaR_b_fromTop_obj2 > 0.8)
		       {
			 if(applyCut)
                           {
			     h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq",1);
			     h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",1);

			     h_Hs_BoostedW_deltaR_qq      -> Fill(deltaR_obj1_obj2);
			     h_Hs_BoostedW_Prob_deltaR_lt -> Fill("< 0.4",0);
			     if(deltaR_obj1_obj2 < 0.4) h_Hs_BoostedW_Prob_deltaR_lt -> Fill("< 0.4",1);
			     else                       h_Hs_BoostedW_Prob_deltaR_lt -> Fill("> 0.4",1);
			     
                             h_Hs_otherBclose_BoostedW -> Fill("qq->bqq",0);
                             if(otherBintoBoostW.at(0) || otherBintoBoostW.at(1)) h_Hs_otherBclose_BoostedW -> Fill("qq->bqq",1);
                             else h_Hs_otherBclose_BoostedW -> Fill("qq",1);
			     
			     h_Hs_Bc_qq_Top_pT                            -> Fill(Top_fromH_p4.pt());
			     h_Hs_Bc_qq_W_pT                              -> Fill(W_fromTop_fromH_p4.pt());
			   } //if apply cuts, boosted W
		     
			 else // if the ptCut on the quarks hasn't been applied
			   {
			     h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("qq",1);
			     h_Hs_isbQuarkintoBaryCenter           -> Fill("No b",1);
			   }
		       }
		     else if(deltaR_obj1_obj2 < 0.8 && (deltaR_b_fromTop_obj1 < 0.8 || deltaR_b_fromTop_obj2 < 0.8))
		       {
			 if(applyCut)
                           {
			     h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq-bq",1);
			     h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("?",1);
			     h_Hs_Bc_qq_bq_Top_pT                         -> Fill(Top_fromH_p4.pt());
			     h_Hs_Bc_qq_bq_W_pT                           -> Fill(W_fromTop_fromH_p4.pt());
			   }
			 else// if the ptCut on the quarks hasn't been applied
			   {
			     h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("qq-bq",1);
			     h_Hs_isbQuarkintoBaryCenter           -> Fill("?",1);
			   }
		       }
		     else if(deltaR_obj1_obj2 > 0.8 && (deltaR_b_fromTop_obj1 < 0.8 || deltaR_b_fromTop_obj2 < 0.8)) 
		       {
			 if(applyCut)
                           {
			     h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bq",1);
			     h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("b",1);
			     h_Hs_Bc_bq_Top_pT                            -> Fill(Top_fromH_p4.pt());
			     h_Hs_Bc_bq_W_pT                              -> Fill(W_fromTop_fromH_p4.pt());
			     
			     math::XYZTLorentzVector bc_bq_sump4(0,0,0,0);
			     h_Hs_bqcase_Prob_deltaR_lt -> Fill("< 0.4",0);
			     if(deltaR_b_fromTop_obj1 < 0.8) 
			       {
				 bc_bq_sump4 = b_fromTop_fromH_p4 + obj1_fromW_fromTop_p4; 
				 h_Hs_bqcase_deltaR_bq -> Fill(deltaR_b_fromTop_obj1);
				 if(deltaR_b_fromTop_obj1 < 0.4) h_Hs_bqcase_Prob_deltaR_lt -> Fill("< 0.4",1);
				 else                            h_Hs_bqcase_Prob_deltaR_lt -> Fill("> 0.4",1);
			       }
			     else
			       {
				 bc_bq_sump4 = b_fromTop_fromH_p4 + obj2_fromW_fromTop_p4;
				 h_Hs_bqcase_deltaR_bq -> Fill(deltaR_b_fromTop_obj2);
                                 if(deltaR_b_fromTop_obj2 < 0.4) h_Hs_bqcase_Prob_deltaR_lt -> Fill("< 0.4",1);
                                 else                            h_Hs_bqcase_Prob_deltaR_lt -> Fill("> 0.4",1);				 
			       }
			     double MassOf_bq = bc_bq_sump4.M();
			     h_Hs_Bc_MassOf_bq            -> Fill(MassOf_bq);
			   }
			 else// if the ptCut on the quarks hasn't been applied
			   {
			     h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("bq",1);
			     h_Hs_isbQuarkintoBaryCenter           -> Fill("b",1); 
			   }
		       }
		     else // if there is no pair of top decay products with dR < 0.8
		       {
			 if(applyCut)
                           {
			     h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("0",1);
			     h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",1);
			   }
			 else// if the ptCut on the quarks hasn't been applied
			   {
			     h_Hs_QuarksintoBaryCenterMultiplicity -> Fill("0",1);
			     h_Hs_isbQuarkintoBaryCenter           -> Fill("No b",1);
			   }
		       }
		     
		     if(applyCut) h_Hs_OnlyQQ_dR_less_p7_pTcuts -> Fill("True",0);
		     else         h_Hs_OnlyQQ_dR_less_p7        -> Fill("True",0);

		     // OnlyQQ -> Something like BoostedW case
		     if ((deltaR_b_fromTop_obj1 > 0.8) && (deltaR_b_fromTop_obj2 > 0.8) && (deltaR_obj1_obj2 < 0.8))  // %1
		       {
			 if(applyCut) h_Hs_OnlyQQ_dR_less_p7_pTcuts -> Fill("True",1);
			 else         h_Hs_OnlyQQ_dR_less_p7        ->Fill("True",1);

			 double deltaRmin_fatJet_obj1 = 1e6; //give an initial, non sense value             
			 math::XYZTLorentzVector closest_toObj1_fatJet_p4(0,0,0,0);
			 int closest_toObj1_fatJet_Subjets = 0;
			 bool closest_toObj1_fatJet_hasBsubjet = false;
			 double closest_toObj1_fatJet_CSV = 0.0;
			 int AK8jetSD_index = -1;
			 for(AK8JetsSoftDrop fatjet: fEvent.ak8jetsSoftDrop())
			   {
			     AK8jetSD_index++;
			     math::XYZTLorentzVector fatJet_p4(0,0,0,0);
			     fatJet_p4 = fatjet.p4();
			     if((fatJet_p4.pt() < 100.0) || (fatJet_p4.eta() > 2.4) || !fatjet.IDloose()) continue;
			     double deltaR_fatJet_obj1  = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,obj1_fromW_fromTop_p4);
			     if(deltaR_fatJet_obj1 < deltaRmin_fatJet_obj1)
			       {
				 deltaRmin_fatJet_obj1 = deltaR_fatJet_obj1;
				 closest_toObj1_fatJet_p4 = fatjet.p4();
				 closest_toObj1_fatJet_Subjets    = fatjet.nSubjets();
				 closest_toObj1_fatJet_hasBsubjet = fatjet.hasBTagSubjets();
				 closest_toObj1_fatJet_CSV        = fatjet.pfCombinedInclusiveSecondaryVertexV2BJetTags();
			       }
			   }   //for fat-jets loop         
			 
			 if(deltaRmin_fatJet_obj1 <0.8) //%2
			   {
			     double deltaR_fatJet_obj2  = ROOT::Math::VectorUtil::DeltaR(closest_toObj1_fatJet_p4,obj2_fromW_fromTop_p4);
			     h_Hs_Prob_Diquark_match_with_fj ->Fill("Matched",0);
			     if(deltaR_fatJet_obj2 < 0.8) 
			       {
				 if(applyCut)
				   {
				     h_Hs_Prob_Diquark_match_with_fj_pTcuts  -> Fill("Matched",1);
				     
				     h_Hs_MatchedWithDiquark_fj_NumOf_Subjets-> Fill(closest_toObj1_fatJet_Subjets);
				     if(closest_toObj1_fatJet_hasBsubjet) h_Hs_MatchedWithDiquark_fj_hasBsubjet   -> Fill(1);
				     else                                 h_Hs_MatchedWithDiquark_fj_hasBsubjet   -> Fill(0);
				     h_Hs_MatchedWithDiquark_fj_CSV       ->Fill(closest_toObj1_fatJet_CSV);
				     h_Hs_BoostedW_W_pT_Vs_Fatjet_pT      ->Fill(W_fromTop_fromH_p4.pt(),closest_toObj1_fatJet_p4.pt());

				     h_Hs_MatchedWithDiquark_fj_pT_pTcuts    -> Fill(closest_toObj1_fatJet_p4.pt());
				     h_Hs_MatchedWithDiquark_Prob_fj_pT_pTcuts -> Fill("<300 GeV",0);
				     if(closest_toObj1_fatJet_p4.pt() < 300.0) h_Hs_MatchedWithDiquark_Prob_fj_pT_pTcuts-> Fill("<300 GeV",1);
				     else                                      h_Hs_MatchedWithDiquark_Prob_fj_pT_pTcuts-> Fill(">300 GeV",1);
				   }
				 else
				   {
				     h_Hs_Prob_Diquark_match_with_fj -> Fill("Matched",1);
				     h_Hs_MatchedWithDiquark_fj_pT -> Fill(closest_toObj1_fatJet_p4.pt());
				     h_Hs_MatchedWithDiquark_Prob_fj_pT -> Fill("<300 GeV",0);
                                     if(closest_toObj1_fatJet_p4.pt() < 300.0) h_Hs_MatchedWithDiquark_Prob_fj_pT -> Fill("<300 GeV",1);
                                     else                                      h_Hs_MatchedWithDiquark_Prob_fj_pT -> Fill(">300 GeV",1);
				   }
			       }
			     else if(applyCut) h_Hs_Prob_Diquark_match_with_fj_pTcuts ->Fill("Not Matched",1);
			     else              h_Hs_Prob_Diquark_match_with_fj ->Fill("Not Matched",1);
			   } //%2
			 
		       } // %1
		     
		     else if(applyCut) h_Hs_OnlyQQ_dR_less_p7_pTcuts -> Fill("False",1);
		     else              h_Hs_OnlyQQ_dR_less_p7 -> Fill("False",1);
		     
		     
		     // another approach, first the matching with AK8 and after the dR between the objects
		     if (applyCut) // Cut on pT of the quarks
		       {
			 double deltaRmin_fatJet_obj1 = 1e6, deltaRmin_fatJet_obj2 = 1e6, deltaRmin_fatJet_b = 1e6;
			 int MatchedFj_withobj1_index = 0, MatchedFj_withobj2_index = 0, MatchedFj_withB_index = 0;
			 math::XYZTLorentzVector ldgPtFatjet_p4(0,0,0,0);
			 math::XYZTLorentzVector subldgPtFatjet_p4(0,0,0,0);
			 int ldgPtFatjet_index = 0, subldgPtFatjet_index = 0;
			 math::XYZTLorentzVector closest_toObj1_fatJet_p4(0,0,0,0);
			 math::XYZTLorentzVector closest_toB_fatJet_p4(0,0,0,0);
			 bool closest_toObj1_fatJet_hasBsubjet = false, closest_toB_fatJet_hasBsubjet = false;

			 std::vector<double> closest_toObj1_fatJet_Njettiness(5), closest_toB_fatJet_Njettiness(5);
			 for(size_t jettinessInitializer=0; jettinessInitializer < closest_toObj1_fatJet_Njettiness.size(); jettinessInitializer++)
			   {
			     closest_toObj1_fatJet_Njettiness.at(jettinessInitializer) = 0.0;
			     closest_toB_fatJet_Njettiness.at(jettinessInitializer) = 0.0;
			   }


			 int AK8jetSD_index = -1;
			 for(AK8Jet fatjet: fEvent.ak8jets())
			   {
			     AK8jetSD_index++;
			     math::XYZTLorentzVector fatJet_p4(0,0,0,0);
			     fatJet_p4 = fatjet.p4();
			     if((fatJet_p4.pt() < 170.0) || (fatJet_p4.eta() > 2.4) || !fatjet.IDloose()) continue;
			     
			     if(fatJet_p4.pt() > ldgPtFatjet_p4.pt())
			       {
				 subldgPtFatjet_p4    = ldgPtFatjet_p4;
				 subldgPtFatjet_index = ldgPtFatjet_index;   // first the ldg goes to subldg 
				 ldgPtFatjet_p4     = fatJet_p4;          // and then give the new ldg
				 ldgPtFatjet_index  = AK8jetSD_index;
			       }
			     else if(fatJet_p4.pt() > subldgPtFatjet_p4.pt())
			       {
                                 subldgPtFatjet_p4    = fatJet_p4;
                                 subldgPtFatjet_index = AK8jetSD_index;
			       }
			    
			     
			     double deltaR_fatJet_obj1  = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,obj1_fromW_fromTop_p4);
			     double deltaR_fatJet_obj2  = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,obj2_fromW_fromTop_p4);
			     double deltaR_fatJet_b     = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,b_fromTop_fromH_p4);
			     if(deltaR_fatJet_obj1 < deltaRmin_fatJet_obj1)
			       {
				 deltaRmin_fatJet_obj1 = deltaR_fatJet_obj1;
				 closest_toObj1_fatJet_p4 = fatjet.p4();
				 MatchedFj_withobj1_index = AK8jetSD_index;

				 // Njettiness
				 closest_toObj1_fatJet_Njettiness.at(1) = fatjet.NjettinessAK8CHStau1();
				 closest_toObj1_fatJet_Njettiness.at(2) = fatjet.NjettinessAK8CHStau2();
				 closest_toObj1_fatJet_Njettiness.at(3) = fatjet.NjettinessAK8CHStau3();
				 closest_toObj1_fatJet_Njettiness.at(4) = fatjet.NjettinessAK8CHStau4();
			 
			       }
			     if(deltaR_fatJet_obj2 < deltaRmin_fatJet_obj2)
                               {
                                 deltaRmin_fatJet_obj2 = deltaR_fatJet_obj2;
                                 MatchedFj_withobj2_index = AK8jetSD_index;
                               }
			     if(deltaR_fatJet_b < deltaRmin_fatJet_b)
                               {
                                 deltaRmin_fatJet_b = deltaR_fatJet_b;
                                 closest_toB_fatJet_p4 = fatjet.p4();
                                 MatchedFj_withB_index = AK8jetSD_index;
				 //closest_toB_fatJet_hasBsubjet = fatjet.hasBTagSubjets();
				 
				 // Njettiness   
                                 closest_toB_fatJet_Njettiness.at(1) = fatjet.NjettinessAK8CHStau1();
                                 closest_toB_fatJet_Njettiness.at(2) = fatjet.NjettinessAK8CHStau2();
                                 closest_toB_fatJet_Njettiness.at(3) = fatjet.NjettinessAK8CHStau3();
                                 closest_toB_fatJet_Njettiness.at(4) = fatjet.NjettinessAK8CHStau4();

                               }
			   
			     //closest_toObj1_fatJet_hasBsubjet = fatjet.hasBTagSubjets();
			     //closest_toObj1_fatJet_CSV        = fatjet.pfCombinedInclusiveSecondaryVertexV2BJetTags();
			   }   //for fat-jets loop
			 
			 h_Hs_QuarksintoFatJetMultiplicity -> Fill("0",0);  //just to determinate the label of the first bin        
			 h_Hs_QuarksintoFatJetMultiplicity -> Fill("qq",0);
			 h_Hs_QuarksintoFatJetMultiplicity -> Fill("bq",0);
			 
			 h_Hs_BoostedWinFatJet_Prob_deltaR_lt ->Fill("< 0.4",0);

			 if(MatchedFj_withobj1_index == MatchedFj_withobj2_index && MatchedFj_withobj2_index == MatchedFj_withB_index &&
			    deltaRmin_fatJet_obj1 < 0.8 && deltaRmin_fatJet_obj2 < 0.8 && deltaRmin_fatJet_b < 0.8)
			   {
			     h_Hs_QuarksintoFatJetMultiplicity       -> Fill("bqq",1);
			     h_Hs_TopProdInFatJet_fatjet_pT          -> Fill(closest_toB_fatJet_p4.pt());
			     h_Hs_TopProdInFatJet_Top_pT_Vs_fatjet_pT-> Fill(closest_toB_fatJet_p4.pt(), Top_fromH_p4.pt());
			     h_Hs_TopProdInFatJet_Higgs_pT           -> Fill(m.pt());
			     h_Hs_TopProdInFatJet_hasBsubjet         -> Fill(closest_toB_fatJet_hasBsubjet);
			     h_Hs_TopProdInFatJet_Njettinesstau1     -> Fill(closest_toB_fatJet_Njettiness.at(1));
			     h_Hs_TopProdInFatJet_Njettinesstau2     -> Fill(closest_toB_fatJet_Njettiness.at(2));
			     h_Hs_TopProdInFatJet_Njettinesstau3     -> Fill(closest_toB_fatJet_Njettiness.at(3));
                             h_Hs_TopProdInFatJet_Njettinesstau4     -> Fill(closest_toB_fatJet_Njettiness.at(4));
			     double tau2DIVtau1 = (closest_toB_fatJet_Njettiness.at(2))/(closest_toB_fatJet_Njettiness.at(1));
			     double tau3DIVtau2 = (closest_toB_fatJet_Njettiness.at(3))/(closest_toB_fatJet_Njettiness.at(2));
			     h_Hs_TopProdInFatJet_tau2DIVtau1        -> Fill(tau2DIVtau1);
			     h_Hs_TopProdInFatJet_tau3DIVtau2        -> Fill(tau3DIVtau2);
			     h_Hs_TopProdInFatJet_tau2DIVtau1_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau2DIVtau1);
			     h_Hs_TopProdInFatJet_tau3DIVtau2_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau3DIVtau2);
			     
			     bool passtau32cut = false, passtau21cut=false;
			     if(tau2DIVtau1 < 0.6) passtau21cut=true;
			     if(tau3DIVtau2 < 0.67) passtau32cut = true;
			     h_Hs_TopProdInFatJet_TFtau32cut_fatjet_pT -> Fill(passtau32cut, closest_toB_fatJet_p4.pt());
			     h_Hs_TopProdInFatJet_TFtau21cut_fatjet_pT -> Fill(passtau21cut, closest_toB_fatJet_p4.pt());

			     h_Hs_TopProdInFatJet_ldgORsubldg -> Fill("leadingPt",0);
			     h_Hs_TopProdInFatJet_ldgORsubldg -> Fill("subleadingPt",0);
			     if (MatchedFj_withB_index == ldgPtFatjet_index)
			       {
				 h_Hs_TopProdInFatJet_ldgORsubldg -> Fill("leadingPt",1);
				 if(subldgPtFatjet_p4.pt() > 170.0) h_Hs_TopProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT -> Fill(subldgPtFatjet_p4.pt(),ldgPtFatjet_p4.pt()); 
			       }
			     else if (MatchedFj_withB_index == subldgPtFatjet_index)
			       {
				 h_Hs_TopProdInFatJet_ldgORsubldg -> Fill("subleadingPt",1);
				 h_Hs_TopProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT -> Fill(ldgPtFatjet_p4.pt(),subldgPtFatjet_p4.pt());
			       }
			     else h_Hs_TopProdInFatJet_ldgORsubldg -> Fill("NoneOf2",1);
			   }

			 // We could put some "else if" cases into the same "else if" with an "or" operator, but i kept them separately  
			 // in order to be easier to read and understand the different cases
			 else if(MatchedFj_withobj1_index == MatchedFj_withobj2_index && MatchedFj_withobj2_index == MatchedFj_withB_index &&
                                 deltaRmin_fatJet_b > 0.8 && deltaRmin_fatJet_obj1 < 0.8 && deltaRmin_fatJet_obj2 < 0.8)
			   // the case where all 3 are close to the same fatjet but the b is out of the cone of the fatjet (0.8)
                           {
                             h_Hs_QuarksintoFatJetMultiplicity -> Fill("qq",1);
			     h_Hs_BoostedWinFatJet_dR_qq       -> Fill(deltaR_obj1_obj2);
			     if(deltaR_obj1_obj2 < 0.4) h_Hs_BoostedWinFatJet_Prob_deltaR_lt -> Fill("< 0.4",1);
			     else h_Hs_BoostedWinFatJet_Prob_deltaR_lt ->Fill("> 0.4",1);
			     
			     h_Hs_WProdInFatJet_fatjet_pT          -> Fill(closest_toObj1_fatJet_p4.pt());
                             h_Hs_WProdInFatJet_W_pT_Vs_fatjet_pT  -> Fill(closest_toObj1_fatJet_p4.pt(), W_fromTop_fromH_p4.pt());
                             h_Hs_WProdInFatJet_Higgs_pT           -> Fill(m.pt());
                             h_Hs_WProdInFatJet_hasBsubjet         -> Fill(closest_toObj1_fatJet_hasBsubjet);
                             h_Hs_WProdInFatJet_Njettinesstau1     -> Fill(closest_toObj1_fatJet_Njettiness.at(1));
                             h_Hs_WProdInFatJet_Njettinesstau2     -> Fill(closest_toObj1_fatJet_Njettiness.at(2));
                             h_Hs_WProdInFatJet_Njettinesstau3     -> Fill(closest_toObj1_fatJet_Njettiness.at(3));
                             h_Hs_WProdInFatJet_Njettinesstau4     -> Fill(closest_toObj1_fatJet_Njettiness.at(4));
                             double tau2DIVtau1 = (closest_toObj1_fatJet_Njettiness.at(2))/(closest_toObj1_fatJet_Njettiness.at(1));
                             double tau3DIVtau2 = (closest_toObj1_fatJet_Njettiness.at(3))/(closest_toObj1_fatJet_Njettiness.at(2));
                             h_Hs_WProdInFatJet_tau2DIVtau1        -> Fill(tau2DIVtau1);
                             h_Hs_WProdInFatJet_tau3DIVtau2        -> Fill(tau3DIVtau2);
			     h_Hs_WProdInFatJet_tau2DIVtau1_Vs_fatjet_pT -> Fill(closest_toObj1_fatJet_p4.pt(), tau2DIVtau1);
                             h_Hs_WProdInFatJet_tau3DIVtau2_Vs_fatjet_pT -> Fill(closest_toObj1_fatJet_p4.pt(), tau3DIVtau2);

			     bool passtau32cut = false, passtau21cut=false;
                             if(tau2DIVtau1 < 0.6) passtau21cut=true;
                             if(tau3DIVtau2 < 0.67) passtau32cut = true;
                             h_Hs_WProdInFatJet_TFtau32cut_fatjet_pT -> Fill(passtau32cut, closest_toObj1_fatJet_p4.pt());
                             h_Hs_WProdInFatJet_TFtau21cut_fatjet_pT -> Fill(passtau21cut, closest_toObj1_fatJet_p4.pt());

			     h_Hs_WProdInFatJet_ldgORsubldg -> Fill("leadingPt",0);
                             h_Hs_WProdInFatJet_ldgORsubldg -> Fill("subleadingPt",0);
                             if (MatchedFj_withobj1_index == ldgPtFatjet_index)
			       {
				 h_Hs_WProdInFatJet_ldgORsubldg -> Fill("leadingPt",1);
				 if(subldgPtFatjet_p4.pt() > 170.0) h_Hs_WProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT -> Fill(subldgPtFatjet_p4.pt(),ldgPtFatjet_p4.pt());
			       }
                             else if (MatchedFj_withobj1_index == subldgPtFatjet_index)
			       {
				 h_Hs_WProdInFatJet_ldgORsubldg -> Fill("subleadingPt",1);
				 h_Hs_WProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT -> Fill(ldgPtFatjet_p4.pt(),subldgPtFatjet_p4.pt());
			       }
                             else h_Hs_WProdInFatJet_ldgORsubldg -> Fill("NoneOf2",1);

                           }
			 else if(MatchedFj_withobj1_index == MatchedFj_withobj2_index && MatchedFj_withobj2_index != MatchedFj_withB_index &&
				 deltaRmin_fatJet_obj1 < 0.8 && deltaRmin_fatJet_obj2 < 0.8)
			   // the case where only the 2 Wdecay products are close to the same fatjet, into the cone of the fatjet (0.8)
			   {
                             h_Hs_QuarksintoFatJetMultiplicity -> Fill("qq",1);
			     h_Hs_BoostedWinFatJet_dR_qq       -> Fill(deltaR_obj1_obj2);
			     if(deltaR_obj1_obj2 < 0.4) h_Hs_BoostedWinFatJet_Prob_deltaR_lt ->Fill("< 0.4",1);
                             else h_Hs_BoostedWinFatJet_Prob_deltaR_lt ->Fill("> 0.4",1);

			     h_Hs_WProdInFatJet_fatjet_pT          -> Fill(closest_toObj1_fatJet_p4.pt());
                             h_Hs_WProdInFatJet_W_pT_Vs_fatjet_pT  -> Fill(closest_toObj1_fatJet_p4.pt(), W_fromTop_fromH_p4.pt());
                             h_Hs_WProdInFatJet_Higgs_pT           -> Fill(m.pt());
                             h_Hs_WProdInFatJet_hasBsubjet         -> Fill(closest_toObj1_fatJet_hasBsubjet);
                             h_Hs_WProdInFatJet_Njettinesstau1     -> Fill(closest_toObj1_fatJet_Njettiness.at(1));
                             h_Hs_WProdInFatJet_Njettinesstau2     -> Fill(closest_toObj1_fatJet_Njettiness.at(2));
                             h_Hs_WProdInFatJet_Njettinesstau3     -> Fill(closest_toObj1_fatJet_Njettiness.at(3));
                             h_Hs_WProdInFatJet_Njettinesstau4     -> Fill(closest_toObj1_fatJet_Njettiness.at(4));
                             double tau2DIVtau1 = (closest_toObj1_fatJet_Njettiness.at(2))/(closest_toObj1_fatJet_Njettiness.at(1));
                             double tau3DIVtau2 = (closest_toObj1_fatJet_Njettiness.at(3))/(closest_toObj1_fatJet_Njettiness.at(2));
                             h_Hs_WProdInFatJet_tau2DIVtau1        -> Fill(tau2DIVtau1);
                             h_Hs_WProdInFatJet_tau3DIVtau2        -> Fill(tau3DIVtau2);
			     h_Hs_WProdInFatJet_tau2DIVtau1_Vs_fatjet_pT -> Fill(closest_toObj1_fatJet_p4.pt(), tau2DIVtau1);
                             h_Hs_WProdInFatJet_tau3DIVtau2_Vs_fatjet_pT -> Fill(closest_toObj1_fatJet_p4.pt(), tau3DIVtau2);

			     bool passtau32cut = false, passtau21cut=false;
                             if(tau2DIVtau1 < 0.6) passtau21cut=true;
                             if(tau3DIVtau2 < 0.67) passtau32cut = true;
                             h_Hs_WProdInFatJet_TFtau32cut_fatjet_pT -> Fill(passtau32cut, closest_toObj1_fatJet_p4.pt());
                             h_Hs_WProdInFatJet_TFtau21cut_fatjet_pT -> Fill(passtau21cut, closest_toObj1_fatJet_p4.pt());

			     h_Hs_WProdInFatJet_ldgORsubldg -> Fill("leadingPt",0);
                             h_Hs_WProdInFatJet_ldgORsubldg -> Fill("subleadingPt",0);
			     if (MatchedFj_withobj1_index == ldgPtFatjet_index)
                               {
                                 h_Hs_WProdInFatJet_ldgORsubldg -> Fill("leadingPt",1);
                                 if(subldgPtFatjet_p4.pt() > 170.0) h_Hs_WProdInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT -> Fill(subldgPtFatjet_p4.pt(),ldgPtFatjet_p4.pt());
                               }
                             else if (MatchedFj_withobj1_index == subldgPtFatjet_index)
                               {
                                 h_Hs_WProdInFatJet_ldgORsubldg -> Fill("subleadingPt",1);
                                 h_Hs_WProdInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT -> Fill(ldgPtFatjet_p4.pt(),subldgPtFatjet_p4.pt());
                               }
                             else h_Hs_WProdInFatJet_ldgORsubldg -> Fill("NoneOf2",1);

                           }

			 else if(MatchedFj_withobj1_index == MatchedFj_withB_index && MatchedFj_withobj2_index == MatchedFj_withB_index &&
                                 deltaRmin_fatJet_obj2 > 0.8 && deltaRmin_fatJet_obj1 < 0.8 && deltaRmin_fatJet_b < 0.8)
			   // the case where all 3 are close to the same fatjet but obj2 is out of the cone of the fatjet (0.8)
                           {
                             h_Hs_QuarksintoFatJetMultiplicity -> Fill("bq",1);
			     
                             h_Hs_bqInFatJet_fatjet_pT          -> Fill(closest_toB_fatJet_p4.pt());
                             h_Hs_bqInFatJet_Top_pT_Vs_fatjet_pT-> Fill(closest_toB_fatJet_p4.pt(), Top_fromH_p4.pt());
			     h_Hs_bqInFatJet_W_pT_Vs_fatjet_pT  -> Fill(closest_toB_fatJet_p4.pt(), W_fromTop_fromH_p4.pt());
                             h_Hs_bqInFatJet_Higgs_pT           -> Fill(m.pt());
                             h_Hs_bqInFatJet_hasBsubjet         -> Fill(closest_toB_fatJet_hasBsubjet);
                             h_Hs_bqInFatJet_Njettinesstau1     -> Fill(closest_toB_fatJet_Njettiness.at(1));
                             h_Hs_bqInFatJet_Njettinesstau2     -> Fill(closest_toB_fatJet_Njettiness.at(2));
                             h_Hs_bqInFatJet_Njettinesstau3     -> Fill(closest_toB_fatJet_Njettiness.at(3));
                             h_Hs_bqInFatJet_Njettinesstau4     -> Fill(closest_toB_fatJet_Njettiness.at(4));
                             double tau2DIVtau1 = (closest_toB_fatJet_Njettiness.at(2))/(closest_toB_fatJet_Njettiness.at(1));
                             double tau3DIVtau2 = (closest_toB_fatJet_Njettiness.at(3))/(closest_toB_fatJet_Njettiness.at(2));
                             h_Hs_bqInFatJet_tau2DIVtau1        -> Fill(tau2DIVtau1);
                             h_Hs_bqInFatJet_tau3DIVtau2        -> Fill(tau3DIVtau2);
			     h_Hs_bqInFatJet_tau2DIVtau1_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau2DIVtau1);
                             h_Hs_bqInFatJet_tau3DIVtau2_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau3DIVtau2);

			     bool passtau32cut = false, passtau21cut=false;
                             if(tau2DIVtau1 < 0.6) passtau21cut=true;
                             if(tau3DIVtau2 < 0.67) passtau32cut = true;
                             h_Hs_bqInFatJet_TFtau32cut_fatjet_pT -> Fill(passtau32cut, closest_toB_fatJet_p4.pt());
                             h_Hs_bqInFatJet_TFtau21cut_fatjet_pT -> Fill(passtau21cut, closest_toB_fatJet_p4.pt());

			     h_Hs_bqInFatJet_ldgORsubldg -> Fill("leadingPt",0);
                             h_Hs_bqInFatJet_ldgORsubldg -> Fill("subleadingPt",0);
                             if (MatchedFj_withB_index == ldgPtFatjet_index)
			       {
				 h_Hs_bqInFatJet_ldgORsubldg -> Fill("leadingPt",1);
				 if(subldgPtFatjet_p4.pt() > 170.0) h_Hs_bqInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT -> Fill(subldgPtFatjet_p4.pt(),ldgPtFatjet_p4.pt());
			       }
                             else if (MatchedFj_withB_index == subldgPtFatjet_index)
			       {
				 h_Hs_bqInFatJet_ldgORsubldg -> Fill("subleadingPt",1);
				 h_Hs_bqInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT -> Fill(ldgPtFatjet_p4.pt(),subldgPtFatjet_p4.pt());
			       }
                             else h_Hs_bqInFatJet_ldgORsubldg -> Fill("NoneOf2",1);
			     
                           }
                         else if(MatchedFj_withobj1_index == MatchedFj_withB_index && MatchedFj_withobj2_index != MatchedFj_withB_index &&
                                 deltaRmin_fatJet_obj1 < 0.8 && deltaRmin_fatJet_b < 0.8)
			   // the case where only the b and 1 q (obj1) from W are close to the same fatjet, into the cone of the fatjet (0.8)
                           {
                             h_Hs_QuarksintoFatJetMultiplicity -> Fill("bq",1);
			     
			     h_Hs_bqInFatJet_fatjet_pT          -> Fill(closest_toB_fatJet_p4.pt());
                             h_Hs_bqInFatJet_Top_pT_Vs_fatjet_pT-> Fill(closest_toB_fatJet_p4.pt(), Top_fromH_p4.pt());
                             h_Hs_bqInFatJet_W_pT_Vs_fatjet_pT  -> Fill(closest_toB_fatJet_p4.pt(), W_fromTop_fromH_p4.pt());
                             h_Hs_bqInFatJet_Higgs_pT           -> Fill(m.pt());
                             h_Hs_bqInFatJet_hasBsubjet         -> Fill(closest_toB_fatJet_hasBsubjet);
                             h_Hs_bqInFatJet_Njettinesstau1     -> Fill(closest_toB_fatJet_Njettiness.at(1));
                             h_Hs_bqInFatJet_Njettinesstau2     -> Fill(closest_toB_fatJet_Njettiness.at(2));
                             h_Hs_bqInFatJet_Njettinesstau3     -> Fill(closest_toB_fatJet_Njettiness.at(3));
                             h_Hs_bqInFatJet_Njettinesstau4     -> Fill(closest_toB_fatJet_Njettiness.at(4));
                             double tau2DIVtau1 = (closest_toB_fatJet_Njettiness.at(2))/(closest_toB_fatJet_Njettiness.at(1));
                             double tau3DIVtau2 = (closest_toB_fatJet_Njettiness.at(3))/(closest_toB_fatJet_Njettiness.at(2));
                             h_Hs_bqInFatJet_tau2DIVtau1        -> Fill(tau2DIVtau1);
                             h_Hs_bqInFatJet_tau3DIVtau2        -> Fill(tau3DIVtau2);
			     h_Hs_bqInFatJet_tau2DIVtau1_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau2DIVtau1);
                             h_Hs_bqInFatJet_tau3DIVtau2_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau3DIVtau2);

			     bool passtau32cut = false, passtau21cut=false;
                             if(tau2DIVtau1 < 0.6) passtau21cut=true;
                             if(tau3DIVtau2 < 0.67) passtau32cut = true;
                             h_Hs_bqInFatJet_TFtau32cut_fatjet_pT -> Fill(passtau32cut, closest_toB_fatJet_p4.pt());
                             h_Hs_bqInFatJet_TFtau21cut_fatjet_pT -> Fill(passtau21cut, closest_toB_fatJet_p4.pt());

			     h_Hs_bqInFatJet_ldgORsubldg -> Fill("leadingPt",0);
                             h_Hs_bqInFatJet_ldgORsubldg -> Fill("subleadingPt",0);
			     if (MatchedFj_withB_index == ldgPtFatjet_index)
                               {
                                 h_Hs_bqInFatJet_ldgORsubldg -> Fill("leadingPt",1);
                                 if(subldgPtFatjet_p4.pt() > 170.0) h_Hs_bqInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT -> Fill(subldgPtFatjet_p4.pt(),ldgPtFatjet_p4.pt());
                               }
                             else if (MatchedFj_withB_index == subldgPtFatjet_index)
                               {
                                 h_Hs_bqInFatJet_ldgORsubldg -> Fill("subleadingPt",1);
                                 h_Hs_bqInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT -> Fill(ldgPtFatjet_p4.pt(),subldgPtFatjet_p4.pt());
                               }
                             else h_Hs_bqInFatJet_ldgORsubldg -> Fill("NoneOf2",1);

                           }
			 else if(MatchedFj_withobj2_index == MatchedFj_withB_index && MatchedFj_withobj1_index == MatchedFj_withB_index &&
                                 deltaRmin_fatJet_obj1 > 0.8 && deltaRmin_fatJet_obj2 < 0.8 && deltaRmin_fatJet_b < 0.8)
			   // the case where all 3 are close to the same fatjet but obj1 from W is out of the cone of the fatjet (0.8)
                           {
                             h_Hs_QuarksintoFatJetMultiplicity -> Fill("bq",1);
			     
			     h_Hs_bqInFatJet_fatjet_pT          -> Fill(closest_toB_fatJet_p4.pt());
                             h_Hs_bqInFatJet_Top_pT_Vs_fatjet_pT-> Fill(closest_toB_fatJet_p4.pt(), Top_fromH_p4.pt());
                             h_Hs_bqInFatJet_W_pT_Vs_fatjet_pT  -> Fill(closest_toB_fatJet_p4.pt(), W_fromTop_fromH_p4.pt());
                             h_Hs_bqInFatJet_Higgs_pT           -> Fill(m.pt());
                             h_Hs_bqInFatJet_hasBsubjet         -> Fill(closest_toB_fatJet_hasBsubjet);
                             h_Hs_bqInFatJet_Njettinesstau1     -> Fill(closest_toB_fatJet_Njettiness.at(1));
                             h_Hs_bqInFatJet_Njettinesstau2     -> Fill(closest_toB_fatJet_Njettiness.at(2));
                             h_Hs_bqInFatJet_Njettinesstau3     -> Fill(closest_toB_fatJet_Njettiness.at(3));
                             h_Hs_bqInFatJet_Njettinesstau4     -> Fill(closest_toB_fatJet_Njettiness.at(4));
                             double tau2DIVtau1 = (closest_toB_fatJet_Njettiness.at(2))/(closest_toB_fatJet_Njettiness.at(1));
                             double tau3DIVtau2 = (closest_toB_fatJet_Njettiness.at(3))/(closest_toB_fatJet_Njettiness.at(2));
                             h_Hs_bqInFatJet_tau2DIVtau1        -> Fill(tau2DIVtau1);
                             h_Hs_bqInFatJet_tau3DIVtau2        -> Fill(tau3DIVtau2);
			     h_Hs_bqInFatJet_tau2DIVtau1_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau2DIVtau1);
                             h_Hs_bqInFatJet_tau3DIVtau2_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau3DIVtau2);
			     
			     bool passtau32cut = false, passtau21cut=false;
                             if(tau2DIVtau1 < 0.6) passtau21cut=true;
                             if(tau3DIVtau2 < 0.67) passtau32cut = true;
                             h_Hs_bqInFatJet_TFtau32cut_fatjet_pT -> Fill(passtau32cut, closest_toB_fatJet_p4.pt());
                             h_Hs_bqInFatJet_TFtau21cut_fatjet_pT -> Fill(passtau21cut, closest_toB_fatJet_p4.pt());

			     h_Hs_bqInFatJet_ldgORsubldg -> Fill("leadingPt",0);
                             h_Hs_bqInFatJet_ldgORsubldg -> Fill("subleadingPt",0);
			     if (MatchedFj_withB_index == ldgPtFatjet_index)
                               {
                                 h_Hs_bqInFatJet_ldgORsubldg -> Fill("leadingPt",1);
                                 if(subldgPtFatjet_p4.pt() > 170.0) h_Hs_bqInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT -> Fill(subldgPtFatjet_p4.pt(),ldgPtFatjet_p4.pt());
                               }
                             else if (MatchedFj_withB_index == subldgPtFatjet_index)
                               {
                                 h_Hs_bqInFatJet_ldgORsubldg -> Fill("subleadingPt",1);
                                 h_Hs_bqInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT -> Fill(ldgPtFatjet_p4.pt(),subldgPtFatjet_p4.pt());
                               }
                             else h_Hs_bqInFatJet_ldgORsubldg -> Fill("NoneOf2",1);
			     
                           }
                         else if(MatchedFj_withobj2_index == MatchedFj_withB_index && MatchedFj_withobj1_index != MatchedFj_withB_index &&
                                 deltaRmin_fatJet_obj2 < 0.8 && deltaRmin_fatJet_b < 0.8)
			   // the case where only the b and 1 q (obj2) from W are close to the same fatjet, into the cone of the fatjet (0.8)
                           {
                             h_Hs_QuarksintoFatJetMultiplicity -> Fill("bq",1);
			     
			     h_Hs_bqInFatJet_fatjet_pT          -> Fill(closest_toB_fatJet_p4.pt());
                             h_Hs_bqInFatJet_Top_pT_Vs_fatjet_pT-> Fill(closest_toB_fatJet_p4.pt(), Top_fromH_p4.pt());
                             h_Hs_bqInFatJet_W_pT_Vs_fatjet_pT  -> Fill(closest_toB_fatJet_p4.pt(), W_fromTop_fromH_p4.pt());
                             h_Hs_bqInFatJet_Higgs_pT           -> Fill(m.pt());
                             h_Hs_bqInFatJet_hasBsubjet         -> Fill(closest_toB_fatJet_hasBsubjet);
                             h_Hs_bqInFatJet_Njettinesstau1     -> Fill(closest_toB_fatJet_Njettiness.at(1));
                             h_Hs_bqInFatJet_Njettinesstau2     -> Fill(closest_toB_fatJet_Njettiness.at(2));
                             h_Hs_bqInFatJet_Njettinesstau3     -> Fill(closest_toB_fatJet_Njettiness.at(3));
                             h_Hs_bqInFatJet_Njettinesstau4     -> Fill(closest_toB_fatJet_Njettiness.at(4));
                             double tau2DIVtau1 = (closest_toB_fatJet_Njettiness.at(2))/(closest_toB_fatJet_Njettiness.at(1));
                             double tau3DIVtau2 = (closest_toB_fatJet_Njettiness.at(3))/(closest_toB_fatJet_Njettiness.at(2));
                             h_Hs_bqInFatJet_tau2DIVtau1        -> Fill(tau2DIVtau1);
                             h_Hs_bqInFatJet_tau3DIVtau2        -> Fill(tau3DIVtau2);
			     h_Hs_bqInFatJet_tau2DIVtau1_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau2DIVtau1);
                             h_Hs_bqInFatJet_tau3DIVtau2_Vs_fatjet_pT -> Fill(closest_toB_fatJet_p4.pt(), tau3DIVtau2);

			     bool passtau32cut = false, passtau21cut=false;
                             if(tau2DIVtau1 < 0.6) passtau21cut=true;
                             if(tau3DIVtau2 < 0.67) passtau32cut = true;
                             h_Hs_bqInFatJet_TFtau32cut_fatjet_pT -> Fill(passtau32cut, closest_toB_fatJet_p4.pt());
                             h_Hs_bqInFatJet_TFtau21cut_fatjet_pT -> Fill(passtau21cut, closest_toB_fatJet_p4.pt());

			     h_Hs_bqInFatJet_ldgORsubldg -> Fill("leadingPt",0);
                             h_Hs_bqInFatJet_ldgORsubldg -> Fill("subleadingPt",0);
                             if (MatchedFj_withB_index == ldgPtFatjet_index)
                               {
                                 h_Hs_bqInFatJet_ldgORsubldg -> Fill("leadingPt",1);
                                 if(subldgPtFatjet_p4.pt() > 170.0) h_Hs_bqInFatJet_matchedFatJetIsldg_pT_Vs_subldg_pT -> Fill(subldgPtFatjet_p4.pt(),ldgPtFatjet_p4.pt());
                               }
                             else if (MatchedFj_withB_index == subldgPtFatjet_index)
                               {
                                 h_Hs_bqInFatJet_ldgORsubldg -> Fill("subleadingPt",1);
                                 h_Hs_bqInFatJet_matchedFatJetIssubldg_pT_Vs_ldg_pT -> Fill(ldgPtFatjet_p4.pt(),subldgPtFatjet_p4.pt());
                               }
                             else h_Hs_bqInFatJet_ldgORsubldg -> Fill("NoneOf2",1);

                           }
			 else h_Hs_QuarksintoFatJetMultiplicity -> Fill("0",1);

			 
		       }//if apply cuts
		     
		   }// for applying cuts
		 
		 
	       }   // if top from Higgs
	     
	     else //if top is not from Higgs//////////////////////////////////////////////////////////////////////////////////////////////
	       {
		 math::XYZTLorentzVector b_fromTop_NOTfromH_p4(0,0,0,0), obj1_fromW_fromTop_p4(0,0,0,0), obj2_fromW_fromTop_p4(0,0,0,0);
		 math::XYZTLorentzVector W_fromTop_NOTfromH_p4(0,0,0,0), Top_NOTfromH_p4(0,0,0,0);
                 bool isTheRightTop = true;
                 bool obj1_exist = false;
                 bool bQuarkIntoFatJet = false;

                 Top_NOTfromH_p4 = p.p4();
		 unsigned int lastcopy = GetTheLastCopy(p.index()); // we need to find the last copy in order to find the real daughters      
                 genParticle lastT; //create the particle object of the last copy of top                                                      
                 lastT =  fEvent.genparticles().getGenParticles()[lastcopy];
		 std::vector<short> lastT_daughters = lastT.daughters();

                 for (unsigned int i = 0; i < lastT_daughters.size() ; i++)
                   {
                     genParticle d; //create the particle object   
                     d =  fEvent.genparticles().getGenParticles()[lastT_daughters.at(i)];
                     if (std::abs(d.pdgId()) == 24) // if W from top             
                       {
                         W_fromTop_NOTfromH_p4 = d.p4();
                         genParticle lastW;
                         unsigned int lastWcopy = GetTheLastCopy(d.index()); // we need to find the last copy in order to find the real daughters                                                                                                                                           
			 lastW = fEvent.genparticles().getGenParticles()[lastWcopy];
			 std::vector<short> W_daughters = lastW.daughters();
			 for (size_t k = 0; k < W_daughters.size() ; k++)
			   {
			     genParticle top_grand_d; //create the particle object, i define this as grand_d, the daus of the W, so the grand-daughters of Top
			     top_grand_d =  fEvent.genparticles().getGenParticles()[W_daughters.at(k)];
			     if ( !mcTools.IsQuark(top_grand_d.pdgId()) ) //if the daus of W are not quarks (f.e. leptons), we reject the top 
                               {
                                 isTheRightTop = false;
                                 continue;
                               }

                             if (obj1_exist == false)
                               {
                                 obj1_exist = true;
                                 obj1_fromW_fromTop_p4 = top_grand_d.p4();
                               }
                             else
                               {
                                 obj2_fromW_fromTop_p4 = top_grand_d.p4();
                               }
                           }      //for loop, W's daus             
                         if (!isTheRightTop) continue; // move out from the checking loop, reject this top                                    
                       } // if W from top
		     
		     if (std::abs(d.pdgId()) == 5) // if b from top                                                                           
                       {
                         b_fromTop_NOTfromH_p4 = d.p4();
                       }
                   }//for top NOT from Higgs daus

                 if (!isTheRightTop) continue; // move out from the checking loop, reject this top

                 bool TheTopQuarkIsMatched_WithFatJet = false;
                 bool TheWQuarkIsMatched_WithFatJet = false;

                 double deltaRmin_fatJet_topQuark = 1e6; //give an initial, non sense value    
                 double deltaRmin_fatJet_WQuark = 1e6; //give an initial, non sense value                
		 math::XYZTLorentzVector closest_toTop_fatJet_p4(0,0,0,0), closest_toW_fatJet_p4(0,0,0,0);
                 int AK8jetSD_index = -1;
		 for(AK8JetsSoftDrop fatjet: fEvent.ak8jetsSoftDrop())
                   {
                     AK8jetSD_index++;
		     math::XYZTLorentzVector fatJet_p4(0,0,0,0);
                     fatJet_p4 = fatjet.p4();
                     double deltaR_fatJet_topQuark  = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,Top_NOTfromH_p4);
                     double deltaR_fatJet_WQuark    = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,W_fromTop_NOTfromH_p4);
                     if(deltaR_fatJet_topQuark < deltaRmin_fatJet_topQuark)
                       {
                         deltaRmin_fatJet_topQuark = deltaR_fatJet_topQuark;
                         closest_toTop_fatJet_p4 = fatjet.p4();
                       }
                     if(deltaR_fatJet_WQuark < deltaRmin_fatJet_WQuark)
                       {
                         deltaRmin_fatJet_WQuark = deltaR_fatJet_WQuark;
                         closest_toW_fatJet_p4 = fatjet.p4();
		       }
                   }   //for fat-jets loop  

                 if(deltaRmin_fatJet_topQuark <0.8) TheTopQuarkIsMatched_WithFatJet = true;
                 if(deltaRmin_fatJet_WQuark <0.8)   TheWQuarkIsMatched_WithFatJet = true;
                 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		 double deltaR_b_fromTop_obj1 = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, obj1_fromW_fromTop_p4);
                 double deltaR_b_fromTop_obj2 = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, obj2_fromW_fromTop_p4);
                 double deltaR_obj1_obj2      = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, obj2_fromW_fromTop_p4);
                 double deltaRmax = -1; // give an initial, non sense value
                 double deltaRmin = 1e6; // give an initial, non sense value       

                 // find the dRmin                   
                 if(deltaR_b_fromTop_obj1 > deltaR_b_fromTop_obj2 && deltaR_b_fromTop_obj1 > deltaR_obj1_obj2)
                   {
                     deltaRmax = deltaR_b_fromTop_obj1;
                   }
                 else if(deltaR_b_fromTop_obj2 > deltaR_b_fromTop_obj1 && deltaR_b_fromTop_obj2 > deltaR_obj1_obj2)
                   {
                     deltaRmax = deltaR_b_fromTop_obj2;
                   }
                 else deltaRmax = deltaR_obj1_obj2;
		 
		 h_objectsNOTfromHiggstop_maxdR             -> Fill(deltaRmax);
		 
		 //find the dRmin                                  
                 if(deltaR_b_fromTop_obj1 < deltaR_b_fromTop_obj2 && deltaR_b_fromTop_obj1 < deltaR_obj1_obj2)   
		   {
		     deltaRmin = deltaR_b_fromTop_obj1;
		   }
                 else if(deltaR_b_fromTop_obj2 < deltaR_b_fromTop_obj1 && deltaR_b_fromTop_obj2 < deltaR_obj1_obj2) 
		   {
		     deltaRmin = deltaR_b_fromTop_obj2;
		   }
                 else deltaRmin = deltaR_obj1_obj2;

                 h_objectsNOTfromHiggstop_mindR                                -> Fill(deltaRmin);
                 h_objectsNOTfromHiggstop_Prob_mindR_lt_p8                     -> Fill("<0.8",0);

                 if(deltaRmin < 0.8)
                   {
                     h_NotHs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj -> Fill("Matched",0);
                     if(TheTopQuarkIsMatched_WithFatJet)
                       {
                         h_NotHs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj -> Fill("Matched",1);
                         h_NotHs_objectsfromtop_mindR_ltp8_matchedfj_pt -> Fill(closest_toTop_fatJet_p4.pt());
                       }
                     else h_NotHs_objectsfromtop_mindR_ltp8_Prob_top_match_with_fj -> Fill("Not Matched",1);

                     h_objectsNOTfromHiggstop_Prob_mindR_lt_p8 -> Fill("<0.8",1);
                     h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("< 400",0);
                     h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("400<pT<500",0);
                     if     (Top_NOTfromH_p4.pt() > 500.0) 
		       {
			 h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("> 500",1);
			 if(TheTopQuarkIsMatched_WithFatJet)
			   {
			     h_NotHs_objectsfromtop_mindR_ltp8__topPt_more_500__matchedfj_pt -> Fill(closest_toTop_fatJet_p4.pt());
			   }
		       }
                     else if(Top_NOTfromH_p4.pt() > 400.0)
		       {
			 h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("400<pT<500",1);
		       }
                     else h_objectsNOTfromHiggstop_Prob_mindR_lt_p8_and_TopPt_more_than -> Fill("< 400",1);
		     
		     if       (Top_NOTfromH_p4.pt() > 400.0)
		       {
			 if(TheTopQuarkIsMatched_WithFatJet)
                           {
                             h_NotHs_objectsfromtop_mindR_ltp8__topPt_more_400__matchedfj_pt -> Fill(closest_toTop_fatJet_p4.pt());
                           }
		       }
                   }

                 else  h_objectsNOTfromHiggstop_Prob_mindR_lt_p8               -> Fill(">0.8",1);
		 
		 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                 /////////////////////////////////////////////////BaryCenter//////////////////////////////////////////////////////////////////
		 /////////////assuming that is the top-direction ////////////////////////////////////////////////////////////////////////////
                 double deltaR_b_fromTop_Top    = ROOT::Math::VectorUtil::DeltaR(b_fromTop_NOTfromH_p4, Top_NOTfromH_p4);
                 double deltaR_obj1_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(obj1_fromW_fromTop_p4, Top_NOTfromH_p4);
                 double deltaR_obj2_fromTop_Top = ROOT::Math::VectorUtil::DeltaR(obj2_fromW_fromTop_p4, Top_NOTfromH_p4);
                 unsigned int decayproductsintofatjet = 0;
                 double objtop_deltaRmax = -1;  // give an initial, non sense value    
                 double objtop_deltaRmin = 999;  // give an initial, non sense value    

                 h_NotHs_which_objectfromtop_maxdR -> Fill ("W'obj",0) ; //just to determinate the first label   

                 if(deltaR_b_fromTop_Top > deltaR_obj1_fromTop_Top && deltaR_b_fromTop_Top > deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmax = deltaR_b_fromTop_Top;
                     h_NotHs_which_objectfromtop_maxdR -> Fill("b",1) ;
		     
		     if(deltaR_b_fromTop_Top > 0.8)
                       {
                         h_NotHs_mostdistantfromtop_isb__top_pT      -> Fill(Top_NOTfromH_p4.pt());
                         h_NotHs_mostdistantfromtop_isb__W_pT        -> Fill(W_fromTop_NOTfromH_p4.pt());
                         h_NotHs_mostdistantfromtop_isb__b_pT        -> Fill(b_fromTop_NOTfromH_p4.pt());
                         h_NotHs_mostdistantfromtop_isb__dRmax_b_top -> Fill(deltaR_b_fromTop_Top);

                         if(deltaR_b_fromTop_obj1 < deltaR_b_fromTop_obj2) 
			   {
			     h_NotHs_mostdistantfromtop_isb_dRmin_b_objfromW -> Fill(deltaR_b_fromTop_obj1);
			   }
                         else h_NotHs_mostdistantfromtop_isb_dRmin_b_objfromW -> Fill(deltaR_b_fromTop_obj2);
                       }

		     if(deltaR_b_fromTop_Top > 0.8 && deltaR_obj1_obj2 < 0.7)
                       {
                         h_NotHs_mostdistantfromtop_isb__dRqq_Vs_W_pT   -> Fill(deltaR_obj1_obj2,W_fromTop_NOTfromH_p4.pt());
                         h_NotHs_mostdistantfromtop_isb__dRqq_Vs_top_pT -> Fill(deltaR_obj1_obj2,Top_NOTfromH_p4.pt());

                         h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj -> Fill("Matched",0);
                         if(TheTopQuarkIsMatched_WithFatJet)
                           {
                             h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj -> Fill("Matched",1);
                             h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Topmatchedfj_pt -> Fill(closest_toTop_fatJet_p4.pt());
                           }
                         else h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_top_match_with_fj -> Fill("Not Matched",1);

                         h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj -> Fill("Matched",0);
                         if(TheWQuarkIsMatched_WithFatJet)
                           {
                             h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj -> Fill("Matched",1);
                             h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Wmatchedfj_pt -> Fill(closest_toW_fatJet_p4.pt());
                           }
                         else h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_Prob_W_match_with_fj -> Fill("Not Matched",1);

                         h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("< 400",0);
                         h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("400<pT<500",0);
                         h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("< 200",0);
                         h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("200<pT<300",0);

                         if     (Top_NOTfromH_p4.pt() > 500.0) 
			   {
			     h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("> 500",1);
			     if(TheTopQuarkIsMatched_WithFatJet)
			       {
				 h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more500_matchedfj_pt->Fill(closest_toTop_fatJet_p4.pt());
			       }
			   }
                         else if(Top_NOTfromH_p4.pt() > 400.0)
			   {
			     h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("400<pT<500",1);
			   }
                         else h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_TopPt_more_than -> Fill("< 400",1);

			 if(Top_NOTfromH_p4.pt() > 400.0)
			   {
			     h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__TopPt_more400_matchedfj_pt->Fill(closest_toTop_fatJet_p4.pt());
			   }


                         if     (W_fromTop_NOTfromH_p4.pt() > 300.0)
			   {
			     h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill(">300",1);
			     if(TheWQuarkIsMatched_WithFatJet)
                               {
                                 h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more300_matchedfj_pt->Fill(closest_toW_fatJet_p4.pt());
                               }
			   }
                         else if(W_fromTop_NOTfromH_p4.pt() > 200.0) 
			   {
			     h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("200<pT<300",1);
			   }
                         else h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__and_WPt_more_than -> Fill("< 200",1);
			 
			 if (W_fromTop_NOTfromH_p4.pt() > 200.0)
			   {
			     if(TheWQuarkIsMatched_WithFatJet)
			       {
				 h_NotHs_mostdistantfromtop_isb__dRqq_less_p7__WPt_more200_matchedfj_pt->Fill(closest_toW_fatJet_p4.pt());
			       }
			   }
                       }
		     
		     h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("True",0);
                     if(deltaR_b_fromTop_Top > 0.8 && deltaR_obj1_obj2 < 0.7) 
		       {
			 h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("True",1);
		       }
                     else h_NotHs_mostdistantfromtop_isb__dRqq_less_p7_tf -> Fill("False",1);
		   }

                 else if(deltaR_obj1_fromTop_Top > deltaR_b_fromTop_Top && deltaR_obj1_fromTop_Top > deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmax = deltaR_obj1_fromTop_Top;
                     h_NotHs_which_objectfromtop_maxdR-> Fill("W'obj",1);
                   }
                 else
                   {
                     objtop_deltaRmax = deltaR_obj2_fromTop_Top;
                     h_NotHs_which_objectfromtop_maxdR-> Fill("W'obj",1);
                   }
                 h_NotHs_objectsfromtop_top_maxdR -> Fill(objtop_deltaRmax);

                 if(deltaR_b_fromTop_Top < deltaR_obj1_fromTop_Top && deltaR_b_fromTop_Top < deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmin = deltaR_b_fromTop_Top;
                   }
                 else if(deltaR_obj1_fromTop_Top < deltaR_b_fromTop_Top && deltaR_obj1_fromTop_Top < deltaR_obj2_fromTop_Top)
                   {
                     objtop_deltaRmin = deltaR_obj1_fromTop_Top;
                   }
                 else objtop_deltaRmin = deltaR_obj2_fromTop_Top;
                 h_NotHs_objectsfromtop_top_mindR -> Fill(objtop_deltaRmin);


                 if(deltaR_b_fromTop_Top < 0.8)  // b into fatjet             
                   {
                     decayproductsintofatjet++;
                     bQuarkIntoFatJet = true;
                   }
                 else bQuarkIntoFatJet = false;
		 
		 if(deltaR_obj1_fromTop_Top < 0.8) decayproductsintofatjet++;  // obj1 into fatjet        
                 if(deltaR_obj2_fromTop_Top < 0.8) decayproductsintofatjet++;  // obj2 into fatjet        

                 h_Probdecayproductsintofatjet_NOTHs          -> Fill("0",0);  //just to determinate the label of the first bin   
                 h_Probdecayproductsintofatjet_NOTHs          -> Fill("1",0);  //just to determinate the label of the second bin 
                 h_Probdecayproductsintofatjet_NOTHs          -> Fill("2",0);  //just to determinate the label of the third bin   

                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("0",0);  //just to determinate the label of the first bin      
                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("q",0);
                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("b",0);
                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("qq",0);
                 h_NotHs_QuarksintofatjetMultiplicity -> Fill("bq",0);

                 h_NotHs_isbQuarkintofatjet           -> Fill("No b",0);  //just to determinate the label of the first bin
		 
		 if      (decayproductsintofatjet == 0)
                   {
                     h_Probdecayproductsintofatjet_NOTHs  -> Fill("0",1);
                     h_NotHs_QuarksintofatjetMultiplicity -> Fill("0",1);
                     h_NotHs_isbQuarkintofatjet           -> Fill("No b",1);
                   }
                 else if (decayproductsintofatjet == 1)
                   {
                     h_Probdecayproductsintofatjet_NOTHs -> Fill("1",1);
                     if(bQuarkIntoFatJet)
                       {
                         h_NotHs_QuarksintofatjetMultiplicity -> Fill("b",1);
                         h_NotHs_isbQuarkintofatjet           -> Fill("b",1);
                       }
                     else
                       {
                         h_NotHs_QuarksintofatjetMultiplicity -> Fill("q",1);
                         h_NotHs_isbQuarkintofatjet           -> Fill("No b",1);
                       }
                   }
		 else if (decayproductsintofatjet == 2)
                   {
                     h_Probdecayproductsintofatjet_NOTHs -> Fill("2",1);
                     if(bQuarkIntoFatJet)
                       {
                         h_NotHs_QuarksintofatjetMultiplicity -> Fill("bq",1);
                         h_NotHs_isbQuarkintofatjet           -> Fill("b",1);
                       }
                     else
                       {
                         h_NotHs_QuarksintofatjetMultiplicity -> Fill("qq",1);
                         h_NotHs_isbQuarkintofatjet           -> Fill("No b",1);
                       }
                   }
                 else
                   {
                     h_Probdecayproductsintofatjet_NOTHs  -> Fill("All 3",1);
                     h_NotHs_QuarksintofatjetMultiplicity -> Fill("bqq",1);
                     h_NotHs_isbQuarkintofatjet           -> Fill("b",1);
                   }
		 
                 if (decayproductsintofatjet == 3)
                   {
                     unsigned int Num_pairs_withenough_deltaR_amongthemselves = 0;
                     if      (deltaR_b_fromTop_obj1 < 0.4 && deltaR_b_fromTop_obj2 < 0.4 && deltaR_obj1_obj2 < 0.4)
                       {
                         Num_pairs_withenough_deltaR_amongthemselves = 0;
                       }

                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 < 0.4 && deltaR_b_fromTop_obj2 < 0.4)
                       {
                         Num_pairs_withenough_deltaR_amongthemselves = 1;
                       }
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 < 0.4 && deltaR_b_fromTop_obj1 < 0.4)
                       {
                         Num_pairs_withenough_deltaR_amongthemselves = 1;
                       }
                     else if (deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj1 < 0.4 && deltaR_b_fromTop_obj2 < 0.4)
                       {
                         Num_pairs_withenough_deltaR_amongthemselves = 1;
                       }

                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj2 < 0.4)
                       {
                         Num_pairs_withenough_deltaR_amongthemselves = 2;
                       }
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 > 0.4 && deltaR_b_fromTop_obj1 < 0.4)
                       {
                         Num_pairs_withenough_deltaR_amongthemselves = 2;
                       }
                     else if (deltaR_b_fromTop_obj2 > 0.4 && deltaR_b_fromTop_obj1 > 0.4 && deltaR_obj1_obj2 <0.4 )
                       {
                         Num_pairs_withenough_deltaR_amongthemselves = 2;
                       }

                     else if (deltaR_b_fromTop_obj1 > 0.4 && deltaR_b_fromTop_obj2 > 0.4 && deltaR_obj1_obj2 > 0.4)
                       {
                         Num_pairs_withenough_deltaR_amongthemselves = 3;
                       }
		     
		     h_Pairsinbarycenter_enoughdeltaR_NOTHs         -> Fill("0",0);  //just to determinate the label of the first bin
                     h_Pairsinbarycenter_enoughdeltaR_NOTHs         -> Fill("1 pair",0);  //just to determinate the label of the second bin
                     h_Pairsinbarycenter_enoughdeltaR_NOTHs         -> Fill("2 pairs",0);  //just to determinate the label of the third bin 

                     if      (Num_pairs_withenough_deltaR_amongthemselves == 0)   h_Pairsinbarycenter_enoughdeltaR_NOTHs -> Fill("0",1);
                     else if (Num_pairs_withenough_deltaR_amongthemselves == 1)   h_Pairsinbarycenter_enoughdeltaR_NOTHs -> Fill("1 pair",1);
                     else if (Num_pairs_withenough_deltaR_amongthemselves == 2)   h_Pairsinbarycenter_enoughdeltaR_NOTHs -> Fill("2 pairs",1);
                     else                                                         h_Pairsinbarycenter_enoughdeltaR_NOTHs -> Fill("All 3",1);

		     h_Iffatjet_NotHs_Top_pT         -> Fill (Top_NOTfromH_p4.pt());
		     h_Iffatjet_NotHs_EventsWithHighTop_pT -> Fill("top.pT < 400 GeV", 0);
		     if(Top_NOTfromH_p4.pt() > 400.0)       h_Iffatjet_NotHs_EventsWithHighTop_pT -> Fill("top.pT > 400 GeV", 1);
                     else                                   h_Iffatjet_NotHs_EventsWithHighTop_pT -> Fill("top.pT < 400 GeV", 1);
		     
		   } // All 3 decay products into 0.8
		 
		 /////////////////////////////////////////////////BaryCenter//////////////////////////////////////////////////////////////////
		 // barycenter not the top-direction.  Actually here take place a different approach.  If at least 2 objects have dr<0.8,
                 // had been assumed that they are into the barycenter.  After that, comes the matching with AK8 jets/////////////////////////
		 h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("0",0);  //just to determinate the label of the first bin
                 h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("qq",0);
                 h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("bq",0);
                 h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("qq-bq",0);

                 h_NotHs_isbQuarkintoBaryCenter           -> Fill("No b",0);  //just to determinate the label of the first bin
                 h_NotHs_isbQuarkintoBaryCenter           -> Fill("?",0);

                 h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("0",0);  //just to determinate the label of the first bin 
                 h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq",0);
                 h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bq",0);
                 h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq-bq",0);

                 h_NotHs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",0);  //just to determinate the label of the first bin  
                 h_NotHs_isbQuarkintoBaryCenter_pTcuts           -> Fill("?",0);
		 
		 for (unsigned int applyCut = 0; applyCut < 2 ; applyCut++)
                   {
                     double quarksPtCut = 30.0;
                     if(applyCut)
                       {
                         h_NotHs_QuarksFromTop_Passed_pTcuts -> Fill("Passed",0);
                         if((obj1_fromW_fromTop_p4.pt() < quarksPtCut) || (obj2_fromW_fromTop_p4.pt() < quarksPtCut) ||  (b_fromTop_NOTfromH_p4.pt() < quarksPtCut))
                           {
                             h_NotHs_QuarksFromTop_Passed_pTcuts -> Fill("Not Passed",1);
                             continue;
                           }
                         else h_NotHs_QuarksFromTop_Passed_pTcuts -> Fill("Passed",1);
                       }

                     if(applyCut)
                       {
			 h_NotHs_objectsfromtop_dRmax_pTcuts        -> Fill(deltaRmax);
                         h_NotHs_objectsfromtop_Prob_dRmax_pTcuts        -> Fill("< 0.8",0);
                         if(deltaRmax < 0.8) h_NotHs_objectsfromtop_Prob_dRmax_pTcuts -> Fill("< 0.8",1);
                         else          h_NotHs_objectsfromtop_Prob_dRmax_pTcuts -> Fill("> 0.8",1);

                         h_NotHs_objectsfromtop_dRmin_pTcuts        -> Fill(deltaRmin);
                         h_NotHs_objectsfromtop_Prob_dRmin_pTcuts        -> Fill("< 0.8",0);
                         if(deltaRmin < 0.8) h_NotHs_objectsfromtop_Prob_dRmin_pTcuts -> Fill("< 0.8",1);
                         else          h_NotHs_objectsfromtop_Prob_dRmin_pTcuts -> Fill("> 0.8",1);
			 
                         h_NotHs_QuarksFromW_deltaR_pTcuts          -> Fill(deltaR_obj1_obj2);
                         h_NotHs_QuarksFromW_Prob_deltaR_pTcuts     -> Fill("< 0.8",0);
                         if(deltaR_obj1_obj2 < 0.8) h_NotHs_QuarksFromW_Prob_deltaR_pTcuts -> Fill("< 0.8",1);
                         else                       h_NotHs_QuarksFromW_Prob_deltaR_pTcuts -> Fill("> 0.8",1);
                       }
		     else
                       {
                         h_NotHs_QuarksFromW_deltaR                                -> Fill(deltaR_obj1_obj2);
                         h_NotHs_QuarksFromW_Prob_deltaR-> Fill("< 0.8",0);
                         if(deltaR_obj1_obj2 < 0.8) h_NotHs_QuarksFromW_Prob_deltaR-> Fill("< 0.8",1);
                         else                       h_NotHs_QuarksFromW_Prob_deltaR-> Fill("> 0.8",1);
                       }

                     if (deltaRmax < 0.8)
                       {
                         if(applyCut)
                           {
                             h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bqq",1);
                             h_NotHs_isbQuarkintoBaryCenter_pTcuts           -> Fill("b",1);
                           }
                         else
                           {
                             h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("bqq",1);
                             h_NotHs_isbQuarkintoBaryCenter           -> Fill("b",1);
                           }
		       }
		     else if(deltaR_obj1_obj2 < 0.8 && deltaR_b_fromTop_obj1 > 0.8 && deltaR_b_fromTop_obj2 > 0.8)
		       {
			 if(applyCut)
			   {
			     h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq",1);
			     h_NotHs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",1);
			   }
			 else
			   {
			     h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("qq",1);
			     h_NotHs_isbQuarkintoBaryCenter           -> Fill("No b",1);
			   }
		       }
		     else if(deltaR_obj1_obj2 < 0.8 && (deltaR_b_fromTop_obj1 < 0.8 || deltaR_b_fromTop_obj2 < 0.8))
                       {
                         if(applyCut)
                           {
                             h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq-bq",1);
                             h_NotHs_isbQuarkintoBaryCenter_pTcuts           -> Fill("?",1);
                           }
                         else
                           {
                             h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("qq-bq",1);
                             h_NotHs_isbQuarkintoBaryCenter           -> Fill("?",1);
                           }
                       }
                     else if(deltaR_obj1_obj2 > 0.8 && (deltaR_b_fromTop_obj1 < 0.8 || deltaR_b_fromTop_obj2 < 0.8))
                       {
                         if(applyCut)
                           {
                             h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bq",1);
                             h_NotHs_isbQuarkintoBaryCenter_pTcuts           -> Fill("b",1);
                           }
                         else
                           {
                             h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("bq",1);
                             h_NotHs_isbQuarkintoBaryCenter           -> Fill("b",1);
                           }
                       }
		     else
                       {
                         if(applyCut)
                           {
                             h_NotHs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("0",1);
                             h_NotHs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",1);
                           }
                         else
                           {
                             h_NotHs_QuarksintoBaryCenterMultiplicity -> Fill("0",1);
                             h_NotHs_isbQuarkintoBaryCenter           -> Fill("No b",1);
                           }
                       }
		     
		     if(applyCut) h_NotHs_OnlyQQ_dR_less_p7_pTcuts -> Fill("True",0);
                     else         h_NotHs_OnlyQQ_dR_less_p7        -> Fill("True",0);

                     if ((deltaR_b_fromTop_obj1 > 0.8) && (deltaR_b_fromTop_obj2 > 0.8) && (deltaR_obj1_obj2 < 0.7))  // %1   
                       {
                         if(applyCut) h_NotHs_OnlyQQ_dR_less_p7_pTcuts -> Fill("True",1);
                         else         h_NotHs_OnlyQQ_dR_less_p7        ->Fill("True",1);

                         double deltaRmin_fatJet_obj1 = 1e6; //give an initial, non sense value    
			 math::XYZTLorentzVector closest_toObj1_fatJet_p4(0,0,0,0);
                         int AK8jetSD_index = -1;
			 for(AK8JetsSoftDrop fatjet: fEvent.ak8jetsSoftDrop())
                           {
                             AK8jetSD_index++;
			     math::XYZTLorentzVector fatJet_p4(0,0,0,0);
                             fatJet_p4 = fatjet.p4();
                             if((fatJet_p4.pt() < 100.0) || (fatJet_p4.eta() > 2.4) || !fatjet.IDloose()) continue;
                             double deltaR_fatJet_obj1  = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,obj1_fromW_fromTop_p4);
                             if(deltaR_fatJet_obj1 < deltaRmin_fatJet_obj1)
                               {
                                 deltaRmin_fatJet_obj1 = deltaR_fatJet_obj1;
                                 closest_toObj1_fatJet_p4 = fatjet.p4();
                               }
                           }   //for fat-jets loop

			 if(deltaRmin_fatJet_obj1 <0.8) //%2     
                           {
                             double deltaR_fatJet_obj2  = ROOT::Math::VectorUtil::DeltaR(closest_toObj1_fatJet_p4,obj2_fromW_fromTop_p4);
                             h_NotHs_Prob_Diquark_match_with_fj ->Fill("Matched",0);
                             if(deltaR_fatJet_obj2 < 0.8)
                               {
                                 if(applyCut)
                                   {
                                     h_NotHs_Prob_Diquark_match_with_fj_pTcuts -> Fill("Matched",1);
                                     h_NotHs_MatchedWithDiquark_fj_pT_pTcuts   -> Fill(closest_toObj1_fatJet_p4.pt());
				     h_NotHs_MatchedWithDiquark_Prob_fj_pT_pTcuts -> Fill("<300 GeV",0);
                                     if(closest_toObj1_fatJet_p4.pt() < 300.0) h_NotHs_MatchedWithDiquark_Prob_fj_pT_pTcuts-> Fill("<300 GeV",1);
                                     else                                      h_NotHs_MatchedWithDiquark_Prob_fj_pT_pTcuts-> Fill(">300 GeV",1);
                                   }
                                 else
                                   {
                                     h_NotHs_Prob_Diquark_match_with_fj -> Fill("Matched",1);
                                     h_NotHs_MatchedWithDiquark_fj_pT -> Fill(closest_toObj1_fatJet_p4.pt());
				     h_NotHs_MatchedWithDiquark_Prob_fj_pT -> Fill("<300 GeV",0);
                                     if(closest_toObj1_fatJet_p4.pt() < 300.0) h_NotHs_MatchedWithDiquark_Prob_fj_pT-> Fill("<300 GeV",1);
                                     else                                      h_NotHs_MatchedWithDiquark_Prob_fj_pT-> Fill(">300 GeV",1);
                                   }
                               }
                             else if(applyCut) h_NotHs_Prob_Diquark_match_with_fj_pTcuts ->Fill("Not Matched",1);
                             else         h_NotHs_Prob_Diquark_match_with_fj        ->Fill("Not Matched",1);
                           } //%2
			 
		       } // %1  

                     else if(applyCut) h_NotHs_OnlyQQ_dR_less_p7_pTcuts -> Fill("False",1);
                     else              h_NotHs_OnlyQQ_dR_less_p7 -> Fill("False",1);
                   }// for applying cuts

		 
	       } //else the top is not from Higgs                 
           }     //if top               
       }         //for gen particles  
   }

 return;

}
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////kchristo/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////FInd the Last Copy Of a Particle////////////////////////////////////////////////
// by index/////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int BoostedTopologies::GetTheLastCopy(unsigned int firstcopy_index)
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

unsigned int BoostedTopologies::GetTheFirstCopy(unsigned int lastcopy_index)
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

vector<float> BoostedTopologies::GetMomentumTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
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


//vector<float> BoostedTopologies::GetMomentumTensorEigenValues2D(std::vector<math::XYZTLorentzVector> jets, float &Circularity){

  // For Circularity, the momentum tensor is the 2¡Á2 submatrix of Mjk, normalized by the sum of pT instead by the sum of       
  // This matrix has two eigenvalues Qi with 0 < Q1 < Q2. The following definition for the circularity C has been used:      
  //  C = 2 ¡Á min (Q1,Q2) / (Q1 +Q2)                                                                                         
  // The circularity C is especially interesting for hadron colliders because it only uses the momentum values in x and y direction transverse
  // to the beam line. So C is a two dimensional event shape variable and is therefore independent from a boost along z.                  
  // In addition, the normalization by the sum of the particle momenta makes C highly independent from energy calibration effects  (systematic uncertainty).                                                                                                                                              
  // C takes small values for linear and high values for circular events.                                       

  // Find the Momentum-Tensor EigenValues (E1, E2)                                                 
  //TMatrixDSym MomentumTensor = ComputeMomentumTensor2D(jets);
  //TMatrixDSymEigen eigen(MomentumTensor);
  //TVectorD eigenvals = eigen.GetEigenValues();

  // Store & Sort the eigenvalues                                                                                       
  //vector<float> eigenvalues(2);
  //eigenvalues.at(0) = eigenvals[0]; // Q1                                                                
  //eigenvalues.at(1) = eigenvals[1]; // Q2                                                   
  //sort( eigenvalues.begin(), eigenvalues.end(), std::greater<float>() );
  
  // Save the final eigenvalues      
  //float Q1 = eigenvalues.at(0);
  //float Q2 = eigenvalues.at(1);

  // Sanity check on eigenvalues: (Q1 > Q2)                                                                        
  //if (Q1 == 0) return eigenvalues;

  //bool bInequality = (Q1 > 0 && Q1 > Q2);
  //if ( !(bInequality) )
  // {
  //    throw hplus::Exception("LogicError") << "Failure of requirement that eigenvalues are ordered as Q1 >= Q2. Found Q1 = " << Q1 << ", Q2 " << Q2;
  //  }

  // Calculate circularity                                                                           
  //Circularity = 2*std::min(Q1, Q2)/(Q1+Q2); // is this definition correct?                  

  //if (0)
  //  {
  //    Table vars("Variable | Value | Allowed Range | Definition", "Text"); //LaTeX or Text        
  //    vars.AddRowColumn(0, "Circularity");
  //    vars.AddRowColumn(0, auxTools.ToString(Circularity) );
  //    vars.AddRowColumn(0, "0.0 <= C <= 1.0 ");
  //    vars.AddRowColumn(0, "C = 2 ¡Á min (Q1,Q2)/(Q1+Q2)");
  //    vars.Print();
  //  }

  //return eigenvalues;
//}



vector<float> BoostedTopologies::GetSphericityTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
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

TMatrixDSym BoostedTopologies::ComputeMomentumTensor(std::vector<math::XYZTLorentzVector> jets, double r)
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



TMatrixDSym BoostedTopologies::ComputeMomentumTensor2D(std::vector<math::XYZTLorentzVector> jets)
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


double BoostedTopologies::GetAlphaT(std::vector<math::XYZTLorentzVector> jets,
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

vector<GenJet> BoostedTopologies::GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCuts, std::vector<float> etaCuts, vector<genParticle> genParticlesToMatch)

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

vector<GenJet> BoostedTopologies::GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCuts, std::vector<float> etaCuts)
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

vector<genParticle> BoostedTopologies::GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy, const bool hasNoDaughters)
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
