// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"
// User
//kchristo--------------------------
#include "Tools/interface/MCTools.h"
//----------------------------------
#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"
//ROOT
#include "TDirectory.h"
//kchristo---------------------
#include "Math/VectorUtil.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
//-----------------------------

//kchristo--------------------------
//typedef Particle<ParticleCollection<double> > genParticle;

//struct PtComparator
//{
//  bool operator() (const genParticle p1, const genParticle p2) const { return ( p1.pt() > p2.pt() ); }
//  bool operator() (const math::XYZTLorentzVector p1, const math::XYZTLorentzVector p2) const { return ( p1.pt() > p2.pt() ); }
//};

//-----------------------------------

class kcHplus2tbAnalysis: public BaseSelector {
public:
  explicit kcHplus2tbAnalysis(const ParameterSet& config, const TH1* skimCounters);
  virtual ~kcHplus2tbAnalysis() {}

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;
  //kchristo ------------------------------------
  //virtual vector<genParticle> GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, 
  //					      const bool isLastCopy=true, const bool hasNoDaughters=false);
  //---------------------------------------------
  //kchristo ----------------------------------------------------------------------------------------------------------------------- 
  virtual unsigned int GetTheLastCopy(unsigned int firstcopy_index);
  virtual unsigned int GetTheFirstCopy(unsigned int lastcopy_index);
  //--------------------------------------------------------------------------------------------------------------------------------

  
  

private:
  // Input parameters
  // const DirectionalCut<float> cfg_PrelimTopFitChiSqr;
  const DirectionalCut<double> cfg_PrelimTopMVACut;

  // Common plots
  CommonPlots fCommonPlots;
  // Event selection classes and event counters (in same order like they are applied)
  Count cAllEvents;
  Count cTrigger;
  METFilterSelection fMETFilterSelection;
  Count cVertexSelection;
  ElectronSelection fElectronSelection;
  MuonSelection fMuonSelection;
  TauSelection fTauSelection;
  JetSelection fJetSelection;
  BJetSelection fBJetSelection;
  Count cBTaggingSFCounter;
  METSelection fMETSelection;
  //TopologySelection fTopologySelection;
  TopSelectionBDT fTopSelection;
  Count cSelected;
  
  // Marina 
  FatJetSelection fFatJetSelection;
  FatJetSoftDropSelection fFatJetSoftDropSelection;
    
  // Non-common histograms
  // WrappedTH1 *hAssociatedTop_Pt;
  // WrappedTH1 *hAssociatedTop_Eta;
  
  //kchristo-----------------------------------------------------------------
  WrappedTH1 *h_IsFatJetBtagged;
  // b-tag fatjets
  WrappedTH1 *h_BtagFatJet_ClosestBjet_dR;
  WrappedTH1 *h_BtagFatJet_bDiscriminator;
  WrappedTH1 *h_BtagFatJet_Prob_dRmin_bjet_less_p8;
  WrappedTH1 *h_BtagFatJet_dRmin_bjet_less_p8_pT;
  WrappedTH1 *h_BtagFatJet_dRmin_bjet_more_p8_pT;
  WrappedTH1 *h_BtagFatJet_NumOfSubjets;

  // not b-tag fatjets
  WrappedTH1 *h_NotBtagFatJet_ClosestBjet_dR;
  WrappedTH1 *h_NotBtagFatJet_bDiscriminator;
  WrappedTH1 *h_NotBtagFatJet_Prob_dRmin_bjet_less_p8;
  WrappedTH1 *h_NotBtagFatJet_dRmin_bjet_less_p8_pT;
  WrappedTH1 *h_NotBtagFatJet_dRmin_bjet_more_p8_pT;
  WrappedTH1 *h_NotBtagFatJet_NumOfSubjets;
  
  // Study boosted topologies
  WrappedTH1 *h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts;
  WrappedTH1 *h_Hs_isbQuarkintoBaryCenter_pTcuts;
  WrappedTH1 *h_Hs_BoostedWcase_matchWithAk4;
  
  //------------------------------------------------------------------------
  
};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(kcHplus2tbAnalysis);

kcHplus2tbAnalysis::kcHplus2tbAnalysis(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    // cfg_PrelimTopFitChiSqr(config, "FakeBMeasurement.prelimTopFitChiSqrCut"),
    cfg_PrelimTopMVACut(config, "FakeBMeasurement.prelimTopMVACut"),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2tbAnalysis, fHistoWrapper),
    cAllEvents(fEventCounter.addCounter("all events")),
    cTrigger(fEventCounter.addCounter("passed trigger")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("passed PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fBJetSelection(config.getParameter<ParameterSet>("BJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cBTaggingSFCounter(fEventCounter.addCounter("b tag SF")),
    fMETSelection(config.getParameter<ParameterSet>("METSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    //fTopologySelection(config.getParameter<ParameterSet>("TopologySelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fTopSelection(config.getParameter<ParameterSet>("TopSelectionBDT"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cSelected(fEventCounter.addCounter("Selected Events")),  // Marina ","      
    // Marina                  
    fFatJetSelection(config.getParameter<ParameterSet>("FatJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fFatJetSoftDropSelection(config.getParameter<ParameterSet>("FatJetSoftDropSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "")
{ }


void kcHplus2tbAnalysis::book(TDirectory *dir) {

  
  // Book common plots histograms
  fCommonPlots.book(dir, isData());

  // Book histograms in event selection classes
  fMETFilterSelection.bookHistograms(dir);
  fElectronSelection.bookHistograms(dir);
  fMuonSelection.bookHistograms(dir);
  fTauSelection.bookHistograms(dir);
  fJetSelection.bookHistograms(dir);
  fBJetSelection.bookHistograms(dir);
  fMETSelection.bookHistograms(dir);
  //fTopologySelection.bookHistograms(dir);
  fTopSelection.bookHistograms(dir);
  
  // Marina
  fFatJetSelection.bookHistograms(dir);
  fFatJetSoftDropSelection.bookHistograms(dir);

  //kchristo-----------------------------------------
  TDirectory* th1 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH1");
  //------------------------------------------------
  
  // Book non-common histograms
  // hAssociatedTop_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "associatedTop_Pt", "Associated t pT;p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  // hAssociatedTop_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "associatedTop_Eta", "Associated t eta;#eta", nBinsEta, minEta, maxEta);
  
  //kchristo ------------------------------------
  h_IsFatJetBtagged          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "IsFatJetBtagged", " ", 2, 0.0, 2.0);

  //b-tag fatjets
  h_BtagFatJet_ClosestBjet_dR= fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "BtagFatJet_ClosestBjet_dR", ";#Delta#R_{min}",                                                                30,0.0, 3.0);
  h_BtagFatJet_bDiscriminator= fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "BtagFatJet_bDiscriminator"
							  , ";b discriminator", 100 , 0.0, 1.0);
  h_BtagFatJet_Prob_dRmin_bjet_less_p8 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "BtagFatJet_Prob_dRmin_bjet_less_p8"
                                                                    , ";#Delta#R_{min}", 2 , 0.0, 2.0);
  h_BtagFatJet_dRmin_bjet_less_p8_pT   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "BtagFatJet_dRmin_bjet_less_p8_pT"
								     , ";p_{T}", 100 , 0.0, 1000.0);
  h_BtagFatJet_dRmin_bjet_more_p8_pT   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "BtagFatJet_dRmin_bjet_more_p8_pT"
                                                                     , ";p_{T}", 100 , 0.0, 1000.0);
  h_BtagFatJet_NumOfSubjets  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "BtagFatJet_NumOfSubjets", ";Subjets",6,-0.5, 5.5);  

  //not b-tag fatjets
  h_NotBtagFatJet_ClosestBjet_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotBtagFatJet_ClosestBjet_dR", ";#Delta#R_{min}"                                                            , 50, 0.0, 5.0);
  h_NotBtagFatJet_bDiscriminator    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotBtagFatJet_bDiscriminator"
								 , ";b discriminator", 100 , 0.0, 1.0);
  h_NotBtagFatJet_Prob_dRmin_bjet_less_p8 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotBtagFatJet_Prob_dRmin_bjet_less_p8"
                                                                    , ";#Delta#R_{min}", 2 , 0.0, 2.0);
  h_NotBtagFatJet_dRmin_bjet_less_p8_pT   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotBtagFatJet_dRmin_bjet_less_p8_pT"
								    , ";p_{T}", 100 , 0.0, 1000.0);
  h_NotBtagFatJet_dRmin_bjet_more_p8_pT   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotBtagFatJet_dRmin_bjet_more_p8_pT"
								    , ";p_{T}", 30 , 0.0, 300.0);
  h_NotBtagFatJet_NumOfSubjets   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "NotBtagFatJet_NumOfSubjets", ";Subjets", 6,-0.5, 5.5);

  ////////////////////////Study Boosted Topologies/////////////////////////////////////
  h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_QuarksintoBaryCenterMultiplicity_pTcuts", " ", 5 , 0.0 , 5.0 );
  h_Hs_isbQuarkintoBaryCenter_pTcuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_isbQuarkintoBaryCenter_pTcuts", " " , 3 ,0.0 , 3.0 );
  h_Hs_BoostedWcase_matchWithAk4 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "Hs_BoostedWcase_matchWithAk4"," ", 4, 0.0, 4.0);
  //---------------------------------------------

  return;
}


void kcHplus2tbAnalysis::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void kcHplus2tbAnalysis::process(Long64_t entry) {
  //====== Initialize
  fCommonPlots.initialize();
  fCommonPlots.setFactorisationBinForEvent(std::vector<float> {});
  cAllEvents.increment();

  // Initialise MCTools object , kchristo 
  MCTools mcTools(fEvent);
  //

  //================================================================================================   
  // 1) Apply trigger 
  //================================================================================================   
  if (0) std::cout << "=== Trigger" << std::endl;
  if ( !(fEvent.passTriggerDecision()) ) return;
  
  cTrigger.increment();
  int nVertices = fEvent.vertexInfo().value();
  fCommonPlots.setNvertices(nVertices);
  fCommonPlots.fillControlPlotsAfterTrigger(fEvent);

  //================================================================================================   
  // 2) MET filters (to remove events with spurious sources of fake MET)
  //================================================================================================   
  if (0) std::cout << "=== MET Filter" << std::endl;
  const METFilterSelection::Data metFilterData = fMETFilterSelection.analyze(fEvent);
  if (!metFilterData.passedSelection()) return;
  fCommonPlots.fillControlPlotsAfterMETFilter(fEvent);  

  //================================================================================================   
  // 3) Primarty Vertex (Check that a PV exists)
  //================================================================================================   
  if (0) std::cout << "=== Vertices" << std::endl;
  if (nVertices < 1) return;
  cVertexSelection.increment();
  fCommonPlots.fillControlPlotsAtVertexSelection(fEvent);
  
  //================================================================================================   
  // 4) Electron veto (Fully hadronic + orthogonality)
  //================================================================================================   
  if (0) std::cout << "=== Electron veto" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  if (eData.hasIdentifiedElectrons()) return;

  //================================================================================================
  // 5) Muon veto (Fully hadronic + orthogonality)
  //================================================================================================
  if (0) std::cout << "=== Muon veto" << std::endl;
  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  if (muData.hasIdentifiedMuons()) return;

  //================================================================================================   
  // 6) Tau Veto (HToTauNu Orthogonality)
  //================================================================================================   
  if (0) std::cout << "=== Tau Veto" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  if (tauData.hasIdentifiedTaus() ) return;

  //================================================================================================
  // 7) Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyzeWithoutTau(fEvent);
  if (!jetData.passedSelection()) return;
  fCommonPlots.fillControlPlotsAfterTopologicalSelections(fEvent, true);
 
  //================================================================================================  
  // 8) BJet selection
  //================================================================================================
  if (0) std::cout << "=== BJet selection" << std::endl;
  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);
  if (!bjetData.passedSelection()) return;
  // fCommonPlots.fillControlPlotsAfterBJetSelection(fEvent, bjetData);

  //================================================================================================  
  // 9) BJet SF  
  //================================================================================================
  if (0) std::cout << "=== BJet SF" << std::endl;
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
    }
  cBTaggingSFCounter.increment();

  //================================================================================================
  // 10) MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  const METSelection::Data METData = fMETSelection.analyze(fEvent, nVertices);
  // if (!METData.passedSelection()) return;

  //================================================================================================
  // 11) Topology selection
  //================================================================================================
  //if (0) std::cout << "=== Topology selection" << std::endl;
  //const TopologySelection::Data topologyData = fTopologySelection.analyze(fEvent, jetData);
  // if (!topologyData.passedSelection()) return; 

  //================================================================================================
  // 12) Top selection
  //================================================================================================
  /*
    if (0) std::cout << "=== Top (ChiSq) selection" << std::endl;
    const TopSelection::Data topData = fTopSelection.analyze(fEvent, jetData, bjetData);
    // Apply preliminary chiSq cut
    bool passPrelimChiSq = cfg_PrelimTopFitChiSqr.passedCut(topData.ChiSqr());
    if (!passPrelimChiSq) return;
  */
  if (0) std::cout << "=== Top (BDT) selection" << std::endl;
  const TopSelectionBDT::Data topData = fTopSelection.analyze(fEvent, jetData, bjetData, true);
  bool passPrelimMVACut = cfg_PrelimTopMVACut.passedCut( std::max(topData.getMVAmax1(), topData.getMVAmax2()) ); //fixme?
  bool hasFreeBJet      = topData.hasFreeBJet();
  if (!hasFreeBJet) return;
  if (!passPrelimMVACut) return;

  //================================================================================================
  // Standard Selections
  //================================================================================================
  if (0) std::cout << "=== Standard Selections" << std::endl;
  //fCommonPlots.fillControlPlotsAfterStandardSelections(fEvent, jetData, bjetData, METData, jetData, topData, bjetData.isGenuineB());
  fCommonPlots.fillControlPlotsAfterStandardSelections(fEvent, jetData, bjetData, METData, TopologySelection::Data(), topData, bjetData.isGenuineB());

  //================================================================================================
  // All Selections
  //================================================================================================
  //if (!topologyData.passedSelection()) return;
  if (!topData.passedSelection()) return;

  if (0) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();


  //================================kchristo========================================================  
  // AK4-AK8 Analysis
  //================================================================================================ 
  // only bjets---------------------------------------
  //int index = -1;
  // For-loop: All selected bjets                                                                                                                                                                                   
  //for (auto bjet: bjetData.getSelectedBJets())
  //  {
  //    index++;
  //  }
  //--------------------------------------------------

  int AK8jetSD_index = -1;
  for(AK8JetsSoftDrop fatjet: fEvent.ak8jetsSoftDrop())
    {
      AK8jetSD_index++;
      math::XYZTLorentzVector fatJet_p4(0,0,0,0);
      fatJet_p4 = fatjet.p4();
      math::XYZTLorentzVector bJet_p4(0,0,0,0);
      double deltaRmin_fatJet_bjet = 1e6; //give an initial, non sense value
      bool FJhasBTagSubjet = fatjet.hasBTagSubjets();
      
      int AK4bjet_index = -1;                                                    
      // For-loop: All selected bjets                   
      
      for (auto bjet: bjetData.getSelectedBJets())                                
	{                                                                     
	  AK4bjet_index++;
	  bJet_p4 = bjet.p4();
	  double deltaR_fatJet_bjet  = ROOT::Math::VectorUtil::DeltaR(fatJet_p4,bJet_p4);
	  if(deltaR_fatJet_bjet < deltaRmin_fatJet_bjet) deltaRmin_fatJet_bjet = deltaR_fatJet_bjet;
	}  
      
      h_IsFatJetBtagged                       -> Fill("b-tag",0);
      h_BtagFatJet_Prob_dRmin_bjet_less_p8    -> Fill("< 0.8",0);
      h_NotBtagFatJet_Prob_dRmin_bjet_less_p8 -> Fill("< 0.8",0);
      
      if(FJhasBTagSubjet) 
	{
	  h_IsFatJetBtagged              -> Fill("b-tag",1);
	  h_BtagFatJet_bDiscriminator    -> Fill(fatjet.pfCombinedInclusiveSecondaryVertexV2BJetTags());
	  h_BtagFatJet_ClosestBjet_dR    -> Fill(deltaRmin_fatJet_bjet);
	  if(deltaRmin_fatJet_bjet < 0.8 )  
	    {
	      h_BtagFatJet_Prob_dRmin_bjet_less_p8 -> Fill("< 0.8",1);
	      h_BtagFatJet_dRmin_bjet_less_p8_pT   -> Fill(fatJet_p4.pt());
	    }
	  else
	    {
	      h_BtagFatJet_Prob_dRmin_bjet_less_p8 -> Fill("> 0.8",1);
	      h_BtagFatJet_dRmin_bjet_more_p8_pT   -> Fill(fatJet_p4.pt());
	    }
	  h_BtagFatJet_NumOfSubjets      -> Fill(fatjet.nSubjets());
	}
      
      else // if does not have btagsubjet
	{
	  h_IsFatJetBtagged              -> Fill("Not b-tag",1);
	  h_NotBtagFatJet_bDiscriminator -> Fill(fatjet.pfCombinedInclusiveSecondaryVertexV2BJetTags());
	  h_NotBtagFatJet_ClosestBjet_dR ->Fill(deltaRmin_fatJet_bjet);
	  if(deltaRmin_fatJet_bjet < 0.8 )
	    {
	      h_NotBtagFatJet_Prob_dRmin_bjet_less_p8 -> Fill("< 0.8",1);
	      h_NotBtagFatJet_dRmin_bjet_less_p8_pT   -> Fill(fatJet_p4.pt());
	    }
	  else     
	    {
	      h_NotBtagFatJet_Prob_dRmin_bjet_less_p8 -> Fill("> 0.8",1);
	      h_NotBtagFatJet_dRmin_bjet_more_p8_pT   -> Fill(fatJet_p4.pt());
	    }
	  h_NotBtagFatJet_NumOfSubjets   -> Fill(fatjet.nSubjets());
	}
       

      //////kchristo/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////// Study Boosted Topologies//////////////////////////////////////////////////////// 
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
      if(0)
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
				  genParticle top_grand_d;
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
		      ///////////////////////////////////////////////////////////////////////////////////////
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
		      
		      /////////////////////////////////////////////////BaryCenter////////////////////////////////////////////////////////
		      // barycenter not the top-direction            
		      h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("0",0);  //just to determinate the label of the first bin
		      h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq",0);
		      h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bq",0);
		      h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq-bq",0);

		      h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",0);  //just to determinate the label of the first bin
		      h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("?",0);
		      
		      double quarksPtCut = 30.0;
		      if((obj1_fromW_fromTop_p4.pt() < quarksPtCut) || (obj2_fromW_fromTop_p4.pt() < quarksPtCut) ||  (b_fromTop_fromH_p4.pt() < quarksPtCut)) continue;

		      if (deltaRmax < 0.8)
			{
			  h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bqq",1);
			  h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("b",1);
			}
		      else if(deltaR_obj1_obj2 < 0.8 && deltaR_b_fromTop_obj1 > 0.8 && deltaR_b_fromTop_obj2 > 0.8)
			{
			  h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq",1);
			  h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",1);

			  // how many times can I match the 2 with different ak4                                                            
			  int AK4jet_index = -1;
			  //std::vector<math::XYZTLorentzVector> Matched_Ak4jets_p4(2);
			  int MatchedAk4_withobj1_index = 0, MatchedAk4_withobj2_index = 0;
			  double deltaRmin_Ak4jet_obj1 = 1e6, deltaRmin_Ak4jet_obj2 = 1e6; //give an initial, non sense value               
			  // For-loop: All selected jets                                                                                    
			  for (const Jet& Ak4jet: jetData.getSelectedJets())
			    {
			      AK4jet_index++;
			      math::XYZTLorentzVector Ak4Jet_p4;
			      Ak4Jet_p4 = Ak4jet.p4();
			      double deltaR_Ak4jet_obj1 = ROOT::Math::VectorUtil::DeltaR(Ak4Jet_p4,obj1_fromW_fromTop_p4);
			      double deltaR_Ak4jet_obj2 = ROOT::Math::VectorUtil::DeltaR(Ak4Jet_p4,obj2_fromW_fromTop_p4);
			      if(deltaR_Ak4jet_obj1 < deltaRmin_Ak4jet_obj1)
				{
				  deltaRmin_Ak4jet_obj1 = deltaR_Ak4jet_obj1;
				  MatchedAk4_withobj1_index = AK4jet_index++;
				}
			      if(deltaR_Ak4jet_obj2 < deltaRmin_Ak4jet_obj2)
				{
				  deltaRmin_Ak4jet_obj2 = deltaR_Ak4jet_obj2;
				  MatchedAk4_withobj2_index = AK4jet_index++;
				}
			    }// for AK4 loop
			  h_Hs_BoostedWcase_matchWithAk4 -> Fill("2 Ak4 jets",0);
			  h_Hs_BoostedWcase_matchWithAk4 -> Fill("Same Ak4 jet",0);
			  h_Hs_BoostedWcase_matchWithAk4 -> Fill("Single match",0);

			  if(deltaRmin_Ak4jet_obj1 < 0.4 && deltaRmin_Ak4jet_obj2 < 0.4 && MatchedAk4_withobj1_index == MatchedAk4_withobj2_index)                         h_Hs_BoostedWcase_matchWithAk4 -> Fill("Same Ak4 jet",1);
			  else if(deltaRmin_Ak4jet_obj1 < 0.4 && deltaRmin_Ak4jet_obj2 < 0.4 && MatchedAk4_withobj1_index != MatchedAk4_withobj2_index)                    h_Hs_BoostedWcase_matchWithAk4 -> Fill("2 Ak4 jets",1);
			  else if(deltaRmin_Ak4jet_obj1< 0.4 || deltaRmin_Ak4jet_obj2 < 0.4)
			     h_Hs_BoostedWcase_matchWithAk4 -> Fill("Single match",1);
			  else h_Hs_BoostedWcase_matchWithAk4 -> Fill("No match",1);
			}
		      
		      else if(deltaR_obj1_obj2 < 0.8 && (deltaR_b_fromTop_obj1 < 0.8 || deltaR_b_fromTop_obj2 < 0.8))
			{
			  h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("qq-bq",1);
			  h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("?",1); 
			}
		      else if(deltaR_obj1_obj2 > 0.8 && (deltaR_b_fromTop_obj1 < 0.8 || deltaR_b_fromTop_obj2 < 0.8))
			{
			  h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("bq",1);
			  h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("b",1);
			}
		      else
			{
			  h_Hs_QuarksintoBaryCenterMultiplicity_pTcuts -> Fill("0",1);
			  h_Hs_isbQuarkintoBaryCenter_pTcuts           -> Fill("No b",1);
			}
		    }// if top from Higgs
		}//     find the top
	    }    //     for gen particls
	} 
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      //      math::XYZTLorentzVector top_fromH_p4(0,0,0,0), top_NotfromH_p4(0,0,0,0);
      //for (auto& p: fEvent.genparticles().getGenParticles())
      //{
      //  if(std::abs(p.pdgId()) == 6 && p.isFirstCopy()) // find the top
	    //    {
      //std::vector<short> top_mothers = p.mothers();
      //      genParticle m = fEvent.genparticles().getGenParticles()[top_mothers.at(0)]; //create the particle object
      //      if(std::abs(m.pdgId()) == 37)//if the top comes from Higgs
		//	{
      //top_fromH_p4 = p.p4();
      //	}
      //      else //if the top is not coming from H
      //	{
      //	  top_NotfromH_p4 = p.p4();
      //	}
      //    }      //find the top
      // }         // for gen Particles
      //      std::cout << "--------Event Ends----------" << std::endl;
	  
      //loop fro all jets (bjets or not)
      //int AK4jet_index = -1;
      // Loop over selected jets      
      //for(const Jet& jet: jetData.getSelectedJets())
      //	{
      //	  AK4jet_index++;
      //	}
    }
  
  //================================================================================================
  // Fill final plots
  //===============================================================================================
  fCommonPlots.fillControlPlotsAfterAllSelections(fEvent, 1);
 
  //================================================================================================
  // Finalize
  //================================================================================================
  fEventSaver.save();

  return;
}

///////kchristo/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////FInd the Last Copy Of a Particle////////////////////////////////////////////////   
// by index/////////////////////////////////////////////////////////////////////////////////////////////////////   

unsigned int kcHplus2tbAnalysis::GetTheLastCopy(unsigned int firstcopy_index)
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

///////kchristo/////////////////////////////////////////////////////////////////////////////////////////////////  
////////////////////////////////FInd the First Copy Of a Particle////////////////////////////////////////////////   
// by index/////////////////////////////////////////////////////////////////////////////////////////////////////    

unsigned int kcHplus2tbAnalysis::GetTheFirstCopy(unsigned int lastcopy_index)
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
