// -*- c++ -*-
#include "EventSelection/interface/TauSelection.h"

#include "Framework/interface/Exception.h"
#include "Framework/interface/ParameterSet.h"
#include "EventSelection/interface/CommonPlots.h"
#include "DataFormat/interface/Event.h"
#include "Framework/interface/HistoWrapper.h"
#include "DataFormat/interface/HLTTau.h"
//#include "Framework/interface/makeTH.h"

#include "Math/VectorUtil.h"

#include <cmath>

TauSelection::Data::Data() 
: fRtau(-1.0) { }

TauSelection::Data::~Data() { }

const Tau& TauSelection::Data::getSelectedTau() const { 
  if (!hasIdentifiedTaus())
    throw hplus::Exception("Assert") << "You forgot to check if taus exist (hasIdentifiedTaus()), this message occurs when none exist!";
  return fSelectedTaus[0];
}

TauSelection::TauSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, commonPlots, postfix),
  bApplyTriggerMatching(config.getParameter<bool>("applyTriggerMatching")),
  fTriggerTauMatchingCone(config.getParameter<float>("triggerMatchingCone")),
  fTauPtCut(config.getParameter<float>("tauPtCut")),
  fTauEtaCut(config.getParameter<float>("tauEtaCut")),
  fTauLdgTrkPtCut(config.getParameter<float>("tauLdgTrkPtCut")),
  fTauNprongs(config.getParameter<int>("prongs")),
  fTauRtauCut(config.getParameter<float>("rtau")),
  bInvertTauIsolation(config.getParameter<bool>("invertTauIsolation")),
  // Event counter for passing selection
  cPassedTauSelection(eventCounter.addCounter("passed tau selection ("+postfix+")")),
  cPassedTauSelectionMultipleTaus(eventCounter.addCounter("multiple taus ("+postfix+")")),
  // Sub counters
  cSubAll(eventCounter.addSubCounter("tau selection ("+postfix+")", "All events")),
  cSubPassedTriggerMatching(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed trigger matching")),
  cSubPassedDecayMode(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed decay mode")),
  cSubPassedGenericDiscriminators(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed generic discriminators")),
  cSubPassedElectronDiscr(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed e discr")),
  cSubPassedMuonDiscr(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed mu discr")),
  cSubPassedPt(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed pt cut")),
  cSubPassedEta(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed eta cut")),
  cSubPassedLdgTrk(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed ldg.trk pt cut")),
  cSubPassedNprongs(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed nprongs")),
  cSubPassedIsolation(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed isolation")),
  cSubPassedRtau(eventCounter.addSubCounter("tau selection ("+postfix+")", "Passed Rtau"))
{ }

TauSelection::~TauSelection() { }

void TauSelection::bookHistograms(TDirectory* dir) {
  TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kDebug, dir, "tauSelection_"+sPostfix);
  hTriggerMatchDeltaR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "triggerMatchDeltaR", "Trigger match #DeltaR", 60, 0, 3.);
  hTauPtTriggerMatched = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "tauPtAll", "Tau pT, all", 40, 0, 400);
  hTauEtaTriggerMatched = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "tauEtaAll", "Tau eta, all", 50, -2.5, 2.5);
  hNPassed = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "tauNpassed", "Tau eta, all", 20, 0, 20);
}

TauSelection::Data TauSelection::silentAnalyze(const Event& event) {
  ensureSilentAnalyzeAllowed(event.eventID());
  // Disable histogram filling and counter
  disableHistogramsAndCounters();
  Data myData = privateAnalyze(event);
  enableHistogramsAndCounters();
  return myData;
}

TauSelection::Data TauSelection::analyze(const Event& event) {
  ensureAnalyzeAllowed(event.eventID());
  return privateAnalyze(event);
}

TauSelection::Data TauSelection::privateAnalyze(const Event& event) {
  Data output;
  cSubAll.increment();
  bool passedTriggerMatching = false;
  bool passedDecayMode = false;
  bool passedGenericDiscr = false;
  bool passedEdiscr = false;
  bool passedMuDiscr = false;
  bool passedPt = false;
  bool passedEta = false;
  bool passedLdgTrkPt = false;
  bool passedNprongs = false;
  bool passedIsol = false;
  bool passedRtau = false;
  
  // Cache vector of trigger tau 4-momenta
  std::vector<math::LorentzVectorT<double>> myTriggerTauMomenta;
  for (HLTTau p: event.triggerTaus()) {
    myTriggerTauMomenta.push_back(p.p4());
  }
  // Loop over taus
  for (Tau tau: event.taus()) {
    // Apply trigger matching
    if (!this->passTrgMatching(tau, myTriggerTauMomenta))
      continue;
    passedTriggerMatching = true;
    // Apply cut on decay mode
    if (!this->passDecayModeFinding(tau))
      continue;
    passedDecayMode = true;
    hTauPtTriggerMatched->Fill(tau.pt());
    hTauEtaTriggerMatched->Fill(tau.eta());
    // Generic discriminators
    if (!this->passGenericDiscriminators(tau))
      continue;
    passedGenericDiscr = true;
    // Electron discrimator
    if (!this->passElectronDiscriminator(tau))
      continue;
    passedEdiscr = true;
    // Muon discriminator
    if (!this->passMuonDiscriminator(tau))
      continue;
    passedMuDiscr = true;
    // Pt cut on tau
    if (!this->passPtCut(tau))
      continue;
    passedPt = true;
    // Eta cut on tau
    if (!this->passEtaCut(tau))
      continue;
    passedEta = true;
    // Ldg. track pt cut
    if (!this->passLdgTrkPtCut(tau))
      continue;
    passedLdgTrkPt = true;
    // Number of prongs
    if (!this->passNprongsCut(tau))
      continue;
    passedNprongs = true;
    // Apply tau isolation
    if (bInvertTauIsolation) {
      if (this->passIsolationDiscriminator(tau))
        passedIsol = true;
      if (passedIsol)
        continue;
    } else {
      if (!this->passIsolationDiscriminator(tau))
        continue;
      passedIsol = true;
    }
    // Apply cut on Rtau
    if (!this->passRtauCut(tau))
      continue;
    passedRtau = true;
    
    output.fSelectedTaus.push_back(tau);
  }
  // If there are multiple taus, choose the one with highest pT
  std::sort(output.fSelectedTaus.begin(), output.fSelectedTaus.end());
  
  // Fill data object
  if (output.fSelectedTaus.size())
    output.fRtau = getRtau(output.getSelectedTau());

  // Fill counters
  if (passedTriggerMatching)
    cSubPassedTriggerMatching.increment();
  if (passedDecayMode)
    cSubPassedDecayMode.increment();
  if (passedGenericDiscr)
    cSubPassedGenericDiscriminators.increment();
  if (passedEdiscr)
    cSubPassedElectronDiscr.increment();
  if (passedMuDiscr)
    cSubPassedMuonDiscr.increment();
  if (passedPt)
    cSubPassedPt.increment();
  if (passedEta)
    cSubPassedEta.increment();
  if (passedLdgTrkPt)
    cSubPassedLdgTrk.increment();
  if (passedNprongs)
    cSubPassedNprongs.increment();
  if (bInvertTauIsolation) {
    if (!passedIsol)
      cSubPassedIsolation.increment();
  } else {
    if (passedIsol)
      cSubPassedIsolation.increment();
  }  
  if (passedRtau)
    cSubPassedRtau.increment();
  if (output.fSelectedTaus.size() > 0)
    cPassedTauSelection.increment();
  if (output.fSelectedTaus.size() > 1)
    cPassedTauSelectionMultipleTaus.increment();
  
  // Return data object
  return output;
}

bool TauSelection::passTrgMatching(const Tau& tau, std::vector<math::LorentzVectorT<double>>& trgTaus) const {
  if (!bApplyTriggerMatching)
    return true;
  double myMinDeltaR = 9999.0;
  for (auto& p: trgTaus) {
    double myDeltaR = ROOT::Math::VectorUtil::DeltaR(p, tau.p4());
    myMinDeltaR = std::min(myMinDeltaR, myDeltaR);
  }
  hTriggerMatchDeltaR->Fill(myMinDeltaR);
  return (myMinDeltaR < fTriggerTauMatchingCone);
}

bool TauSelection::passNprongsCut(const Tau& tau) const {
  if (fTauNprongs > 0) {
    if (fTauNprongs == 13) {
      return (tau.nProngs() == 1 || tau.nProngs() == 3);
    } else {
      return (tau.nProngs() == fTauNprongs);
    }
  }
  return true;
}

double TauSelection::getRtau(const Tau& tau) const {
  // Calculate p_z of leading track
  double pz = std::sinh(static_cast<double>(tau.lTrkEta()))*static_cast<double>(tau.lTrkPt());
  // Calcualte p of leading track
  double p = std::sqrt(pz*pz + static_cast<double>(tau.lTrkPt()*tau.lTrkPt()));
  // Calculate Rtau
  double rtau = -1.0;
  double taup = tau.p4().P();
  if (taup > 0.0)
    rtau = p / taup;
  return rtau;
  
  
}