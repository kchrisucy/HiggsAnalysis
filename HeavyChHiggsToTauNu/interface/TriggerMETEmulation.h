// -*- c++ -*-
#ifndef HiggsAnalysis_HeavyChHiggsToTauNu_TriggerMETEmulation_h
#define HiggsAnalysis_HeavyChHiggsToTauNu_TriggerMETEmulation_h

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/METReco/interface/MET.h"

#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/EventCounter.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

namespace HPlus {
  class HistoWrapper;
  class WrappedTH1;
  class WrappedTH2;

  class TriggerMETEmulation {
  public:
    /**
     * Class to encapsulate the access to the data members of
     * TauSelection. If you want to add a new accessor, add it here
     * and keep all the data of TauSelection private.
     */
    class Data {
    public:
      // The reason for pointer instead of reference is that const
      // reference allows temporaries, while const pointer does not.
      // Here the object pointed-to must live longer than this object.
      Data(const TriggerMETEmulation* triggerMETEmulation, bool passedEvent);
      ~Data();

      bool passedEvent() const { return fPassedEvent; }

    private:
      const TriggerMETEmulation *fTriggerMETEmulation;
      const bool fPassedEvent;
    };
    
    TriggerMETEmulation(const edm::ParameterSet& iConfig, EventCounter& eventCounter, HistoWrapper& histoWrapper);
    ~TriggerMETEmulation();

    Data analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

    const edm::Ptr<reco::MET> getSelectedMET() const {
      return fSelectedTriggerMET;
    }

  private:
    // Input parameters
    edm::InputTag fSrc;
    const double fmetEmulationCut;

    // Counters
    Count fmetEmulationCutCount;

    // Histograms
    WrappedTH1 *hMetBeforeEmulation;
    WrappedTH1 *hMetAfterEmulation;

    // Selected jets
    edm::Ptr<reco::MET> fSelectedTriggerMET;
  };
}

#endif
