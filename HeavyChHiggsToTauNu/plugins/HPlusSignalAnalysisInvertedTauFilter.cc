#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/EventCounter.h"
#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/EventWeight.h"
#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/SignalAnalysisInvertedTau.h"

class HPlusSignalAnalysisInvertedTauFilter: public edm::EDFilter {
 public:

  explicit HPlusSignalAnalysisInvertedTauFilter(const edm::ParameterSet&);
  ~HPlusSignalAnalysisInvertedTauFilter();

 private:
  virtual void beginJob();
  virtual bool filter(edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob();

  virtual bool endLuminosityBlock(edm::LuminosityBlock& iBlock, const edm::EventSetup & iSetup);

  HPlus::EventWeight eventWeight;
  HPlus::EventCounter eventCounter;
  HPlus::SignalAnalysisInvertedTau analysis;
};

HPlusSignalAnalysisInvertedTauFilter::HPlusSignalAnalysisInvertedTauFilter(const edm::ParameterSet& pset):
  eventWeight(pset), eventCounter(pset, eventWeight), analysis(pset, eventCounter, eventWeight)
{
  analysis.produces(this);
}
HPlusSignalAnalysisInvertedTauFilter::~HPlusSignalAnalysisInvertedTauFilter() {}
void HPlusSignalAnalysisInvertedTauFilter::beginJob() {}

bool HPlusSignalAnalysisInvertedTauFilter::endLuminosityBlock(edm::LuminosityBlock& iBlock, const edm::EventSetup & iSetup) {
  eventCounter.endLuminosityBlock(iBlock, iSetup);
  return true;
}

bool HPlusSignalAnalysisInvertedTauFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  return analysis.filter(iEvent, iSetup);
}

void HPlusSignalAnalysisInvertedTauFilter::endJob() {
  eventCounter.endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPlusSignalAnalysisInvertedTauFilter);
