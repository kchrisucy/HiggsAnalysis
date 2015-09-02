#include "DataFormat/interface/EventNPU.h"
#include "Framework/interface/BranchManager.h"

EventNPU::EventNPU():
  fNPU(nullptr)
{}
EventNPU::~EventNPU() {}

void EventNPU::setupBranches(BranchManager& mgr) {
  //mgr.book("nPU", &fNPU); // The MC number of PU vertices is not available for data
  mgr.book("nGoodOfflinePV", &fNPU);
  mgr.book("nPU", &fSimulatedNPU);
}
