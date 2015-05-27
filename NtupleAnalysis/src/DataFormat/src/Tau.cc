#include "DataFormat/interface/Tau.h"

#include "Framework/interface/BranchManager.h"

void TauCollection::setupBranches(BranchManager& mgr) {
  TauGeneratedCollection::setupBranches(mgr);

  fConfigurableDiscriminators.resize(fConfigurableDiscriminatorNames.size());
  size_t i=0;
  for(const auto& name: fConfigurableDiscriminatorNames) {
    mgr.book(prefix()+"_"+name, &(fConfigurableDiscriminators[i]));
    ++i;
  }
}

