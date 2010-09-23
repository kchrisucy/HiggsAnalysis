#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include "FWCore/Utilities/interface/InputTag.h"

/* class CaloRecoTauDiscriminationByNProngs
 * created : September 23 2010,
 * contributors : Sami Lehti (sami.lehti@cern.ch ; HIP, Helsinki)
 * based on H+ tau ID by Lauri Wendland
 */

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include "TLorentzVector.h"

using namespace reco;
using namespace std;

class CaloRecoTauDiscriminationByNProngs : public CaloTauDiscriminationProducerBase  {
    public:
	explicit CaloRecoTauDiscriminationByNProngs(const ParameterSet& iConfig):CaloTauDiscriminationProducerBase(iConfig){

		nprongs			= iConfig.getParameter<uint32_t>("nProngs");
		booleanOutput = iConfig.getParameter<bool>("BooleanOutput");
	}

      	~CaloRecoTauDiscriminationByNProngs(){}

	void beginEvent(const edm::Event&, const edm::EventSetup&);
	double discriminate(const reco::CaloTauRef&);

    private:

//	CaloTauQualityCutWrapper qualityCuts_;

	uint32_t nprongs;
	bool booleanOutput;
};

void CaloRecoTauDiscriminationByNProngs::beginEvent(const Event& iEvent, const EventSetup& iSetup){}

double CaloRecoTauDiscriminationByNProngs::discriminate(const CaloTauRef& tau){

	bool accepted = false;
	int np = tau->signalTracks().size();

	if((np == 1 && (nprongs == 1 || nprongs == 0)) ||
           (np == 3 && (nprongs == 3 || nprongs == 0)) ) accepted = true;

	if(!accepted) np = 0;
	if(booleanOutput) return accepted;
	return np;
}

DEFINE_FWK_MODULE(CaloRecoTauDiscriminationByNProngs);

