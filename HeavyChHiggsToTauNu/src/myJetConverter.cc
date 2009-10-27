
#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/MyEventConverter.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

MyJet MyEventConverter::myJetConverter(const Muon& recMuon){

        MyJet muon(recMuon.px(), recMuon.py(), recMuon.pz(), recMuon.p()); // FIXME: should we use .energy() instead of .p()?
        muon.type = 13 * recMuon.charge();

	TrackRef track = recMuon.globalTrack();
	if(track.isNull()) track = recMuon.innerTrack();
	if(track.isNull()) track = recMuon.combinedMuon();

	if(track.isNonnull()){
		const TransientTrack transientTrack = transientTrackBuilder->build(track);

		MyTrack muonTrack = myTrackConverter(transientTrack);
		muonTrack.ip = impactParameter(transientTrack);
		muon.tracks.push_back(muonTrack);

		muon.tracks = getTracks(muon);
	}

	muon.tagInfo = muonTag(recMuon);

	return muon;
}

MyJet MyEventConverter::myJetConverter(const pat::Muon& recMuon){

        MyJet muon(recMuon.px(), recMuon.py(), recMuon.pz(), recMuon.p()); // FIXME: should we use .energy() instead of .p()?
        muon.type = 13 * recMuon.charge();

        TrackRef track = recMuon.globalTrack();
        if(track.isNull()) track = recMuon.innerTrack();

        if(track.isNonnull()){
                const TransientTrack transientTrack = transientTrackBuilder->build(track);

                MyTrack muonTrack = myTrackConverter(transientTrack);
                muonTrack.ip = impactParameter(transientTrack);
                muon.tracks.push_back(muonTrack);

                muon.tracks = getTracks(muon);
        }

	muon.tagInfo = muonTag(recMuon);

        return muon;
}

MyJet MyEventConverter::myJetConverter(const GsfElectron* recElectron){
	GsfTrackRef track = recElectron->gsfTrack();
        const TransientTrack transientTrack = transientTrackBuilder->build(track);

        MyJet electron(recElectron->px(), recElectron->py(), recElectron->pz(), recElectron->p()); // FIXME: should we use .energy() instead of .p()?
        electron.type = 11 * (*track).charge();

	MyTrack electronTrack = myTrackConverter(transientTrack);
	electronTrack.ip = impactParameter(transientTrack);
	electronTrack.trackEcalHitPoint = trackEcalHitPoint(recElectron);
	electron.tracks.push_back(electronTrack);
        electron.tracks = getTracks(electron);

	vector<TLorentzVector> superClusters;
	superClusters.push_back(TLorentzVector(recElectron->superCluster()->x(),
	                                       recElectron->superCluster()->y(),
                                               recElectron->superCluster()->z(),
                                               recElectron->superCluster()->energy()));
	electron.clusters = superClusters;

        return electron;
}

MyJet MyEventConverter::myJetConverter(const pat::Electron& recElectron){
        GsfTrackRef track = recElectron.gsfTrack();
        const TransientTrack transientTrack = transientTrackBuilder->build(track);

	MyJet electron(recElectron.px(), recElectron.py(), recElectron.pz(), recElectron.p()); // FIXME: should we use .energy() instead of .p()?
	electron.type = 11 * (*track).charge();

        MyTrack electronTrack = myTrackConverter(transientTrack);
        electronTrack.ip = impactParameter(transientTrack);
	electronTrack.trackEcalHitPoint = trackEcalHitPoint(&recElectron);
        electron.tracks.push_back(electronTrack);
        electron.tracks = getTracks(electron);

        vector<TLorentzVector> superClusters;
        superClusters.push_back(TLorentzVector(recElectron.superCluster()->x(),
                                               recElectron.superCluster()->y(),
                                               recElectron.superCluster()->z(),
                                               recElectron.superCluster()->energy()));
        electron.clusters = superClusters;

	electron.tagInfo = etag(recElectron);

	return electron;
}

MyJet MyEventConverter::myJetConverter(const Photon* recPhoton){

        MyJet photon(recPhoton->px(), recPhoton->py(), recPhoton->pz(), recPhoton->p()); // FIXME: should we use .energy() instead of .p()?
        photon.type = 0; //unconverted

        photon.tracks = getTracks(photon);

	photon.tagInfo = photontag(recPhoton);

        return photon;
}

MyJet MyEventConverter::myJetConverter(const Conversion* recPhoton){

        const GlobalVector& mom(recPhoton->pairMomentum());
        MyJet photon(mom.x(), mom.y(), mom.z(), mom.mag());
	photon.type = 1; //converted

        vector<MyTrack> tracks;
	vector<TrackRef> associatedTracks = recPhoton->tracks();
	vector<TrackRef>::const_iterator iTrack;
        for(iTrack = associatedTracks.begin(); iTrack!= associatedTracks.end(); ++iTrack){

                const TransientTrack transientTrack = transientTrackBuilder->build(**iTrack);

                MyTrack track           = myTrackConverter(transientTrack);
                track.ip                = impactParameter(transientTrack,recPhoton);
                track.trackEcalHitPoint = trackEcalHitPoint(transientTrack,recPhoton);
                tracks.push_back(track);
        }
        photon.tracks = tracks;

        photon.tagInfo = photontag(recPhoton);

        return photon;
}

MyJet MyEventConverter::myJetConverter(const CaloJet* caloJet){

        MyJet jet(caloJet->px(), caloJet->py(), caloJet->pz(), caloJet->energy());
        jet.tracks = getTracks(jet);

        // Jet energy corrections
        for(unsigned int i = 0; i < jetEnergyCorrectionTypes.size(); ++i){
                double jetEnergyCorrectionFactor = jetEnergyCorrections[i]->correction(*caloJet);
                string jetEnergyCorrectionName = jetEnergyCorrectionTypes[i].label();
                jet.addEnergyCorrection(jetEnergyCorrectionName,jetEnergyCorrectionFactor);
		cout << "    jet correction " << jetEnergyCorrectionName << " " 
                                              << jetEnergyCorrectionFactor << endl;
        }

        return jet;
}

MyJet MyEventConverter::myJetConverter(const pat::Jet* recoJet){

        MyJet jet(recoJet->px(), recoJet->py(), recoJet->pz(), recoJet->energy());
        jet.tracks = getTracks(jet);

	return jet;
}

MyJet MyEventConverter::myJetConverter(const JetTag& recJet){
        const CaloJet* caloJet = dynamic_cast<const CaloJet*>(recJet.first.get());
        return myJetConverter(caloJet);
}

MyJet MyEventConverter::myJetConverter(const IsolatedTauTagInfo& recTau){

	const CaloJet* caloJet = dynamic_cast<const CaloJet*>(recTau.jet().get());

        MyJet tau(caloJet->px(), caloJet->py(), caloJet->pz(), caloJet->energy());

	const TrackRefVector associatedTracks = recTau.allTracks();
	RefVector<TrackCollection>::const_iterator iTrack;
	vector<MyTrack> tracks;
	for(iTrack = associatedTracks.begin(); iTrack!= associatedTracks.end(); iTrack++){

                const TransientTrack transientTrack = transientTrackBuilder->build(*iTrack);

		MyTrack track           = myTrackConverter(transientTrack);
		track.ip                = impactParameter(transientTrack,caloJet);
		track.trackEcalHitPoint = trackEcalHitPoint(**iTrack, caloJet);
		tracks.push_back(track);
	}
	tau.tracks = tracks;

        tau.tagInfo = tauTag(recTau);

	tau.caloInfo = caloTowers(*caloJet);

	addECALClusters(&tau);

        return tau;
}

MyJet MyEventConverter::myJetConverter(const CaloTau& recTau){


        const CaloJet* caloJet = recTau.caloTauTagInfoRef()->calojetRef().get();

	MyJet tau(recTau.px(), recTau.py(), recTau.pz(), recTau.energy());

	vector<MyTrack> tracks;
	vector<MyHit> hits;

	vector<TransientTrack> transientTracks;
	if(trackCollectionSelection.label() == "iterativeTracks"){
	  vector<Trajectory> associatedTrajectories;
	  vector<Track> associatedTracks = tracksInCone(recTau.p4(),0.5,&associatedTrajectories);
	  vector<Track>::const_iterator iTrack;
	  vector<Trajectory>::const_iterator iTrajectory = associatedTrajectories.begin();
          // Make sure, that each track has a trajectory; only this guarantees one to one correspondence
          bool myTrajectoryStatus = (associatedTracks.size() == associatedTrajectories.size());
	  int trackCounter = 0;
          for(iTrack = associatedTracks.begin(); iTrack!= associatedTracks.end(); iTrack++){

                const TransientTrack transientTrack = transientTrackBuilder->build(*iTrack);
		transientTracks.push_back(transientTrack);

		MyTrack track = myTrackConverter(transientTrack);

                if (myTrajectoryStatus) {
			vector<MyHit> assocHits = getHits(*iTrajectory,trackCounter);
			hits.insert(hits.end(),assocHits.begin(),assocHits.end());
		}
                track.ip                = impactParameter(transientTrack,caloJet);
                track.trackEcalHitPoint = trackEcalHitPoint(*iTrack, caloJet);
                tracks.push_back(track);
		++iTrajectory;
		++trackCounter;
          }
	}else{
	  // at this point, adding MyHit information is not implemented for calotau data
	  const TrackRefVector associatedTracks = recTau.caloTauTagInfoRef()->Tracks();
	  RefVector<TrackCollection>::const_iterator iTrack;
	  for(iTrack = associatedTracks.begin(); iTrack!= associatedTracks.end(); iTrack++){

                const TransientTrack transientTrack = transientTrackBuilder->build(*iTrack);

                MyTrack track           = myTrackConverter(transientTrack);
                track.ip                = impactParameter(transientTrack,caloJet);
                track.trackEcalHitPoint = trackEcalHitPoint(**iTrack, caloJet);
                tracks.push_back(track);
          }
	}

        tau.tracks = tracks;

	tau.hits   = hits;

        tau.tagInfo = tauTag(recTau);

        // Jet energy correction
        double jetEnergyCorrectionFactor = tauJetCorrection->correction(recTau.p4());
        tau.addEnergyCorrection("TauJet",jetEnergyCorrectionFactor);

        tau.caloInfo = caloTowers(*caloJet);

	tau.secVertices = secondaryVertices(transientTracks);

	addECALClusters(&tau);

        return tau;
}

MyJet MyEventConverter::myJetConverter(const pat::Tau& recTau){

        MyJet tau(recTau.px(), recTau.py(), recTau.pz(), recTau.energy());

        vector<MyTrack> tracks;
        const PFCandidateRefVector pfSignalCandidates = recTau.signalPFCands();

	vector<TransientTrack> transientTracks;
        RefVector<PFCandidateCollection>::const_iterator iTrack;
        for(iTrack = pfSignalCandidates.begin(); iTrack!= pfSignalCandidates.end(); iTrack++){

                const PFCandidate* pfCand = iTrack->get();
		const TransientTrack transientTrack = transientTrackBuilder->build(pfCand->trackRef());
                transientTracks.push_back(transientTrack);

                MyTrack track = myTrackConverter(pfCand);
                tracks.push_back(track);
        }

        const PFCandidateRefVector pfIsolCandidates = recTau.isolationPFCands();
        for(iTrack = pfIsolCandidates.begin(); iTrack!= pfIsolCandidates.end(); iTrack++){

                const PFCandidate* pfCand = iTrack->get();
		const TransientTrack transientTrack = transientTrackBuilder->build(pfCand->trackRef());
                transientTracks.push_back(transientTrack);

                MyTrack track = myTrackConverter(pfCand);
                tracks.push_back(track);
        }

        tau.tracks = tracks;

        tau.tagInfo = tauTag(recTau);

	tau.secVertices = secondaryVertices(transientTracks);
	addECALClusters(&tau);

	return tau;
}

MyJet MyEventConverter::myJetConverter(const PFTau& recTau){

	MyJet tau(recTau.px(), recTau.py(), recTau.pz(), recTau.energy());

	vector<MyTrack> tracks;
        const PFCandidateRefVector pfSignalCandidates = recTau.signalPFCands();

	vector<TransientTrack> transientTracks;
        RefVector<PFCandidateCollection>::const_iterator iTrack;
        for(iTrack = pfSignalCandidates.begin(); iTrack!= pfSignalCandidates.end(); iTrack++){

		const PFCandidate* pfCand = iTrack->get();
		if(pfCand->trackRef().isNonnull()){
                  const TransientTrack transientTrack = transientTrackBuilder->build(pfCand->trackRef());
                  transientTracks.push_back(transientTrack);
		}
		MyTrack track = myTrackConverter(pfCand);
		track.trackEcalHitPoint = trackEcalHitPoint(pfCand);
                tracks.push_back(track);
        }

        const PFCandidateRefVector pfIsolCandidates = recTau.isolationPFCands();
        for(iTrack = pfIsolCandidates.begin(); iTrack!= pfIsolCandidates.end(); iTrack++){

		const PFCandidate* pfCand = iTrack->get();
		if(pfCand->trackRef().isNonnull()){
                  const TransientTrack transientTrack = transientTrackBuilder->build(pfCand->trackRef());
                  transientTracks.push_back(transientTrack);
		}
                MyTrack track = myTrackConverter(pfCand);
		track.trackEcalHitPoint = trackEcalHitPoint(pfCand);
                tracks.push_back(track);
        }

        tau.tracks = tracks;

        tau.tagInfo = tauTag(recTau);

	tau.secVertices = secondaryVertices(transientTracks);
	addECALClusters(&tau);

	return tau;	
}

void MyEventConverter::addECALClusters(MyJet* jet) {
  // Loops over barrel and endcap ECAL cluster
  // and stores to jet those, which are within specified DR to
  // leading track hit point on ECAL surface

  const MyTrack *myLeadingTrack = jet->leadingTrack();
  if (!myLeadingTrack || myLeadingTrack->Pt() < 0.0001) return;
  MyGlobalPoint myECALHitPoint = myLeadingTrack->ecalHitPoint();
  //double myLdgEta = myECALHitPoint.Eta();
  //double myLdgPhi = myECALHitPoint.Phi();

  // Loop over barrel ECAL clusters
  unsigned int myBarrelCollectionSize = theBarrelBCCollection->size();
  for(unsigned int i_BC=0; i_BC != myBarrelCollectionSize; ++i_BC) { 
    BasicClusterRef theBasicClusterRef(theBarrelBCCollection, i_BC);    
    if (theBasicClusterRef.isNull()) continue;  
    if (ROOT::Math::VectorUtil::DeltaR(math::XYZPoint(myECALHitPoint), (*theBasicClusterRef).position()) <= 0.7) {
      TLorentzVector myCluster((*theBasicClusterRef).position().x(),
			       (*theBasicClusterRef).position().y(),
			       (*theBasicClusterRef).position().z(),
			       (*theBasicClusterRef).energy());
      jet->clusters.push_back(myCluster);
    }
  }
  // Loop over endcap ECAL clusters
  unsigned int myEndcapCollectionSize = theEndcapBCCollection->size();
  for(unsigned int i_BC=0; i_BC != myEndcapCollectionSize; ++i_BC) { 
    BasicClusterRef theBasicClusterRef(theEndcapBCCollection, i_BC);    
    if (theBasicClusterRef.isNull()) continue;  
    if (ROOT::Math::VectorUtil::DeltaR(math::XYZPoint(myECALHitPoint), (*theBasicClusterRef).position()) <= 0.7) {
      TLorentzVector myCluster((*theBasicClusterRef).position().x(),
			       (*theBasicClusterRef).position().y(),
			       (*theBasicClusterRef).position().z(),
			       (*theBasicClusterRef).energy());
      jet->clusters.push_back(myCluster);
    }
  }
}
