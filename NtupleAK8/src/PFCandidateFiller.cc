/*
 * PFCandidateFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleAK8/interface/PFCandidateFiller.h"
#include "DeepNTuples/NtupleCommons/interface/sorting_modules.h"

namespace deepntuples {

void PFCandidateFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void PFCandidateFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(svToken_, SVs);
}

void PFCandidateFiller::book() {

  data.add<int>("n_pfcands", 0);
  data.add<float>("npfcands", 0);

  // basic kinematics, valid for both charged and neutral
  data.addMulti<float>("pfcand_ptrel");
  data.addMulti<float>("pfcand_erel");
  data.addMulti<float>("pfcand_phirel");
  data.addMulti<float>("pfcand_etarel");
  data.addMulti<float>("pfcand_deltaR");
  data.addMulti<float>("pfcand_puppiw");
  data.addMulti<float>("pfcand_mass");

  data.addMulti<float>("pfcand_drminsv");
  data.addMulti<float>("pfcand_drsubjet1");
  data.addMulti<float>("pfcand_drsubjet2");

  data.addMulti<float>("pfcand_charge");
  data.addMulti<float>("pfcand_isMu");
  data.addMulti<float>("pfcand_isEl");
  data.addMulti<float>("pfcand_isGamma");
  data.addMulti<float>("pfcand_isChargedHad");
  data.addMulti<float>("pfcand_isNeutralHad");

  // for neutral
  data.addMulti<float>("pfcand_hcalFrac");

  // for charged
  data.addMulti<float>("pfcand_VTX_ass");
  data.addMulti<float>("pfcand_fromPV");
  data.addMulti<float>("pfcand_lostInnerHits");

  // impact parameters
  data.addMulti<float>("pfcand_dz");
  data.addMulti<float>("pfcand_dzsig");
  data.addMulti<float>("pfcand_dxy");
  data.addMulti<float>("pfcand_dxysig");

}

bool PFCandidateFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  std::vector<const pat::PackedCandidate*> neutralPFCands;
  std::unordered_map<const pat::PackedCandidate*, double> drMinSvMap;
  //std::unordered_map<const pat::PackedCandidate*, double> ptRelMap;

  std::vector<sorting::sortingClass<size_t>> sortedneutral;

  const float jet_uncorr_pt=jet.correctedJet("Uncorrected").pt();

  unsigned int i = 0;
  for (const auto * pfcand : jet_helper.getJetConstituents()){
    if (!pfcand) continue;
    if (pfcand->pt() < minPt_) continue;
    if (pfcand->charge() == 0) {
      neutralPFCands.push_back(pfcand);
      drMinSvMap[pfcand];
      double minDR = 0.8;
      for (const auto &sv : *SVs){
        double dr = reco::deltaR(*pfcand, sv);
        if (dr < minDR) minDR = dr;
      }
      drMinSvMap[pfcand] = minDR;
      sortedneutral.push_back(sorting::sortingClass<size_t>
			      (i, -1,
			       -minDR, pfcand->pt()/jet_uncorr_pt));
      i++;
    }
  }

  // sort by ABCInv
  std::sort(sortedneutral.begin(),sortedneutral.end(),sorting::sortingClass<size_t>::compareByABCInv);
  std::vector<size_t> sortedneutralindices;
  sortedneutralindices=sorting::invertSortingVector(sortedneutral);

  data.fill<int>("n_pfcands", neutralPFCands.size());
  data.fill<float>("npfcands", neutralPFCands.size());

  float etasign = jet.eta()>0 ? 1 : -1;


  for (i = 0; i < neutralPFCands.size(); i++){

    const auto *pfcand = neutralPFCands.at(sortedneutral.at(i).get());    
    //  const auto *pfcand = neutralPFCands.at(i);

    if (pfcand->pt() < minPt_) continue;
    // basic kinematics, valid for both charged and neutral
    data.fillMulti<float>("pfcand_ptrel", pfcand->pt()/jet.pt());
    data.fillMulti<float>("pfcand_erel", pfcand->energy()/jet.energy());
    data.fillMulti<float>("pfcand_phirel", reco::deltaPhi(*pfcand, jet));
    data.fillMulti<float>("pfcand_etarel", etasign * (pfcand->eta() - jet.eta()));
    data.fillMulti<float>("pfcand_deltaR", reco::deltaR(*pfcand, jet));
    data.fillMulti<float>("pfcand_puppiw", pfcand->puppiWeight());
    data.fillMulti<float>("pfcand_mass", pfcand->mass());

    double minDR = 0.8;
    for (const auto &sv : *SVs){
      double dr = reco::deltaR(*pfcand, sv);
      if (dr < minDR) minDR = dr;
    }
    data.fillMulti<float>("pfcand_drminsv", catchInfsAndBound(drMinSvMap.at(pfcand),0,-0.8,0,-0.8));

    const auto& subjets = jet_helper.getSubJets();
    data.fillMulti<float>("pfcand_drsubjet1", subjets.size()>0 ? reco::deltaR(*pfcand, *subjets.at(0)) : -1);
    data.fillMulti<float>("pfcand_drsubjet2", subjets.size()>1 ? reco::deltaR(*pfcand, *subjets.at(1)) : -1);

    data.fillMulti<float>("pfcand_charge", pfcand->charge());
    data.fillMulti<float>("pfcand_isEl", std::abs(pfcand->pdgId())==11);
    data.fillMulti<float>("pfcand_isMu", std::abs(pfcand->pdgId())==13);
    data.fillMulti<float>("pfcand_isGamma", std::abs(pfcand->pdgId())==22);
    data.fillMulti<float>("pfcand_isChargedHad", std::abs(pfcand->pdgId())==211);
    data.fillMulti<float>("pfcand_isNeutralHad", std::abs(pfcand->pdgId())==130);

    // for neutral
    data.fillMulti<float>("pfcand_hcalFrac", pfcand->hcalFraction());

    // for charged
    data.fillMulti<float>("pfcand_VTX_ass", pfcand->pvAssociationQuality());
    data.fillMulti<float>("pfcand_fromPV", pfcand->fromPV());
    data.fillMulti<float>("pfcand_lostInnerHits", pfcand->lostInnerHits());

    // impact parameters
    data.fillMulti<float>("pfcand_dz", catchInfs(pfcand->dz()));
    data.fillMulti<float>("pfcand_dxy", catchInfs(pfcand->dxy()));
    if( pfcand->hasTrackDetails() ) {
      data.fillMulti<float>("pfcand_dzsig", catchInfsAndBound(pfcand->dz()/pfcand->dzError(),0,-2000,2000));
      data.fillMulti<float>("pfcand_dxysig", catchInfsAndBound(pfcand->dxy()/pfcand->dxyError(),0,-2000,2000));
    }
    else {
      data.fillMulti<float>("pfcand_dzsig", -1);
      data.fillMulti<float>("pfcand_dxysig", -1);
    }

  }

  return true;
}

} /* namespace deepntuples */
