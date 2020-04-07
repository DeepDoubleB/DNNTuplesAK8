/*
 * PFCompleteFiller.cc
 *
 *  Created on: Sep 25, 2017
 *      Author: hqu
 */

#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DeepNTuples/Ntupler/interface/PFCompleteFiller.h"

namespace deepntuples {

void PFCompleteFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void PFCompleteFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
}

void PFCompleteFiller::book() {

  data.add<int>("n_pfcands", 0);
  data.add<float>("npfcands", 0);

  // basic kinematics
  data.addMulti<float>("pfcand_ptrel");
  data.addMulti<float>("pfcand_erel");
  data.addMulti<float>("pfcand_pt");
  data.addMulti<float>("pfcand_mass");
  data.addMulti<float>("pfcand_ptrel_log");
  data.addMulti<float>("pfcand_erel_log");
  data.addMulti<float>("pfcand_pt_log");
  data.addMulti<float>("pfcand_e_log");

  // no puppi scaled
  data.addMulti<float>("pfcand_ptrel_nopuppi");
  data.addMulti<float>("pfcand_erel_nopuppi");
  data.addMulti<float>("pfcand_pt_nopuppi");
  data.addMulti<float>("pfcand_mass_nopuppi");
  data.addMulti<float>("pfcand_ptrel_log_nopuppi");
  data.addMulti<float>("pfcand_erel_log_nopuppi");
  data.addMulti<float>("pfcand_pt_log_nopuppi");
  data.addMulti<float>("pfcand_e_log_nopuppi");


  data.addMulti<float>("pfcand_phirel");
  data.addMulti<float>("pfcand_etarel");
  data.addMulti<float>("pfcand_deltaR");
  data.addMulti<float>("pfcand_puppiw");
  data.addMulti<float>("pfcand_abseta");

  data.addMulti<float>("pfcand_drminsv");
  data.addMulti<float>("pfcand_drminsvin"); // restricted to within the jet cone
  data.addMulti<float>("pfcand_drsubjet1"); // deprecated
  data.addMulti<float>("pfcand_drsubjet2"); // deprecated

  // use uncorrected pT to order the two subjets
  data.addMulti<float>("pfcand_dr_uncorrsj1");
  data.addMulti<float>("pfcand_dr_uncorrsj2");

  data.addMulti<float>("pfcand_charge");
  data.addMulti<float>("pfcand_isMu");
  data.addMulti<float>("pfcand_isEl");
  data.addMulti<float>("pfcand_isChargedHad");
  data.addMulti<float>("pfcand_isGamma");
  data.addMulti<float>("pfcand_isNeutralHad");

  // for neutral
  data.addMulti<float>("pfcand_hcalFrac");
  data.addMulti<float>("pfcand_hcalFracCalib");

  // for charged
  data.addMulti<float>("pfcand_VTX_ass");
  data.addMulti<float>("pfcand_fromPV");
  data.addMulti<float>("pfcand_lostInnerHits");
  data.addMulti<float>("pfcand_trackHighPurity");

  // impact parameters
  data.addMulti<float>("pfcand_dz");
  data.addMulti<float>("pfcand_dzsig");
  data.addMulti<float>("pfcand_dxy");
  data.addMulti<float>("pfcand_dxysig");

  // track quality
  data.addMulti<float>("pfcand_normchi2");
  data.addMulti<float>("pfcand_quality");

  // track covariance
  data.addMulti<float>("pfcand_dptdpt");
  data.addMulti<float>("pfcand_detadeta");
  data.addMulti<float>("pfcand_dphidphi");
  data.addMulti<float>("pfcand_dxydxy");
  data.addMulti<float>("pfcand_dzdz");
  data.addMulti<float>("pfcand_dxydz");
  data.addMulti<float>("pfcand_dphidxy");
  data.addMulti<float>("pfcand_dlambdadz");

  // track btag info
  data.addMulti<float>("pfcand_btagMomentum");
  data.addMulti<float>("pfcand_btagEta");
  data.addMulti<float>("pfcand_btagEtaRel");
  data.addMulti<float>("pfcand_btagPtRel");
  data.addMulti<float>("pfcand_btagPPar");
  data.addMulti<float>("pfcand_btagDeltaR");
  data.addMulti<float>("pfcand_btagPtRatio");
  data.addMulti<float>("pfcand_btagPParRatio");
  data.addMulti<float>("pfcand_btagSip2dVal");
  data.addMulti<float>("pfcand_btagSip2dSig");
  data.addMulti<float>("pfcand_btagSip3dVal");
  data.addMulti<float>("pfcand_btagSip3dSig");
  data.addMulti<float>("pfcand_btagJetDistVal");
  data.addMulti<float>("pfcand_btagJetDistSig"); // always gives 0?

}

bool PFCompleteFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  const auto& pfCands = jet_helper.getJetConstituents();

  data.fill<int>("n_pfcands", pfCands.size());
  data.fill<float>("npfcands", pfCands.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  for (const auto& cand : pfCands){

    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    
    if (packed_cand->pt() < minPt_) continue;
    // basic kinematics, valid for both charged and neutral
    {
      // puppi weighted
      auto puppiP4 = jet_helper.getPuppiWeight(cand) * cand->p4();
      data.fillMulti<float>("pfcand_ptrel", puppiP4.pt()/jet.pt());
      data.fillMulti<float>("pfcand_erel", puppiP4.energy()/jet.energy());
      data.fillMulti<float>("pfcand_pt", puppiP4.pt());
      data.fillMulti<float>("pfcand_mass", puppiP4.mass());
      data.fillMulti<float>("pfcand_ptrel_log", catchInfs(std::log(puppiP4.pt()/jet.pt()), -99));
      data.fillMulti<float>("pfcand_erel_log", catchInfs(std::log(puppiP4.energy()/jet.energy()), -99));
      data.fillMulti<float>("pfcand_pt_log", catchInfs(std::log(puppiP4.pt()), -99));
      data.fillMulti<float>("pfcand_e_log", catchInfs(std::log(puppiP4.energy()), -99));
    }

    // not puppi weighted
    data.fillMulti<float>("pfcand_ptrel_nopuppi", packed_cand->pt()/jet.pt());
    data.fillMulti<float>("pfcand_erel_nopuppi", packed_cand->energy()/jet.energy());
    data.fillMulti<float>("pfcand_pt_nopuppi", packed_cand->pt());
    data.fillMulti<float>("pfcand_mass_nopuppi", packed_cand->mass());
    data.fillMulti<float>("pfcand_ptrel_log_nopuppi", catchInfs(std::log(packed_cand->pt()/jet.pt()), -99));
    data.fillMulti<float>("pfcand_erel_log_nopuppi", catchInfs(std::log(packed_cand->energy()/jet.energy()), -99));
    data.fillMulti<float>("pfcand_pt_log_nopuppi", catchInfs(std::log(packed_cand->pt()), -99));
    data.fillMulti<float>("pfcand_e_log_nopuppi", catchInfs(std::log(packed_cand->energy()), -99));

    data.fillMulti<float>("pfcand_phirel", reco::deltaPhi(*packed_cand, jet));
    data.fillMulti<float>("pfcand_etarel", etasign * (packed_cand->eta() - jet.eta()));
    data.fillMulti<float>("pfcand_deltaR", reco::deltaR(*packed_cand, jet));
    data.fillMulti<float>("pfcand_abseta", std::abs(packed_cand->eta()));

    data.fillMulti<float>("pfcand_puppiw", jet_helper.getPuppiWeight(cand));

    double minDR = 999;
    double minDRin = 2.*jetR_;
    for (const auto &sv : *SVs){
      double dr = reco::deltaR(*packed_cand, sv);
      if (dr < minDR) minDR = dr;
      if (dr < minDRin && reco::deltaR(jet, sv) < jetR_) minDRin = dr;
    }
    data.fillMulti<float>("pfcand_drminsv", minDR==999 ? -1 : minDR);
    data.fillMulti<float>("pfcand_drminsvin", minDRin);

    { // deprecated
      const auto& subjets = jet_helper.getSubJets();
      data.fillMulti<float>("pfcand_drsubjet1", subjets.size()>0 ? reco::deltaR(*packed_cand, *subjets.at(0)) : -1);
      data.fillMulti<float>("pfcand_drsubjet2", subjets.size()>1 ? reco::deltaR(*packed_cand, *subjets.at(1)) : -1);
    }

    // use uncorrected pT to order the two subjets
    {
      const auto& subjets = jet_helper.getUncorrSubJets();
      data.fillMulti<float>("pfcand_dr_uncorrsj1", subjets.size()>0 ? reco::deltaR(*packed_cand, *subjets.at(0)) : -1);
      data.fillMulti<float>("pfcand_dr_uncorrsj2", subjets.size()>1 ? reco::deltaR(*packed_cand, *subjets.at(1)) : -1);
    }

    data.fillMulti<float>("pfcand_charge", packed_cand->charge());
    data.fillMulti<float>("pfcand_isEl", std::abs(packed_cand->pdgId())==11);
    data.fillMulti<float>("pfcand_isMu", std::abs(packed_cand->pdgId())==13);
    data.fillMulti<float>("pfcand_isChargedHad", std::abs(packed_cand->pdgId())==211);
    data.fillMulti<float>("pfcand_isGamma", std::abs(packed_cand->pdgId())==22);
    data.fillMulti<float>("pfcand_isNeutralHad", std::abs(packed_cand->pdgId())==130);

    // for neutral
    float hcal_fraction = 0.;
    if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130) {
      hcal_fraction = packed_cand->hcalFraction();
    } else if (packed_cand->isIsolatedChargedHadron()) {
      hcal_fraction = packed_cand->rawHcalFraction();
    }
    data.fillMulti<float>("pfcand_hcalFrac", hcal_fraction);
    data.fillMulti<float>("pfcand_hcalFracCalib", packed_cand->hcalFraction());

    // for charged
    data.fillMulti<float>("pfcand_VTX_ass", packed_cand->pvAssociationQuality());
    data.fillMulti<float>("pfcand_fromPV", packed_cand->fromPV());
    data.fillMulti<float>("pfcand_lostInnerHits", packed_cand->lostInnerHits());
    data.fillMulti<float>("pfcand_trackHighPurity", packed_cand->trackHighPurity());

    // impact parameters
    data.fillMulti<float>("pfcand_dz", catchInfs(packed_cand->dz()));
    data.fillMulti<float>("pfcand_dzsig", packed_cand->bestTrack() ? catchInfs(packed_cand->dz()/packed_cand->dzError()) : 0);
    data.fillMulti<float>("pfcand_dxy", catchInfs(packed_cand->dxy()));
    data.fillMulti<float>("pfcand_dxysig", packed_cand->bestTrack() ? catchInfs(packed_cand->dxy()/packed_cand->dxyError()) : 0);

    if (packed_cand->bestTrack()){
      const auto *trk = packed_cand->bestTrack();
      data.fillMulti<float>("pfcand_normchi2", catchInfs(trk->normalizedChi2()));
      data.fillMulti<float>("pfcand_quality", trk->qualityMask());

      // track covariance
      auto cov = [&](unsigned i, unsigned j) {
        return catchInfs(trk->covariance(i, j));
      };
      data.fillMulti<float>("pfcand_dptdpt", cov(0,0));
      data.fillMulti<float>("pfcand_detadeta", cov(1,1));
      data.fillMulti<float>("pfcand_dphidphi", cov(2,2));
      data.fillMulti<float>("pfcand_dxydxy", cov(3,3));
      data.fillMulti<float>("pfcand_dzdz", cov(4,4));
      data.fillMulti<float>("pfcand_dxydz", cov(3,4));
      data.fillMulti<float>("pfcand_dphidxy", cov(2,3));
      data.fillMulti<float>("pfcand_dlambdadz", cov(1,4));
    }else{
      data.fillMulti<float>("pfcand_normchi2", 999);
      data.fillMulti<float>("pfcand_quality", 0);

      data.fillMulti<float>("pfcand_dptdpt", 0);
      data.fillMulti<float>("pfcand_detadeta", 0);
      data.fillMulti<float>("pfcand_dphidphi", 0);
      data.fillMulti<float>("pfcand_dxydxy", 0);
      data.fillMulti<float>("pfcand_dzdz", 0);
      data.fillMulti<float>("pfcand_dxydz", 0);
      data.fillMulti<float>("pfcand_dphidxy", 0);
      data.fillMulti<float>("pfcand_dlambdadz", 0);
    }

    // build track info map
    TrackInfoBuilder trkinfo;
    trkinfo.buildTrackInfo(builder_, *packed_cand, jet, vertices->at(0));

    data.fillMulti<float>("pfcand_btagMomentum", trkinfo.getTrackMomentum());
    data.fillMulti<float>("pfcand_btagEta", trkinfo.getTrackEta());
    data.fillMulti<float>("pfcand_btagEtaRel", trkinfo.getTrackEtaRel());
    data.fillMulti<float>("pfcand_btagPtRel", trkinfo.getTrackPtRel());
    data.fillMulti<float>("pfcand_btagPPar", trkinfo.getTrackPPar());
    data.fillMulti<float>("pfcand_btagDeltaR", trkinfo.getTrackDeltaR());
    data.fillMulti<float>("pfcand_btagPtRatio", trkinfo.getTrackPtRatio());
    data.fillMulti<float>("pfcand_btagPParRatio", trkinfo.getTrackPParRatio());
    data.fillMulti<float>("pfcand_btagSip2dVal", trkinfo.getTrackSip2dVal());
    data.fillMulti<float>("pfcand_btagSip2dSig", trkinfo.getTrackSip2dSig());
    data.fillMulti<float>("pfcand_btagSip3dVal", trkinfo.getTrackSip3dVal());
    data.fillMulti<float>("pfcand_btagSip3dSig", trkinfo.getTrackSip3dSig());
    data.fillMulti<float>("pfcand_btagJetDistVal", trkinfo.getTrackJetDistVal());
    data.fillMulti<float>("pfcand_btagJetDistSig", trkinfo.getTrackJetDistSig());


  }


  return true;
}

} /* namespace deepntuples */
