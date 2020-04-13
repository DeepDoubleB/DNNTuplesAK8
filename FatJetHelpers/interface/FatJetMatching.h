/*
 * FatJetMatching.hh
 *
 *  Created on: Feb 1, 2017
 *      Author: hqu
 */

#ifndef FATJETHELPERS_INTERFACE_FATJETMATCHING_H_
#define FATJETHELPERS_INTERFACE_FATJETMATCHING_H_

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include <unordered_set>
#include <utility>

namespace deepntuples {

namespace ParticleID{
enum PdgId { p_unknown, p_d, p_u, p_s, p_c, p_b, p_t, p_bprime, p_tprime,
  p_eminus = 11, p_nu_e, p_muminus, p_nu_mu, p_tauminus, p_nu_tau,
  p_tauprimeminus, p_nu_tauprime, p_g = 21, p_gamma, p_Z0,
  p_Wplus, p_h0, p_Zprime0 = 32, p_Zpprime0, p_Wprimeplus, p_H0,
  p_A0, p_Hplus, p_G = 39, p_R0 = 41, p_H30 = 45, p_A20 = 46, p_Phi =55,
  p_LQ, p_cluster = 91, p_string,
  p_pi0 = 111, p_rho0 = 113, p_klong = 130, p_piplus = 211, p_rhoplus = 213, p_eta = 221, p_omega = 223,
  p_kshort = 310, p_k0, p_kstar0 = 313, p_kplus = 321, p_kstarplus = 323, p_phi = 333,
  p_dplus = 411, p_d0 = 421, p_dsplus = 431, p_b0 =511, p_bplus = 521,
  p_bs0 = 531, p_bcplus = 541,
  p_neutron = 2112, p_proton = 2212,
  p_sigmaminus = 3112, p_lambda0 = 3122,
  p_sigma0 = 3212, p_sigmaplus = 3222, p_ximinus = 3312, p_xi0 = 3322, p_omegaminus = 3334,
  p_sigmac0 = 4112, p_lambdacplus = 4122, p_xic0 = 4132,
  p_sigmacplus = 4212, p_sigmacpp = 4222, p_xicplus = 4232, p_omegac0 = 4332,
  p_sigmabminus = 5112, p_lambdab0 = 5122, p_xibminus = 5132, p_sigmab0 = 5212, p_sigmabplus = 5222,
  p_xib0 = 5232, p_omegabminus = 5332,
};
}

class FatJetMatching {
public:
  enum FatJetFlavor {
    Default = 0,
    Top = 1,
    W = 2,
    Z = 3,
    H = 4,
  };

  enum FatJetLabel {
    Invalid=0,
    Top_all=100, Top_bcq, Top_bqq, Top_bc, Top_bq, Top_bele, Top_bmu, Top_btau,
    W_all=200, W_cq, W_qq, W_cb, W_qb, W_ud, W_enu, W_munu, W_taunu,
    Z_all=300, Z_bb, Z_cc, Z_qq, Z_dd, Z_ee, Z_mumu, Z_tautau, Z_nunu,
    H_all=400, H_bb, H_cc, H_qq, H_qqqq, H_tautau, 
    H_WW_ud_ud, H_WW_ud_cs, H_WW_cs_ud, H_WW_cs_cs, 
    H_WW_ud_enu, H_WW_ud_munu, H_WW_ud_taunu, 
    H_WW_cs_enu, H_WW_cs_munu, H_WW_cs_taunu, 
    H_WW_enu_ud, H_WW_enu_cs, 
    H_WW_munu_ud, H_WW_munu_cs, 
    H_WW_taunu_ud, H_WW_taunu_cs,  
    H_WW_enu_enu, H_WW_enu_munu, H_WW_enu_taunu, 
    H_WW_munu_enu, H_WW_munu_munu, H_WW_munu_taunu, 
    H_WW_taunu_enu, H_WW_taunu_munu, H_WW_taunu_taunu, 
    H_ZZ_dd_dd, H_ZZ_dd_cc, H_ZZ_dd_bb, 
    H_ZZ_cc_dd, H_ZZ_cc_cc, H_ZZ_cc_bb, 
    H_ZZ_bb_dd, H_ZZ_bb_cc, H_ZZ_bb_bb, 
    H_ZZ_dd_ee, H_ZZ_dd_mumu, H_ZZ_dd_tautau, 
    H_ZZ_cc_ee, H_ZZ_cc_mumu, H_ZZ_cc_tautau, 
    H_ZZ_bb_ee, H_ZZ_bb_mumu, H_ZZ_bb_tautau, 
    H_ZZ_ee_dd, H_ZZ_ee_cc, H_ZZ_ee_bb, 
    H_ZZ_mumu_dd, H_ZZ_mumu_cc, H_ZZ_mumu_bb, 
    H_ZZ_tautau_dd, H_ZZ_tautau_cc, H_ZZ_tautau_bb, 
    H_ZZ_ee_ee, H_ZZ_ee_mumu, H_ZZ_ee_tautau, 
    H_ZZ_mumu_ee, H_ZZ_mumu_mumu, H_ZZ_mumu_tautau, 
    H_ZZ_tautau_ee, H_ZZ_tautau_mumu, H_ZZ_tautau_tautau, 
    H_ZZ_dd_nunu, H_ZZ_cc_nunu, H_ZZ_bb_nunu, 
    H_ZZ_nunu_dd, H_ZZ_nunu_cc, H_ZZ_nunu_bb, 
    H_ZZ_ee_nunu, H_ZZ_mumu_nunu, H_ZZ_tautau_nunu, 
    H_ZZ_nunu_ee, H_ZZ_nunu_mumu, H_ZZ_nunu_tautau, 
    H_ZZ_nunu_nunu,
    QCD_all=500, QCD_bb, QCD_cc, QCD_b, QCD_c, QCD_others,
  };

public:
  FatJetMatching() {}
  FatJetMatching(double jet_R, bool matchQuarks) : jetR_(jet_R), requiresQuarksContained_(matchQuarks) {}

  virtual ~FatJetMatching() {}

  std::pair<FatJetFlavor, const reco::GenParticle* > flavorJMAR(const pat::Jet *jet, const reco::GenParticleCollection& genParticles, double genRadius = 0.6);

  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > flavorLabel(const pat::Jet *jet, const reco::GenParticleCollection& genParticles, double distR);

private:
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > top_label(const pat::Jet *jet, const reco::GenParticle *parton, double distR);
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > w_label(const pat::Jet *jet, const reco::GenParticle *parton, double distR);
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > z_label(const pat::Jet *jet, const reco::GenParticle *parton, double distR);
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > higgs_label(const pat::Jet *jet, const reco::GenParticle *parton, double distR);
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > qcd_label(const pat::Jet *jet, const reco::GenParticleCollection& genParticles, double distR);


private:
  void printGenInfoHeader() const;
  void printGenParticleInfo(const reco::GenParticle* genParticle, const int idx) const;
  const reco::GenParticle* getFinal(const reco::GenParticle* particle);
  bool isHadronic(const reco::GenParticle* particle) const;
  std::vector<const reco::GenParticle*> getDaughterQuarks(const reco::GenParticle* particle);
  FatJetLabel getVLabel(const reco::GenParticle* particle) const;
  FatJetLabel getHVVLabel(const reco::GenParticle* daughter1, const reco::GenParticle* daughter2) const; 
  template <typename T>
  double maxDeltaRToDaughterQuarks(const T *center, const reco::GenParticle* mother) const {
    // mother particle needs to be the final version before decay
    double maxDeltaR = -1;
    for (const auto &q : mother->daughterRefVector()){
      if (std::abs(q->pdgId()) > ParticleID::p_b) continue;
      double deltaR = reco::deltaR(q->p4(), center->p4());
      if (deltaR > maxDeltaR) maxDeltaR = deltaR;
    }
    return maxDeltaR > 0 ? maxDeltaR : 1e9;
  }

private:
  double jetR_ = 0.8;
  bool   requiresQuarksContained_ = true;

  bool debug_ = false;
  std::unordered_set<const reco::GenParticle*> processed_;


};

}

#endif /* FATJETHELPERS_INTERFACE_FATJETMATCHING_H_ */
