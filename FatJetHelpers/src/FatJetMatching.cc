/*
 * FatJetMatching.cc
 *
 *  Created on: Feb 1, 2017
 *      Author: hqu
 */

#include "DeepNTuples/FatJetHelpers/interface/FatJetMatching.h"

#include <unordered_set>
#include "TString.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace deepntuples;

std::pair<FatJetMatching::FatJetFlavor, const reco::GenParticle*> FatJetMatching::flavorJMAR(const pat::Jet* jet,
    const reco::GenParticleCollection& genParticles, double genRadius) {

  processed_.clear();

  if (debug_) {
    std::cout << "\n=======\nJet (energy, pT, eta, phi) = "
        << jet->energy() << ", " << jet->pt() << ", " << jet->eta() << ", " << jet->phi()
        << std::endl << std::endl;
    printGenInfoHeader();
    for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
      printGenParticleInfo(&genParticles[ipart], ipart);
    }
  }

  for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
    const auto *gp = &genParticles[ipart];

    if (processed_.count(gp)) continue;
    processed_.insert(gp);

    auto pdgid = std::abs(gp->pdgId());
    if (pdgid == ParticleID::p_t){
      // top
      auto top = getFinal(gp);
      // find the W and test if it's hadronic
      const reco::GenParticle *w_from_top = nullptr, *b_from_top = nullptr;
      for (const auto &dau : top->daughterRefVector()){
        if (std::abs(dau->pdgId()) == ParticleID::p_Wplus){
          w_from_top = getFinal(&(*dau));
        }else if (std::abs(dau->pdgId()) <= ParticleID::p_b){
          // ! use <= p_b ! -- can also have charms etc.
          // for quarks use the first one in the decay chain
          b_from_top = dynamic_cast<const reco::GenParticle*>(&(*dau));
        }
      }
      if (!w_from_top || !b_from_top) throw std::logic_error("[FatJetMatching::flavor] Cannot find b or W from top decay: "+std::to_string(ipart));
      if (isHadronic(w_from_top)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "top: "; printGenParticleInfo(top, -1);
          cout << "b:   "; printGenParticleInfo(b_from_top, -1);
          cout << "W:   "; printGenParticleInfo(w_from_top, -1);
        }

        double dr_jet_top = reco::deltaR(jet->p4(), top->p4());
        double dr_top_wdaus = maxDeltaRToDaughterQuarks(top, w_from_top);
        double dr_top_b     = reco::deltaR(top->p4(), b_from_top->p4());
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, top)   : " << dr_jet_top << endl;
          cout << "deltaR(top, b)     : " << dr_top_b << endl;
          cout << "deltaR(top, w daus): " << dr_top_wdaus << endl;
        }
        // top
        if (dr_top_wdaus < genRadius && dr_top_b < genRadius && dr_jet_top < genRadius) return std::make_pair(FatJetFlavor::Top, top);

        double dr_jet_w = reco::deltaR(jet->p4(), w_from_top->p4());
        double dr_w = maxDeltaRToDaughterQuarks(w_from_top, w_from_top);
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, w)     : " << dr_jet_w << endl;
          cout << "deltaR(w, w_daus)  : " << dr_w << endl;
        }
        if (dr_w < genRadius && dr_jet_w < genRadius) return std::make_pair(FatJetFlavor::W, w_from_top);
      }
    }else if (pdgid == ParticleID::p_h0 || pdgid ==ParticleID::p_Phi) {
      // Higgs
      auto h = getFinal(gp);
      if (isHadronic(h)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "H:   "; printGenParticleInfo(h, -1);
        }
        double dr_jet_h = reco::deltaR(jet->p4(), h->p4());
        double dr_hdaus = maxDeltaRToDaughterQuarks(h, h); // only works for h->bb??
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, H)   : " << dr_jet_h << endl;
          cout << "deltaR(h, h daus): " << dr_hdaus << endl;
        }
        if (dr_hdaus < genRadius && dr_jet_h < genRadius) return std::make_pair(FatJetFlavor::H, h);
      }
    }else if (pdgid == ParticleID::p_Wplus){
      // W: not from top, or top not in jet cone
      auto w = getFinal(gp);
      if (isHadronic(w)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "W:   "; printGenParticleInfo(w, -1);
        }
        double dr_jet_w = reco::deltaR(jet->p4(), w->p4());
        double dr_wdaus = maxDeltaRToDaughterQuarks(w, w);
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, w)   : " << dr_jet_w << endl;
          cout << "deltaR(w, w daus): " << dr_wdaus << endl;
        }
        if (dr_wdaus < genRadius && dr_jet_w < genRadius) return std::make_pair(FatJetFlavor::W, w);

      }
    }else if (pdgid == ParticleID::p_Z0) {
      // Z
      auto z = getFinal(gp);
      if (isHadronic(z)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "Z:   "; printGenParticleInfo(z, -1);
        }
        double dr_jet_z = reco::deltaR(jet->p4(), z->p4());
        double dr_zdaus = maxDeltaRToDaughterQuarks(z, z);
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, Z)   : " << dr_jet_z << endl;
          cout << "deltaR(Z, Z daus): " << dr_zdaus << endl;
        }
        if (dr_zdaus < genRadius && dr_jet_z < genRadius) return std::make_pair(FatJetFlavor::Z, z);
      }
    }else {
      // ?
    }
  }

  if (genParticles.size() != processed_.size())
    throw std::logic_error("[FatJetMatching::flavor] Not all genParticles are processed!");

  const reco::GenParticle *parton = nullptr;
  double minDR = 999;
  for (const auto &gp : genParticles){
    if (gp.status() != 23) continue;
    auto pdgid = std::abs(gp.pdgId());
    if (!(pdgid<ParticleID::p_t || pdgid==ParticleID::p_g)) continue;
    auto dr = reco::deltaR(gp, *jet);
    if (dr<genRadius && dr<minDR){
      minDR = dr;
      parton = &gp;
    }
  }
  if (debug_){
    using namespace std;
    if (parton){
      cout << "parton"; printGenParticleInfo(parton, -1);
      cout << "dr(jet, parton): " << minDR << endl;
    }
  }

  return std::make_pair(FatJetFlavor::Default, parton);

}



std::pair<FatJetMatching::FatJetLabel, std::vector<const reco::GenParticle*> > FatJetMatching::flavorLabel(const pat::Jet* jet,
    const reco::GenParticleCollection& genParticles, double distR) {

  processed_.clear();

  if (debug_) {
    std::cout << "\n=======\nJet (energy, pT, eta, phi) = "
        << jet->energy() << ", " << jet->pt() << ", " << jet->eta() << ", " << jet->phi()
        << std::endl << std::endl;
    printGenInfoHeader();
    for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
      printGenParticleInfo(&genParticles[ipart], ipart);
    }
  }

  for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
    const auto *gp = &genParticles[ipart];

    if (processed_.count(gp)) continue;
    processed_.insert(gp);

    auto pdgid = std::abs(gp->pdgId());
    if (pdgid == ParticleID::p_t){
      auto result = top_label(jet, gp, distR);
      if (result.first != FatJetLabel::Invalid){
        return result;
      }
    }else if (pdgid == ParticleID::p_h0 || pdgid ==ParticleID::p_Phi){            
      auto result = higgs_label(jet, gp, distR);
      if (result.first != FatJetLabel::Invalid){
        return result;
      } else {
      }
    }else if (pdgid == ParticleID::p_Wplus){
      auto result = w_label(jet, gp, distR);
      if (result.first != FatJetLabel::Invalid){
        return result;
      }
    }else if (pdgid == ParticleID::p_Z0){
      auto result = z_label(jet, gp, distR);
      if (result.first != FatJetLabel::Invalid){
        return result;
      }
    }
  }

  if (genParticles.size() != processed_.size())
    throw std::logic_error("[FatJetMatching::flavor] Not all genParticles are processed!");

  return qcd_label(jet, genParticles, distR);

}


void FatJetMatching::printGenInfoHeader() const {
  using namespace std;
  cout    << right << setw(6) << "#" << " " << setw(10) << "pdgId"
      << "  " << "Chg" << "  " << setw(10) << "Mass" << "  " << setw(48) << " Momentum"
      << left << "  " << setw(10) << "Mothers" << " " << setw(30) << "Daughters" << endl;
}

void FatJetMatching::printGenParticleInfo(const reco::GenParticle* genParticle, const int idx) const {
  using namespace std;
  cout  << right << setw(3) << genParticle->status();
  cout  << right << setw(3) << idx << " " << setw(10) << genParticle->pdgId() << "  ";
  cout  << right << "  " << setw(3) << genParticle->charge() << "  " << TString::Format("%10.3g", genParticle->mass() < 1e-5 ? 0 : genParticle->mass());
  cout  << left << setw(50) << TString::Format("  (E=%6.4g pT=%6.4g eta=%7.3g phi=%7.3g)", genParticle->energy(), genParticle->pt(), genParticle->eta(), genParticle->phi());

  TString                     mothers;
  for (unsigned int iMom = 0; iMom < genParticle->numberOfMothers(); ++iMom) {
    if (mothers.Length())     mothers        += ",";
    mothers   += genParticle->motherRef(iMom).key();
  }
  cout << "  " << setw(10) << mothers;
  TString                     daughters;
  for (unsigned int iDau = 0; iDau < genParticle->numberOfDaughters(); ++iDau) {
    if (daughters.Length())   daughters      += ",";
    daughters += genParticle->daughterRef(iDau).key();
  }
  cout << " " << setw(30) << daughters << endl;
}

const reco::GenParticle* FatJetMatching::getFinal(const reco::GenParticle* particle) {
  // will mark intermediate particles as processed
  if (!particle) return nullptr;
  processed_.insert(particle);
  const reco::GenParticle *final = particle;

  while (final->numberOfDaughters()) {
    const reco::GenParticle *chain = nullptr;
    for (unsigned idau = 0; idau < final->numberOfDaughters(); ++idau){
      if (final->daughter(idau)->pdgId() == particle->pdgId()) {
        chain = dynamic_cast<const reco::GenParticle*>(final->daughter(idau));
        processed_.insert(chain);
        break;
      }
    }
    if (!chain) break;
    final = chain;
  }
  return final;
}

bool FatJetMatching::isHadronic(const reco::GenParticle* particle) const {
  // particle needs to be the final version before decay
  if (!particle) throw std::invalid_argument("[FatJetMatching::isHadronic()] Null particle!");
  for(const auto &dau : particle->daughterRefVector()){
    auto pdgid = std::abs(dau->pdgId());
    if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b) return true;
  }
  return false;
}

std::vector<const reco::GenParticle*> FatJetMatching::getDaughterQuarks(const reco::GenParticle* particle) {
  std::vector<const reco::GenParticle*> daughters;

  for (unsigned i=0; i<particle->numberOfDaughters(); ++i){
    const auto *dau = dynamic_cast<const reco::GenParticle*>(particle->daughter(i));
    auto pdgid = std::abs(dau->pdgId());
    if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b){
      daughters.push_back(dau);
    }
  }

  return daughters;
}

FatJetMatching::FatJetLabel FatJetMatching::getVLabel(const reco::GenParticle* particle) const {
  if (particle->pdgId() == ParticleID::p_Z0) {
    if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_b) return FatJetLabel::Z_bb;
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_c) return FatJetLabel::Z_cc;
    else if (abs(particle->daughter(0)->pdgId()) >= ParticleID::p_d && abs(particle->daughter(0)->pdgId()) <= ParticleID::p_s) return FatJetLabel::Z_dd;
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_eminus) return FatJetLabel::Z_ee;
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_muminus) return FatJetLabel::Z_mumu;
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_tauminus) return FatJetLabel::Z_tautau;
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_nu_e || 
	     abs(particle->daughter(0)->pdgId()) == ParticleID::p_nu_mu || 
	     abs(particle->daughter(0)->pdgId()) == ParticleID::p_nu_tau
	     ) return FatJetLabel::Z_nunu;
  } else {
    if ( (abs(particle->daughter(0)->pdgId()) == ParticleID::p_b && abs(particle->daughter(1)->pdgId()) == ParticleID::p_c)
	 || (abs(particle->daughter(0)->pdgId()) == ParticleID::p_c && abs(particle->daughter(1)->pdgId()) == ParticleID::p_b) ) return FatJetLabel::W_cb;	    
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_b || abs(particle->daughter(1)->pdgId()) == ParticleID::p_b) return FatJetLabel::W_qb;	    
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_c || abs(particle->daughter(1)->pdgId()) == ParticleID::p_c) return FatJetLabel::W_cq;	    
    else if (abs(particle->daughter(0)->pdgId()) >= ParticleID::p_d && abs(particle->daughter(0)->pdgId()) <= ParticleID::p_s) return FatJetLabel::W_ud;	    
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_eminus || abs(particle->daughter(1)->pdgId()) == ParticleID::p_eminus) return FatJetLabel::W_enu;	    
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_muminus || abs(particle->daughter(1)->pdgId()) == ParticleID::p_muminus) return FatJetLabel::W_munu;	    
    else if (abs(particle->daughter(0)->pdgId()) == ParticleID::p_tauminus || abs(particle->daughter(1)->pdgId()) == ParticleID::p_tauminus) return FatJetLabel::W_taunu;	          
  }
  return FatJetLabel::Invalid;
}

FatJetMatching::FatJetLabel FatJetMatching::getHVVLabel(const reco::GenParticle* daughter1, const reco::GenParticle* daughter2) const {
  FatJetLabel daughter1Label = getVLabel(daughter1);
  FatJetLabel daughter2Label = getVLabel(daughter2);
  if (daughter1Label == FatJetLabel::W_ud && daughter2Label == FatJetLabel::W_ud) return FatJetLabel::H_WW_ud_ud;
  else if ( daughter1Label == FatJetLabel::W_ud 
	    && ( daughter2Label == FatJetLabel::W_cq || daughter2Label == FatJetLabel::W_cb || daughter2Label == FatJetLabel::W_qb)	      
	    ) return FatJetLabel::H_WW_ud_cs;
  else if ( daughter1Label == FatJetLabel::W_ud  &&  daughter2Label == FatJetLabel::W_enu ) 
    return FatJetLabel::H_WW_ud_enu;
  else if ( daughter1Label == FatJetLabel::W_ud  &&  daughter2Label == FatJetLabel::W_munu ) 
    return FatJetLabel::H_WW_ud_munu;
  else if ( daughter1Label == FatJetLabel::W_ud  &&  daughter2Label == FatJetLabel::W_taunu ) 
    return FatJetLabel::H_WW_ud_taunu;
  else if ( ( daughter1Label == FatJetLabel::W_cq || daughter1Label == FatJetLabel::W_cb || daughter1Label == FatJetLabel::W_qb) 
	    && daughter2Label == FatJetLabel::W_ud 
	    ) return FatJetLabel::H_WW_cs_ud;
  else if ( ( daughter1Label == FatJetLabel::W_cq || daughter1Label == FatJetLabel::W_cb || daughter1Label == FatJetLabel::W_qb) 
	    && ( daughter2Label == FatJetLabel::W_cq || daughter2Label == FatJetLabel::W_cb || daughter2Label == FatJetLabel::W_qb) 
	    ) return FatJetLabel::H_WW_cs_cs;
  else if ( ( daughter1Label == FatJetLabel::W_cq || daughter1Label == FatJetLabel::W_cb || daughter1Label == FatJetLabel::W_qb) 
	    && daughter2Label == FatJetLabel::W_enu 
	    ) return FatJetLabel::H_WW_cs_enu;
  else if ( ( daughter1Label == FatJetLabel::W_cq || daughter1Label == FatJetLabel::W_cb || daughter1Label == FatJetLabel::W_qb) 
	      && daughter2Label == FatJetLabel::W_munu 
	    ) return FatJetLabel::H_WW_cs_munu;
  else if ( ( daughter1Label == FatJetLabel::W_cq || daughter1Label == FatJetLabel::W_cb || daughter1Label == FatJetLabel::W_qb) 
	    && daughter2Label == FatJetLabel::W_taunu 
	    ) return FatJetLabel::H_WW_cs_taunu;
  else if ( daughter1Label == FatJetLabel::W_enu && daughter2Label == FatJetLabel::W_ud 
	    ) return FatJetLabel::H_WW_enu_ud;
  else if ( daughter1Label == FatJetLabel::W_enu
	    && ( daughter2Label == FatJetLabel::W_cq || daughter2Label == FatJetLabel::W_cb || daughter2Label == FatJetLabel::W_qb) 
	    ) return FatJetLabel::H_WW_enu_cs;
  else if ( daughter1Label == FatJetLabel::W_enu && daughter2Label == FatJetLabel::W_enu 
	    ) return FatJetLabel::H_WW_enu_enu;
  else if ( daughter1Label == FatJetLabel::W_enu && daughter2Label == FatJetLabel::W_munu 
	    ) return FatJetLabel::H_WW_enu_munu;
  else if ( daughter1Label == FatJetLabel::W_enu && daughter2Label == FatJetLabel::W_taunu 
	    ) return FatJetLabel::H_WW_enu_taunu;	      
  else if ( daughter1Label == FatJetLabel::W_munu && daughter2Label == FatJetLabel::W_ud 
	    ) return FatJetLabel::H_WW_munu_ud;
  else if ( daughter1Label == FatJetLabel::W_munu
	    && ( daughter2Label == FatJetLabel::W_cq || daughter2Label == FatJetLabel::W_cb || daughter2Label == FatJetLabel::W_qb) 
	    ) return FatJetLabel::H_WW_munu_cs;
  else if ( daughter1Label == FatJetLabel::W_munu && daughter2Label == FatJetLabel::W_enu 
	    ) return FatJetLabel::H_WW_munu_enu;
  else if ( daughter1Label == FatJetLabel::W_munu && daughter2Label == FatJetLabel::W_munu 
	    ) return FatJetLabel::H_WW_munu_munu;
  else if ( daughter1Label == FatJetLabel::W_munu && daughter2Label == FatJetLabel::W_taunu 
	    ) return FatJetLabel::H_WW_munu_taunu;	      
  else if ( daughter1Label == FatJetLabel::W_taunu && daughter2Label == FatJetLabel::W_ud 
	    ) return FatJetLabel::H_WW_taunu_ud;
  else if ( daughter1Label == FatJetLabel::W_taunu
	    && ( daughter2Label == FatJetLabel::W_cq || daughter2Label == FatJetLabel::W_cb || daughter2Label == FatJetLabel::W_qb) 
	    ) return FatJetLabel::H_WW_taunu_cs;
  else if ( daughter1Label == FatJetLabel::W_taunu && daughter2Label == FatJetLabel::W_enu 
	    ) return FatJetLabel::H_WW_taunu_enu;
  else if ( daughter1Label == FatJetLabel::W_taunu && daughter2Label == FatJetLabel::W_munu 
	    ) return FatJetLabel::H_WW_taunu_munu;
  else if ( daughter1Label == FatJetLabel::W_taunu && daughter2Label == FatJetLabel::W_taunu 
	    ) return FatJetLabel::H_WW_taunu_taunu;
  else if (daughter1Label == FatJetLabel::Z_dd && daughter2Label == FatJetLabel::Z_dd ) return FatJetLabel::H_ZZ_dd_dd;
  else if (daughter1Label == FatJetLabel::Z_dd && daughter2Label == FatJetLabel::Z_cc ) return FatJetLabel::H_ZZ_dd_cc;
  else if (daughter1Label == FatJetLabel::Z_dd && daughter2Label == FatJetLabel::Z_bb ) return FatJetLabel::H_ZZ_dd_bb;
  else if (daughter1Label == FatJetLabel::Z_dd && daughter2Label == FatJetLabel::Z_ee ) return FatJetLabel::H_ZZ_dd_ee;
  else if (daughter1Label == FatJetLabel::Z_dd && daughter2Label == FatJetLabel::Z_mumu ) return FatJetLabel::H_ZZ_dd_mumu;
  else if (daughter1Label == FatJetLabel::Z_dd && daughter2Label == FatJetLabel::Z_tautau ) return FatJetLabel::H_ZZ_dd_tautau;
  else if (daughter1Label == FatJetLabel::Z_dd && daughter2Label == FatJetLabel::Z_nunu ) return FatJetLabel::H_ZZ_dd_nunu;
  else if (daughter1Label == FatJetLabel::Z_cc && daughter2Label == FatJetLabel::Z_dd ) return FatJetLabel::H_ZZ_cc_dd;
  else if (daughter1Label == FatJetLabel::Z_cc && daughter2Label == FatJetLabel::Z_cc ) return FatJetLabel::H_ZZ_cc_cc;
  else if (daughter1Label == FatJetLabel::Z_cc && daughter2Label == FatJetLabel::Z_bb ) return FatJetLabel::H_ZZ_cc_bb;
  else if (daughter1Label == FatJetLabel::Z_cc && daughter2Label == FatJetLabel::Z_ee ) return FatJetLabel::H_ZZ_cc_ee;
  else if (daughter1Label == FatJetLabel::Z_cc && daughter2Label == FatJetLabel::Z_mumu ) return FatJetLabel::H_ZZ_cc_mumu;
  else if (daughter1Label == FatJetLabel::Z_cc && daughter2Label == FatJetLabel::Z_tautau ) return FatJetLabel::H_ZZ_cc_tautau;
  else if (daughter1Label == FatJetLabel::Z_cc && daughter2Label == FatJetLabel::Z_nunu ) return FatJetLabel::H_ZZ_cc_nunu;	    
  else if (daughter1Label == FatJetLabel::Z_bb && daughter2Label == FatJetLabel::Z_dd ) return FatJetLabel::H_ZZ_bb_dd;
  else if (daughter1Label == FatJetLabel::Z_bb && daughter2Label == FatJetLabel::Z_cc ) return FatJetLabel::H_ZZ_bb_cc;
  else if (daughter1Label == FatJetLabel::Z_bb && daughter2Label == FatJetLabel::Z_bb ) return FatJetLabel::H_ZZ_bb_bb;
  else if (daughter1Label == FatJetLabel::Z_bb && daughter2Label == FatJetLabel::Z_ee ) return FatJetLabel::H_ZZ_bb_ee;
  else if (daughter1Label == FatJetLabel::Z_bb && daughter2Label == FatJetLabel::Z_mumu ) return FatJetLabel::H_ZZ_bb_mumu;
  else if (daughter1Label == FatJetLabel::Z_bb && daughter2Label == FatJetLabel::Z_tautau ) return FatJetLabel::H_ZZ_bb_tautau;
  else if (daughter1Label == FatJetLabel::Z_bb && daughter2Label == FatJetLabel::Z_nunu ) return FatJetLabel::H_ZZ_bb_nunu;	     
  else if (daughter1Label == FatJetLabel::Z_ee && daughter2Label == FatJetLabel::Z_dd ) return FatJetLabel::H_ZZ_ee_dd;
  else if (daughter1Label == FatJetLabel::Z_ee && daughter2Label == FatJetLabel::Z_cc ) return FatJetLabel::H_ZZ_ee_cc;
  else if (daughter1Label == FatJetLabel::Z_ee && daughter2Label == FatJetLabel::Z_bb ) return FatJetLabel::H_ZZ_ee_bb;
  else if (daughter1Label == FatJetLabel::Z_ee && daughter2Label == FatJetLabel::Z_ee ) return FatJetLabel::H_ZZ_ee_ee;
  else if (daughter1Label == FatJetLabel::Z_ee && daughter2Label == FatJetLabel::Z_mumu ) return FatJetLabel::H_ZZ_ee_mumu;
  else if (daughter1Label == FatJetLabel::Z_ee && daughter2Label == FatJetLabel::Z_tautau ) return FatJetLabel::H_ZZ_ee_tautau;
  else if (daughter1Label == FatJetLabel::Z_ee && daughter2Label == FatJetLabel::Z_nunu ) return FatJetLabel::H_ZZ_ee_nunu;	     
  else if (daughter1Label == FatJetLabel::Z_mumu && daughter2Label == FatJetLabel::Z_dd ) return FatJetLabel::H_ZZ_mumu_dd;
  else if (daughter1Label == FatJetLabel::Z_mumu && daughter2Label == FatJetLabel::Z_cc ) return FatJetLabel::H_ZZ_mumu_cc;
  else if (daughter1Label == FatJetLabel::Z_mumu && daughter2Label == FatJetLabel::Z_bb ) return FatJetLabel::H_ZZ_mumu_bb;
  else if (daughter1Label == FatJetLabel::Z_mumu && daughter2Label == FatJetLabel::Z_ee ) return FatJetLabel::H_ZZ_mumu_ee;
  else if (daughter1Label == FatJetLabel::Z_mumu && daughter2Label == FatJetLabel::Z_mumu ) return FatJetLabel::H_ZZ_mumu_mumu;
  else if (daughter1Label == FatJetLabel::Z_mumu && daughter2Label == FatJetLabel::Z_tautau ) return FatJetLabel::H_ZZ_mumu_tautau;
  else if (daughter1Label == FatJetLabel::Z_mumu && daughter2Label == FatJetLabel::Z_nunu ) return FatJetLabel::H_ZZ_mumu_nunu;	     
  else if (daughter1Label == FatJetLabel::Z_tautau && daughter2Label == FatJetLabel::Z_dd ) return FatJetLabel::H_ZZ_tautau_dd;
  else if (daughter1Label == FatJetLabel::Z_tautau && daughter2Label == FatJetLabel::Z_cc ) return FatJetLabel::H_ZZ_tautau_cc;
  else if (daughter1Label == FatJetLabel::Z_tautau && daughter2Label == FatJetLabel::Z_bb ) return FatJetLabel::H_ZZ_tautau_bb;
  else if (daughter1Label == FatJetLabel::Z_tautau && daughter2Label == FatJetLabel::Z_ee ) return FatJetLabel::H_ZZ_tautau_ee;
  else if (daughter1Label == FatJetLabel::Z_tautau && daughter2Label == FatJetLabel::Z_mumu ) return FatJetLabel::H_ZZ_tautau_mumu;
  else if (daughter1Label == FatJetLabel::Z_tautau && daughter2Label == FatJetLabel::Z_tautau ) return FatJetLabel::H_ZZ_tautau_tautau;
  else if (daughter1Label == FatJetLabel::Z_tautau && daughter2Label == FatJetLabel::Z_nunu ) return FatJetLabel::H_ZZ_tautau_nunu;	     
  else if (daughter1Label == FatJetLabel::Z_nunu && daughter2Label == FatJetLabel::Z_dd ) return FatJetLabel::H_ZZ_nunu_dd;
  else if (daughter1Label == FatJetLabel::Z_nunu && daughter2Label == FatJetLabel::Z_cc ) return FatJetLabel::H_ZZ_nunu_cc;
  else if (daughter1Label == FatJetLabel::Z_nunu && daughter2Label == FatJetLabel::Z_bb ) return FatJetLabel::H_ZZ_nunu_bb;
  else if (daughter1Label == FatJetLabel::Z_nunu && daughter2Label == FatJetLabel::Z_ee ) return FatJetLabel::H_ZZ_nunu_ee;
  else if (daughter1Label == FatJetLabel::Z_nunu && daughter2Label == FatJetLabel::Z_mumu ) return FatJetLabel::H_ZZ_nunu_mumu;
  else if (daughter1Label == FatJetLabel::Z_nunu && daughter2Label == FatJetLabel::Z_tautau ) return FatJetLabel::H_ZZ_nunu_tautau;
  else if (daughter1Label == FatJetLabel::Z_nunu && daughter2Label == FatJetLabel::Z_nunu ) return FatJetLabel::H_ZZ_nunu_nunu;
  return FatJetLabel::Invalid;
}

std::pair<FatJetMatching::FatJetLabel,std::vector<const reco::GenParticle*> > FatJetMatching::top_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{
  std::vector<const reco::GenParticle*> resultVector;
  // top
  auto top = getFinal(parton);
  resultVector.push_back(top);

  // find the W and test if it's hadronic
  const reco::GenParticle *w_from_top = nullptr, *b_from_top = nullptr;
  for (const auto &dau : top->daughterRefVector()){
    if (std::abs(dau->pdgId()) == ParticleID::p_Wplus){
      w_from_top = getFinal(&(*dau));
    }else if (std::abs(dau->pdgId()) <= ParticleID::p_b){
      // ! use <= p_b ! -- can also have charms etc.
      b_from_top = dynamic_cast<const reco::GenParticle*>(&(*dau));
    }
  }
  if (!w_from_top || !b_from_top) throw std::logic_error("[FatJetMatching::top_label] Cannot find b or W from top decay!");

  if (isHadronic(w_from_top)) {
    if (debug_){
      using namespace std;
      cout << "jet: " << jet->polarP4() << endl;
      cout << "top: "; printGenParticleInfo(top, -1);
      cout << "b:   "; printGenParticleInfo(b_from_top, -1);
      cout << "W:   "; printGenParticleInfo(w_from_top, -1);
    }

    auto wdaus = getDaughterQuarks(w_from_top);
    if (wdaus.size() < 2) throw std::logic_error("[FatJetMatching::top_label] W decay has less than 2 quarks!");
//    if (wdaus.size() >= 2)
    {
      double dr_b     = reco::deltaR(jet->p4(), b_from_top->p4());
      double dr_q1    = reco::deltaR(jet->p4(), wdaus.at(0)->p4());
      double dr_q2    = reco::deltaR(jet->p4(), wdaus.at(1)->p4());
      if (dr_q1 > dr_q2){
        // swap q1 and q2 so that dr_q1<=dr_q2
        std::swap(dr_q1, dr_q2);
        std::swap(wdaus.at(0), wdaus.at(1));
      }

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, b)     : " << dr_b << endl;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
      }

      if (dr_b < distR){
        auto pdgid_q1 = std::abs(wdaus.at(0)->pdgId());
        auto pdgid_q2 = std::abs(wdaus.at(1)->pdgId());
        if (debug_){
          using namespace std;
          cout << "pdgid(q1)        : " << pdgid_q1 << endl;
          cout << "pdgid(q2)        : " << pdgid_q2 << endl;
        }

        if (dr_q1<distR && dr_q2<distR){
          if (pdgid_q1 >= ParticleID::p_c || pdgid_q2 >= ParticleID::p_c) {
            return std::make_pair(FatJetLabel::Top_bcq, resultVector);
          }
          else {
            return std::make_pair(FatJetLabel::Top_bqq, resultVector);
          }
        }else if (dr_q1<distR && dr_q2>=distR){
          if (pdgid_q1 >= ParticleID::p_c){
            return std::make_pair(FatJetLabel::Top_bc, resultVector);
          }else{
            return std::make_pair(FatJetLabel::Top_bq, resultVector);
          }
        }
      }else{
        // test for W if dr(b, jet) > distR
        return w_label(jet, w_from_top, distR);
      }
    }
  } else {
    // leptonic W
    if (debug_){
      using namespace std;
      cout << "jet: " << jet->polarP4() << endl;
      cout << "top: "; printGenParticleInfo(top, -1);
      cout << "b:   "; printGenParticleInfo(b_from_top, -1);
      cout << "W:   "; printGenParticleInfo(w_from_top, -1);
    }

    const reco::GenParticle* lep = nullptr;
    for (unsigned i=0; i<w_from_top->numberOfDaughters(); ++i){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(w_from_top->daughter(i));
      auto pdgid = std::abs(dau->pdgId());
      if (pdgid == ParticleID::p_eminus || pdgid == ParticleID::p_muminus || pdgid == ParticleID::p_tauminus){
        // use final version here!
        lep = getFinal(dau); break;
      }
    }

    if (!lep) throw std::logic_error("[FatJetMatching::top_label] Cannot find charged lepton from leptonic W decay!");

    double dr_b     = reco::deltaR(jet->p4(), b_from_top->p4());
    double dr_l     = reco::deltaR(jet->p4(), lep->p4());
    if (debug_){
      using namespace std;
      cout << "deltaR(jet, b)     : " << dr_b << endl;
      cout << "deltaR(jet, l)     : " << dr_l << endl;
      cout << "pdgid(l)           : " << lep->pdgId() << endl;
    }

    if (dr_b < distR && dr_l < distR){
      auto pdgid = std::abs(lep->pdgId());
      if (pdgid == ParticleID::p_eminus){
        return std::make_pair(FatJetLabel::Top_bele, resultVector);
      } else if (pdgid == ParticleID::p_muminus){
        return std::make_pair(FatJetLabel::Top_bmu, resultVector);
      } else if (pdgid == ParticleID::p_tauminus){
        return std::make_pair(FatJetLabel::Top_btau, resultVector);
      }
    }
  }
  std::vector<const reco::GenParticle*> emptyVector;
  return std::make_pair(FatJetLabel::Invalid, emptyVector );

}

std::pair<FatJetMatching::FatJetLabel,std::vector<const reco::GenParticle*> > FatJetMatching::w_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{
  std::vector<const reco::GenParticle*> resultVector;
  auto w = getFinal(parton);
  resultVector.push_back(w);

  if (isHadronic(w)) {
    if (debug_){
      using namespace std;
      cout << "jet: " << jet->polarP4() << endl;
      cout << "W:   "; printGenParticleInfo(w, -1);
    }

    auto wdaus = getDaughterQuarks(w);
    if (wdaus.size() < 2) throw std::logic_error("[FatJetMatching::w_label] W decay has less than 2 quarks!");
//    if (wdaus.size() >= 2)
    {
      double dr_q1    = reco::deltaR(jet->p4(), wdaus.at(0)->p4());
      double dr_q2    = reco::deltaR(jet->p4(), wdaus.at(1)->p4());
      if (dr_q1 > dr_q2){
        // swap q1 and q2 so that dr_q1<=dr_q2
        std::swap(dr_q1, dr_q2);
        std::swap(wdaus.at(0), wdaus.at(1));
      }
      auto pdgid_q1 = std::abs(wdaus.at(0)->pdgId());
      auto pdgid_q2 = std::abs(wdaus.at(1)->pdgId());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
        cout << "pdgid(q1)        : " << pdgid_q1 << endl;
        cout << "pdgid(q2)        : " << pdgid_q2 << endl;
      }

      if (dr_q1<distR && dr_q2<distR){
        if (pdgid_q1 >= ParticleID::p_c || pdgid_q2 >= ParticleID::p_c) {
          return std::make_pair(FatJetLabel::W_cq, resultVector);
        }
        else {
          return std::make_pair(FatJetLabel::W_qq, resultVector);
        }
      }
    }
  }

  std::vector<const reco::GenParticle*> emptyVector;
  return std::make_pair(FatJetLabel::Invalid, emptyVector);

}

std::pair<FatJetMatching::FatJetLabel,std::vector<const reco::GenParticle*> > FatJetMatching::z_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{
  std::vector<const reco::GenParticle*> resultVector;
  auto z = getFinal(parton);
  resultVector.push_back(z);

  if (isHadronic(z)) {
    if (debug_){
      using namespace std;
      cout << "jet: " << jet->polarP4() << endl;
      cout << "Z:   "; printGenParticleInfo(z, -1);
    }

    auto zdaus = getDaughterQuarks(z);
    if (zdaus.size() < 2) throw std::logic_error("[FatJetMatching::z_label] Z decay has less than 2 quarks!");
//    if (zdaus.size() >= 2)
    {
      double dr_q1    = reco::deltaR(jet->p4(), zdaus.at(0)->p4());
      double dr_q2    = reco::deltaR(jet->p4(), zdaus.at(1)->p4());
      if (dr_q1 > dr_q2){
        // swap q1 and q2 so that dr_q1<=dr_q2
        std::swap(dr_q1, dr_q2);
        std::swap(zdaus.at(0), zdaus.at(1));
      }
      auto pdgid_q1 = std::abs(zdaus.at(0)->pdgId());
      auto pdgid_q2 = std::abs(zdaus.at(1)->pdgId());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
        cout << "pdgid(q1)        : " << pdgid_q1 << endl;
        cout << "pdgid(q2)        : " << pdgid_q2 << endl;
      }

      if (dr_q1<distR && dr_q2<distR){
        if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
          return std::make_pair(FatJetLabel::Z_bb, resultVector);
        }else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
          return std::make_pair(FatJetLabel::Z_cc, resultVector);
        }else {
          return std::make_pair(FatJetLabel::Z_qq, resultVector);
        }
      }
    }
  }

  std::vector<const reco::GenParticle*> emptyVector;
  return std::make_pair(FatJetLabel::Invalid, emptyVector );

}

std::pair<FatJetMatching::FatJetLabel,std::vector<const reco::GenParticle*> > FatJetMatching::higgs_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{
  std::vector<const reco::GenParticle*> resultVector;   
  auto higgs = getFinal(parton);
  resultVector.push_back(higgs);

  if (debug_){
    using namespace std;
    cout << "jet: " << jet->polarP4() << endl;
    cout << "H:   "; printGenParticleInfo(higgs, -1);
  }

  bool is_hVV = false;
  if (higgs->numberOfDaughters() >= 3) {
    // e.g., h->Vqq or h->qqqq
    is_hVV = true;
  }else {
    // e.g., h->VV*
    for (const auto &p : higgs->daughterRefVector()){
      auto pdgid = std::abs(p->pdgId());
      if (pdgid == ParticleID::p_Wplus || pdgid == ParticleID::p_Z0){
        is_hVV = true;
        break;
      }
    }
  }

  if (is_hVV){
    // h->WW or h->ZZ
    std::vector<const reco::GenParticle*> hVV_daus;
    //W- is defined as daughter1, W+ is defined as daughter2
    //For Z's, we preserve the order in which the Z shows up
    const reco::GenParticle* daughter1 = 0; 
    const reco::GenParticle* daughter2 = 0;
    std::vector<const reco::GenParticle*> hVV_daughters;
    std::vector<const reco::GenParticle*> hVV_granddaughters;

    for (unsigned idau=0; idau<higgs->numberOfDaughters(); ++idau){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(higgs->daughter(idau));
      auto pdgid = std::abs(higgs->daughter(idau)->pdgId());

      if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b){
        hVV_daus.push_back(dau);
      }else{
	
	if (pdgid == ParticleID::p_Z0) {
	  if (!daughter1) daughter1 = getFinal(dau);
	  else if (!daughter2) daughter2 = getFinal(dau);
	  else {
	    std::cout << "Warning: seems like the Higgs has more than two Z boson daughters...check it.\n";
	  }
	} else {
	  if (higgs->daughter(idau)->pdgId() == -1*ParticleID::p_Wplus) {
	    daughter1 = getFinal(dau);
	  }
	  if (higgs->daughter(idau)->pdgId() == ParticleID::p_Wplus) {
	    daughter2 = getFinal(dau);
	  }
	}
	hVV_daughters.push_back(daughter1);
	hVV_daughters.push_back(daughter2);
	const auto d = getDaughterQuarks(getFinal(dau));
        hVV_daus.insert(hVV_daus.end(), d.begin(), d.end());	
      }
    }

    if (!(daughter1 && daughter2)) {
      std::cout << "Error: Did not find both higgs daughters\n";
      std::vector<const reco::GenParticle*> emptyVector;
      return std::make_pair(FatJetLabel::Invalid, emptyVector );
    }

    
    if (daughter1->numberOfDaughters() != 2 
	|| daughter2->numberOfDaughters() != 2) {
      std::cout << "Warning: Higgs daughter daughters did not come in pairs. Check it.\n";
    }
    
    FatJetLabel eventLabel = getHVVLabel(daughter1, daughter2);

    hVV_granddaughters.push_back((const reco::GenParticle*)daughter1->daughter(0));
    hVV_granddaughters.push_back((const reco::GenParticle*)daughter1->daughter(1));
    hVV_granddaughters.push_back((const reco::GenParticle*)daughter2->daughter(0));
    hVV_granddaughters.push_back((const reco::GenParticle*)daughter2->daughter(1));
    double drHiggs    = reco::deltaR(jet->p4(), higgs->p4());
    double dr_d1_d1    = reco::deltaR(jet->p4(), daughter1->daughter(0)->p4());
    double dr_d1_d2    = reco::deltaR(jet->p4(), daughter1->daughter(1)->p4());
    double dr_d2_d1    = reco::deltaR(jet->p4(), daughter2->daughter(0)->p4());
    double dr_d2_d2    = reco::deltaR(jet->p4(), daughter2->daughter(1)->p4());


    if (debug_){
      using namespace std;     
      std::cout << "Checking Event Label Classification\n";
      printGenParticleInfo(higgs, -1);
      printGenParticleInfo(daughter1, -1);
      printGenParticleInfo(daughter2, -1);
      for (const auto * gp : hVV_granddaughters){      
	printGenParticleInfo(gp, -1);	 
      }
      std::cout << "Event Label: " << eventLabel << "\n";
      std::cout << "\n\n";
    }

    int NDaughtersInsideCone = 0;
    if (dr_d1_d1 < distR) NDaughtersInsideCone++;
    if (dr_d1_d2 < distR) NDaughtersInsideCone++;
    if (dr_d2_d1 < distR) NDaughtersInsideCone++;
    if (dr_d2_d2 < distR) NDaughtersInsideCone++;

    if (eventLabel != FatJetLabel::Invalid) {
      if ( drHiggs < distR && NDaughtersInsideCone >= 4
	   ) {
	resultVector.push_back(daughter1);
	resultVector.push_back(daughter2);
	return std::make_pair(eventLabel, resultVector);
      } else {
	if (debug_){
	  std::cout << "warning: daughters fell outside of fat jet cone.\n";
	  std::cout << "dr_d1_d1 : " << dr_d1_d1 << "\n";
	  std::cout << "dr_d1_d2 : " << dr_d1_d2 << "\n";
	  std::cout << "dr_d2_d1 : " << dr_d2_d1 << "\n";
	  std::cout << "dr_d2_d2 : " << dr_d2_d2 << "\n";
	  printGenParticleInfo(higgs, -1);
	  printGenParticleInfo(daughter1, -1);
	  printGenParticleInfo(daughter2, -1);
	  for (const auto * gp : hVV_granddaughters){      
	    printGenParticleInfo(gp, -1);	 
	  }
	  std::cout << "Event Label: " << eventLabel << "\n";
	}
      }
    } else {
      std::cout << "sixie error: invalid event label\n";
      printGenParticleInfo(higgs, -1);
      printGenParticleInfo(daughter1, -1);
      printGenParticleInfo(daughter2, -1);
      for (const auto * gp : hVV_granddaughters){      
	printGenParticleInfo(gp, -1);	 
      }
      std::cout << "Event Label: " << eventLabel << "\n";
      std::cout << "\n\n";
    }     
  } else if (isHadronic(higgs)) {
    // direct h->qq

    auto hdaus = getDaughterQuarks(higgs);
    if (hdaus.size() < 2) throw std::logic_error("[FatJetMatching::higgs_label] Higgs decay has less than 2 quarks!");
//    if (zdaus.size() >= 2)
    {
      double dr_q1    = reco::deltaR(jet->p4(), hdaus.at(0)->p4());
      double dr_q2    = reco::deltaR(jet->p4(), hdaus.at(1)->p4());
      if (dr_q1 > dr_q2){
        // swap q1 and q2 so that dr_q1<=dr_q2
        std::swap(dr_q1, dr_q2);
        std::swap(hdaus.at(0), hdaus.at(1));
      }
      auto pdgid_q1 = std::abs(hdaus.at(0)->pdgId());
      auto pdgid_q2 = std::abs(hdaus.at(1)->pdgId());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
        cout << "pdgid(q1)        : " << pdgid_q1 << endl;
        cout << "pdgid(q2)        : " << pdgid_q2 << endl;
      }

      if (dr_q1<distR && dr_q2<distR){
        if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
          return std::make_pair(FatJetLabel::H_bb, resultVector);
        }else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
	  return std::make_pair(FatJetLabel::H_cc, resultVector);
        }else {
          return std::make_pair(FatJetLabel::H_qq, resultVector);
        }
      }
    }
  }else {
    // test h->tautau
    std::vector<const reco::GenParticle*> taus;
    for (unsigned i=0; i<higgs->numberOfDaughters(); ++i){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(higgs->daughter(i));
      if (std::abs(dau->pdgId()) == ParticleID::p_tauminus){
        taus.push_back(dau);
      }
    }
    if (taus.size() == 2){
      // higgs -> tautau
      // use first version or last version of the tau in dr?
      double dr_tau1    = reco::deltaR(jet->p4(), taus.at(0)->p4());
      double dr_tau2    = reco::deltaR(jet->p4(), taus.at(1)->p4());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, tau1)    : " << dr_tau1 << endl;
        cout << "deltaR(jet, tau2)    : " << dr_tau2 << endl;
      }

      auto isHadronicTau = [](const reco::GenParticle* tau){
        for (const auto &dau : tau->daughterRefVector()){
          auto pdgid = std::abs(dau->pdgId());
          if (pdgid==ParticleID::p_eminus || pdgid==ParticleID::p_muminus){
            return false;
          }
        }
        return true;
      };

      auto tau1 = getFinal(taus.at(0));
      auto tau2 = getFinal(taus.at(1));
      if (dr_tau1<distR && dr_tau2<distR){
        if (isHadronicTau(tau1) && isHadronicTau(tau2)) {
          return std::make_pair(FatJetLabel::H_tautau, resultVector);
        }
      }
    }
  }

  std::vector<const reco::GenParticle*> emptyVector;
  return std::make_pair(FatJetLabel::Invalid, emptyVector );

}

std::pair<FatJetMatching::FatJetLabel,std::vector<const reco::GenParticle*> > FatJetMatching::qcd_label(const pat::Jet* jet, const reco::GenParticleCollection& genParticles, double distR)
{
  std::vector<const reco::GenParticle*> resultVector;
  const reco::GenParticle *parton = nullptr;
  double minDR = 999;
  for (const auto &gp : genParticles){
    if (gp.status() != 23) continue;
    auto pdgid = std::abs(gp.pdgId());
    if (!(pdgid<ParticleID::p_t || pdgid==ParticleID::p_g)) continue;
    auto dr = reco::deltaR(gp, *jet);
    if (dr<distR && dr<minDR){
      minDR = dr;
      parton = &gp;
    }
  }
  if (debug_){
    using namespace std;
    if (parton){
      cout << "parton"; printGenParticleInfo(parton, -1);
      cout << "dr(jet, parton): " << minDR << endl;
    }
  }

  auto n_bHadrons = jet->jetFlavourInfo().getbHadrons().size();
  auto n_cHadrons = jet->jetFlavourInfo().getcHadrons().size();

  resultVector.push_back(parton);
  if (n_bHadrons>=2) {
    return std::make_pair(FatJetLabel::QCD_bb, resultVector);
  }else if (n_bHadrons==1){
    return std::make_pair(FatJetLabel::QCD_b, resultVector);
  }else if (n_cHadrons>=2){
    return std::make_pair(FatJetLabel::QCD_cc, resultVector);
  }else if (n_cHadrons==1){
    return std::make_pair(FatJetLabel::QCD_c, resultVector);
  }

  return std::make_pair(FatJetLabel::QCD_others, resultVector);
}


