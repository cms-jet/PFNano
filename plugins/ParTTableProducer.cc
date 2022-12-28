#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "RecoBTag/FeatureTools/interface/sorting_modules.h"


#include "TVector3.h"

#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

using btagbtvdeep::SortingClass;
using btagbtvdeep::mindrsvpfcand;
using btagbtvdeep::invertSortingVector;

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <string> 

// add tag info and a way to go back to the jet reference
#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/DeepFlavourTagInfo.h"


class TrackInfoBuilder{
public:
    TrackInfoBuilder(edm::ESHandle<TransientTrackBuilder> & build):
        builder(build),
        trackMomentum_(0),
        trackEta_(0),
        trackEtaRel_(0),
        trackPtRel_(0),
        trackPPar_(0),
        trackDeltaR_(0),
        trackPtRatio_(0),
        trackPParRatio_(0),
        trackSip2dVal_(0),
        trackSip2dSig_(0),
        trackSip3dVal_(0),
        trackSip3dSig_(0),

        trackJetDistVal_(0),
        trackJetDistSig_(0),
	//trackMinDistSV_(0),
	//trackMinDistXSV_(0),
	//trackMinDistYSV_(0),
	//trackMinDistZSV_(0)
	ttrack_(0)

{


}

    void buildTrackInfo(const pat::PackedCandidate* PackedCandidate_ ,const math::XYZVector&  jetDir, GlobalVector refjetdirection, const reco::Vertex & pv){
			TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());
			if(!PackedCandidate_->hasTrackDetails()) {
				TVector3 trackMom3(
					PackedCandidate_->momentum().x(),
					PackedCandidate_->momentum().y(),
					PackedCandidate_->momentum().z()
					);
				trackMomentum_=PackedCandidate_->p();
				trackEta_= PackedCandidate_->eta();
				trackEtaRel_=reco::btau::etaRel(jetDir, PackedCandidate_->momentum());
				trackPtRel_=trackMom3.Perp(jetDir3);
				trackPPar_=jetDir.Dot(PackedCandidate_->momentum());
				trackDeltaR_=reco::deltaR(PackedCandidate_->momentum(), jetDir);
				trackPtRatio_=trackMom3.Perp(jetDir3) / PackedCandidate_->p();
				trackPParRatio_=jetDir.Dot(PackedCandidate_->momentum()) / PackedCandidate_->p();
				trackSip2dVal_=0.;
				trackSip2dSig_=0.;
				trackSip3dVal_=0.;
				trackSip3dSig_=0.;
				trackJetDistVal_=0.;
				trackJetDistSig_=0.;
				//ttrack_ = 0;
				//trackMinDistSV_=0.;
				//trackMinDistXSV_=0.;
				//trackMinDistYSV_=0.;
				//trackMinDistZSV_=0.;
				return;
			}

        const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();

        reco::TransientTrack transientTrack;
        //reco::TransientTrack* transientTrack2;
        transientTrack=builder->build(PseudoTrack);
        Measurement1D meas_ip2d=IPTools::signedTransverseImpactParameter(transientTrack, refjetdirection, pv).second;
        Measurement1D meas_ip3d=IPTools::signedImpactParameter3D(transientTrack, refjetdirection, pv).second;
        Measurement1D jetdist=IPTools::jetTrackDistance(transientTrack, refjetdirection, pv).second;
        math::XYZVector trackMom = PseudoTrack.momentum();
        double trackMag = std::sqrt(trackMom.Mag2());
        TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());
	//float mindistsv = mindistsvpfcand(transientTrack);
	//GlobalPoint mindistgpsv = mingpsvpfcand(transientTrack);


        trackMomentum_=std::sqrt(trackMom.Mag2());
        trackEta_= trackMom.Eta();
        trackEtaRel_=reco::btau::etaRel(jetDir, trackMom);
        trackPtRel_=trackMom3.Perp(jetDir3);
        trackPPar_=jetDir.Dot(trackMom);
        trackDeltaR_=reco::deltaR(trackMom, jetDir);
        trackPtRatio_=trackMom3.Perp(jetDir3) / trackMag;
        trackPParRatio_=jetDir.Dot(trackMom) / trackMag;
        trackSip2dVal_=(meas_ip2d.value());

        trackSip2dSig_=(meas_ip2d.significance());
        trackSip3dVal_=(meas_ip3d.value());


        trackSip3dSig_= meas_ip3d.significance();
        trackJetDistVal_= jetdist.value();
        trackJetDistSig_= jetdist.significance();

	//trackMinDistSV_= mindistsv;
	//trackMinDistXSV_= mindistgpsv.x();
	//trackMinDistYSV_= mindistgpsv.y();
	//trackMinDistZSV_= mindistgpsv.z();
	ttrack_ = transientTrack;


    }

    const float& getTrackDeltaR() const {return trackDeltaR_;}
    const float& getTrackEta() const {return trackEta_;}
    const float& getTrackEtaRel() const {return trackEtaRel_;}
    const float& getTrackJetDistSig() const {return trackJetDistSig_;}
    const float& getTrackJetDistVal() const {return trackJetDistVal_;}
    const float& getTrackMomentum() const {return trackMomentum_;}
    const float& getTrackPPar() const {return trackPPar_;}
    const float& getTrackPParRatio() const {return trackPParRatio_;}
    const float& getTrackPtRatio() const {return trackPtRatio_;}
    const float& getTrackPtRel() const {return trackPtRel_;}
    const float& getTrackSip2dSig() const {return trackSip2dSig_;}
    const float& getTrackSip2dVal() const {return trackSip2dVal_;}
    const float& getTrackSip3dSig() const {return trackSip3dSig_;}
    const float& getTrackSip3dVal() const {return trackSip3dVal_;}

    const reco::TransientTrack getTTrack() const {return ttrack_;}

  //const float& getrackMinDistSV() const {return trackMinDistSV_;}
  //const float& getTrackMinDistXSV() const {return trackMinDistXSV_;}
  //const float& getTrackMinDistYSV() const {return trackMinDistYSV_;}
  //const float& getTrackMinDistZSV() const {return trackMinDistZSV_;}

private:

    edm::ESHandle<TransientTrackBuilder>& builder;

    float trackMomentum_;
    float trackEta_;
    float trackEtaRel_;
    float trackPtRel_;
    float trackPPar_;
    float trackDeltaR_;
    float trackPtRatio_;
    float trackPParRatio_;
    float trackSip2dVal_;
    float trackSip2dSig_;
    float trackSip3dVal_;
    float trackSip3dSig_;

    float trackJetDistVal_;
    float trackJetDistSig_;
    reco::TransientTrack ttrack_;
  //float trackMinDistSV_;
  //float trackMinDistXSV_;
  //float trackMinDistYSV_;
  //float trackMinDistZSV_;

};


template<typename T>
class ParTTableProducer : public edm::stream::EDProducer<> {
public:
  explicit ParTTableProducer(const edm::ParameterSet &);
  ~ParTTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  static float mindistsvpfcand(const reco::TransientTrack track, const reco::VertexCompositePtrCandidateCollection svs_unsorted);
    
  
private:
  void produce(edm::Event &, const edm::EventSetup &) override;
    
  typedef reco::VertexCollection VertexCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
    
    
    
  const std::string nameParT_;
  const std::string idx_nameParT_;

  edm::EDGetTokenT<edm::View<T>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<reco::CandidateView> cand_token_;
  edm::EDGetTokenT<SVCollection> sv_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<reco::CandidateView> cands_;
  edm::Handle<SVCollection> svs_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;

  const reco::Vertex *pv_ = nullptr;
    
  const float min_candidate_pt_ = 0.95;
    
    
  constexpr static unsigned n_cpf_ = 25;
  constexpr static unsigned n_npf_ = 25;
  constexpr static unsigned n_sv_ = 5;
    
  constexpr static double jetR_ = 0.4;    
    
  
};

//
// constructors and destructor
//
template< typename T>
ParTTableProducer<T>::ParTTableProducer(const edm::ParameterSet &iConfig)
    : nameParT_(iConfig.getParameter<std::string>("nameParT")),
      idx_nameParT_(iConfig.getParameter<std::string>("idx_nameParT")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))){
  produces<nanoaod::FlatTable>(nameParT_);
}

template< typename T>
ParTTableProducer<T>::~ParTTableProducer() {}

template< typename T>
void ParTTableProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
    
  // only necessary to explicitly check correct matching of jets
  std::vector<int> jetIdx_dj;
    
  auto jets = iEvent.getHandle(jet_token_);
  iEvent.getByToken(vtx_token_, vtxs_);
  iEvent.getByToken(cand_token_, cands_);
  iEvent.getByToken(sv_token_, svs_);
    
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);
    
    
  unsigned nJets = jets->size();
    
  // new variables for ParT
    
  std::vector<std::vector<float>> Cpfcan_distminsv_nCpf(n_cpf_, std::vector<float>(nJets));
  // 4 vectors
  std::vector<std::vector<float>> Cpfcan_pt_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_eta_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_phi_nCpf(n_cpf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Cpfcan_e_nCpf(n_cpf_, std::vector<float>(nJets)); 
    

  std::vector<std::vector<float>> Npfcan_etarel_nNpf(n_npf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Npfcan_phirel_nNpf(n_npf_, std::vector<float>(nJets));
  // 4 vectors
  std::vector<std::vector<float>> Npfcan_pt_nNpf(n_npf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Npfcan_eta_nNpf(n_npf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Npfcan_phi_nNpf(n_npf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Npfcan_e_nNpf(n_npf_, std::vector<float>(nJets)); 

  
  std::vector<std::vector<float>> sv_etarel_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_phirel_nSV(n_sv_, std::vector<float>(nJets));
  // 4 vectors (pt was already included in DeepJet)
  std::vector<std::vector<float>> sv_eta_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_phi_nSV(n_sv_, std::vector<float>(nJets));
  std::vector<std::vector<float>> sv_e_nSV(n_sv_, std::vector<float>(nJets)); 
    
    
  // old variables, but use a fifth SV as well for ParT
  std::vector<float> sv_mass_4(nJets);
  std::vector<float> sv_pt_4(nJets); 
  std::vector<float> sv_ntracks_4(nJets); 
  std::vector<float> sv_chi2_4(nJets); 
  std::vector<float> sv_normchi2_4(nJets); 
  std::vector<float> sv_dxy_4(nJets); 
  std::vector<float> sv_dxysig_4(nJets); 
  std::vector<float> sv_d3d_4(nJets); 
  std::vector<float> sv_d3dsig_4(nJets); 
  std::vector<float> sv_costhetasvpv_4(nJets);
  std::vector<float> sv_deltaR_4(nJets); 
  std::vector<float> sv_enratio_4(nJets); 
    
    
  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet = jets->at(i_jet);
      
    const float jet_uncorr_pt=jet.correctedJet("Uncorrected").pt();
    const float jet_uncorr_e=jet.correctedJet("Uncorrected").energy();
      
    math::XYZVector jet_dir = jet.momentum().Unit();
    GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
    VertexDistance3D vdist;

    pv_ = &vtxs_->at(0);
    
    //////////////////////
    // Secondary Vertices
    std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
    std::vector<const reco::VertexCompositePtrCandidate *> allSVs;  
      
    
    jetIdx_dj.push_back(i_jet);
      
    for (const auto &sv : *svs_) {
      // Factor in cuts in NanoAOD for indexing
      Measurement1D dl= vdist.distance(vtxs_->front(), VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
      if(dl.significance() > 3){
        allSVs.push_back(&sv);
      }
      if (reco::deltaR2(sv, jet) < jetR_ * jetR_) {
        jetSVs.push_back(&sv);
      }
    }
      
      
    // sort by dxy significance
    std::sort(jetSVs.begin(),
              jetSVs.end(),
              [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) {
                return btagbtvdeep::sv_vertex_comparator(*sva, *svb, *pv_);
              });
    
    // counter to get flat info per jet for SVs
    unsigned i_sv_in_jet = 0;
    
    for (const auto &sv : jetSVs) {
      //if (readBtag_ && !vtxs_->empty()) {
      if (!vtxs_->empty()) {
          
        if (i_sv_in_jet < 4) {
            
          // this is new for ParT:
          sv_etarel_nSV[i_sv_in_jet][i_jet] = btagbtvdeep::catch_infs_and_bound(std::fabs(sv->eta()-jet.eta())-0.5,0,-2,0);
          sv_phirel_nSV[i_sv_in_jet][i_jet] = btagbtvdeep::catch_infs_and_bound(std::fabs(reco::deltaPhi(sv->phi(),jet.phi()))-0.5,0,-2,0);
            
          sv_eta_nSV[i_sv_in_jet][i_jet] = sv->eta();
          sv_phi_nSV[i_sv_in_jet][i_jet] = sv->phi();
          sv_e_nSV[i_sv_in_jet][i_jet] = sv->energy();
            
        } else {
            if (i_sv_in_jet == 4) { // the fifth SV, new for ParT, not yet for DeepJet
                sv_pt_4[i_jet] = sv->pt();
                sv_deltaR_4[i_jet] = btagbtvdeep::catch_infs_and_bound(std::fabs(reco::deltaR(*sv, jet)) - 0.5, 0, -2, 0);
                sv_mass_4[i_jet] = sv->mass();
                sv_ntracks_4[i_jet] = sv->numberOfDaughters();
                sv_chi2_4[i_jet] = sv->vertexChi2();
                sv_normchi2_4[i_jet] = btagbtvdeep::catch_infs_and_bound(sv_chi2_4[i_jet] / sv->vertexNdof(), 1000, -1000, 1000);
                const auto& dxy_meas = btagbtvdeep::vertexDxy(*sv, *pv_);
                sv_dxy_4[i_jet] = dxy_meas.value();
                sv_dxysig_4[i_jet] = btagbtvdeep::catch_infs_and_bound(sv_dxy_4[i_jet] / dxy_meas.error(), 0, -1, 800);
                const auto& d3d_meas = btagbtvdeep::vertexD3d(*sv, *pv_);
                sv_d3d_4[i_jet] = d3d_meas.value();
                sv_d3dsig_4[i_jet] = btagbtvdeep::catch_infs_and_bound(sv_d3d_4[i_jet] / d3d_meas.error(), 0, -1, 800);
                sv_costhetasvpv_4[i_jet] = btagbtvdeep::vertexDdotP(*sv, *pv_);
                sv_enratio_4[i_jet] = sv->energy() / jet_uncorr_e;
                
                
                // also the fifth SV needs new features:
                sv_etarel_nSV[i_sv_in_jet][i_jet] = btagbtvdeep::catch_infs_and_bound(std::fabs(sv->eta()-jet.eta())-0.5,0,-2,0);
                sv_phirel_nSV[i_sv_in_jet][i_jet] = btagbtvdeep::catch_infs_and_bound(std::fabs(reco::deltaPhi(sv->phi(),jet.phi()))-0.5,0,-2,0);
            
                sv_eta_nSV[i_sv_in_jet][i_jet] = sv->eta();
                sv_phi_nSV[i_sv_in_jet][i_jet] = sv->phi();
                sv_e_nSV[i_sv_in_jet][i_jet] = sv->energy();
            } else { // more than 5 SV irrelevant for this tagger
                continue;
            }
        }
          
      } 
      i_sv_in_jet++;
    } // end sv loop

    // PF Cands    
    std::vector<reco::CandidatePtr> const & daughters = jet.daughterPtrVector();

    const auto& svs_unsorted = *svs_;  
      
    std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted, n_sorted;
      
    // first time looping over all pf candidates
    //     to fill sorted indices and get a connection back to the old indices
    for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
      auto cand = jet.daughter(i);
      if (cand) {
        // candidates under 950MeV (configurable) are not considered
        // might change if we use also white-listing
        if (cand->pt() < min_candidate_pt_) continue;

        if (cand->charge() != 0) {
          btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
            // similar to https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h
          trkinfo.buildTrackInfo(cand, jet_dir, jet_ref_track_dir, *pv_);
          c_sorted.emplace_back(
              i, trkinfo.getTrackSip2dSig(), -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet_uncorr_pt);
        } else {
          n_sorted.emplace_back(i, -1, -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet_uncorr_pt);
        }
      }
    }
      
    // sort collections (open the black-box if you please)
    std::sort(c_sorted.begin(), c_sorted.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);
    std::sort(n_sorted.begin(), n_sorted.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);

    std::vector<size_t> c_sortedindices, n_sortedindices;

    // this puts 0 everywhere and the right position in ind
    c_sortedindices = btagbtvdeep::invertSortingVector(c_sorted);
    n_sortedindices = btagbtvdeep::invertSortingVector(n_sorted);
    

    int i_pf_in_jet = 0;
    for (const auto &cand : daughters) {
      
      auto candPtrs = cands_->ptrs();
      auto candInNewList = std::find( candPtrs.begin(), candPtrs.end(), cand );
      if ( candInNewList == candPtrs.end() ) {
        continue;
      }

      if (!vtxs_->empty()) {
        if ( cand.isNull() ) continue;
        auto const *packedCand = dynamic_cast <pat::PackedCandidate const *>(cand.get());
        if ( packedCand == nullptr ) continue;
        if (packedCand->pt() < min_candidate_pt_) continue;
          
        if ( packedCand->charge() != 0 ) {
            
            // is charged candidate
            auto entry_sorted = c_sortedindices.at(i_pf_in_jet);
            // need only the first 25 cpfs for ParT
            if (entry_sorted >= n_cpf_) {
                continue;
            }
           
            Cpfcan_pt_nCpf[entry_sorted][i_jet] = cand->pt();
            Cpfcan_eta_nCpf[entry_sorted][i_jet] = cand->eta();
            Cpfcan_phi_nCpf[entry_sorted][i_jet] = cand->phi();
            Cpfcan_e_nCpf[entry_sorted][i_jet] = cand->energy();
            
            
            if ( packedCand && packedCand->hasTrackDetails()){
                
              TrackInfoBuilder trkinfo(track_builder_);
                // similar to https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h
              trkinfo.buildTrackInfo(&(*packedCand), jet_dir, jet_ref_track_dir, *pv_);
              
              const reco::TransientTrack ttrack = trkinfo.getTTrack();
              float mindistsv = mindistsvpfcand(ttrack, svs_unsorted);

              Cpfcan_distminsv_nCpf[entry_sorted][i_jet] = mindistsv;  
                
                
            } else {
                continue;
                // nothing else new to fill for ParT which does not have track details for cpfs
           }
        } else {
            
            // is neutral candidate
            auto entry_sorted = n_sortedindices.at(i_pf_in_jet);
            // need only the first 25 npfs for ParT
            if (entry_sorted >= n_npf_) {
                continue;
            }
            Npfcan_etarel_nNpf[entry_sorted][i_jet] = btagbtvdeep::catch_infs_and_bound(std::fabs(cand->eta()-jet.eta()), 0, -2, 0, -0.5);
            Npfcan_phirel_nNpf[entry_sorted][i_jet] = btagbtvdeep::catch_infs_and_bound(std::fabs(reco::deltaPhi(cand->phi(),jet.phi())), 0, -2, 0, -0.5);
            Npfcan_pt_nNpf[entry_sorted][i_jet] = cand->pt();
            Npfcan_eta_nNpf[entry_sorted][i_jet] = cand->eta();
            Npfcan_phi_nNpf[entry_sorted][i_jet] = cand->phi();
            Npfcan_e_nNpf[entry_sorted][i_jet] = cand->energy();
            
        }
      }
      i_pf_in_jet++;
    }  // end constituents loop
  } // end jet loop

  
  // ParTInputs table
  auto partTable = std::make_unique<nanoaod::FlatTable>(nJets, nameParT_, false, true);
    
    
  // ============================================================== Cpfs ===================================================================
  for (unsigned int p = 0; p < n_cpf_; p++) {
      auto s = std::to_string(p);
      partTable->addColumn<float>("ParT_Cpfcan_pt_" + s,
                                Cpfcan_pt_nCpf[p],
                                "transverse momentum of the " + s + ". cpf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Cpfcan_eta_" + s,
                                Cpfcan_eta_nCpf[p],
                                "pseudorapidity of the " + s + ". cpf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Cpfcan_phi_" + s,
                                Cpfcan_phi_nCpf[p],
                                "phi of the " + s + ". cpf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Cpfcan_e_" + s,
                                Cpfcan_e_nCpf[p],
                                "energy of the " + s + ". cpf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Cpfcan_distminsv_" + s,
                                Cpfcan_distminsv_nCpf[p],
                                "track 3D cartesian distance (signed) from the closest secondary vertex of the " + s + ". cpf",
                                nanoaod::FlatTable::FloatColumn, 10);
  }
      
  // ============================================================== Npfs ===================================================================
  for (unsigned int p = 0; p < n_npf_; p++) {
      auto s = std::to_string(p);
      partTable->addColumn<float>("ParT_Npfcan_etarel_" + s,
                                Npfcan_etarel_nNpf[p],
                                "pseudorapidity relative to the jet for the " + s + ". npf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Npfcan_phirel_" + s,
                                Npfcan_phirel_nNpf[p],
                                "phi relative to the jet for the " + s + ". npf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Npfcan_pt_" + s,
                                Npfcan_pt_nNpf[p],
                                "transverse momentum of the " + s + ". npf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Npfcan_eta_" + s,
                                Npfcan_eta_nNpf[p],
                                "pseudorapidity of the " + s + ". npf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Npfcan_phi_" + s,
                                Npfcan_phi_nNpf[p],
                                "phi of the " + s + ". npf",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_Npfcan_e_" + s,
                                Npfcan_e_nNpf[p],
                                "energy of the " + s + ". npf",
                                nanoaod::FlatTable::FloatColumn, 10);
  }  

  // ============================================================== SVs ===================================================================
  for (unsigned int p = 0; p < n_sv_; p++) {
      auto s = std::to_string(p);
      // new for ParT
      partTable->addColumn<float>("ParT_sv_etarel_" + s,
                                sv_etarel_nSV[p], 
                                "pseudorapidity relative to parent jet for the " + s + ". SV",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_sv_phirel_" + s,
                                sv_phirel_nSV[p],
                                "DeltaPhi(sv, jet) for the " + s + ". SV",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_sv_eta_" + s,
                                sv_eta_nSV[p], 
                                "pseudorapidity of the " + s + ". SV",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_sv_phi_" + s,
                                sv_phi_nSV[p],
                                "phi of the " + s + ". SV",
                                nanoaod::FlatTable::FloatColumn, 10);
      partTable->addColumn<float>("ParT_sv_e_" + s,
                                sv_e_nSV[p], 
                                "energy of the " + s + ". SV",
                                nanoaod::FlatTable::FloatColumn, 10);
  }
  
  
  partTable->addColumn<float>("ParT_sv_mass_4",
                            sv_mass_4,
                            "SV mass of the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_pt_4",
                            sv_pt_4,
                            "SV pt of the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_ntracks_4",
                            sv_ntracks_4,
                            "Number of tracks asociated to the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_chi2_4", 
                            sv_chi2_4,
                            "chi2 of the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_normchi2_4",
                            sv_normchi2_4,
                            "chi2/dof of the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_dxy_4",
                            sv_dxy_4,
                            "2D impact parameter (flight distance) value of the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_dxysig_4",
                            sv_dxysig_4,
                            "2D impact parameter (flight distance) significance of the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_d3d_4",
                            sv_d3d_4,
                            "3D impact parameter (flight distance) value of the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_d3dsig_4",
                            sv_d3dsig_4,
                            "3D impact parameter (flight distance) significance of the 4. SV",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_costhetasvpv_4",
                            sv_costhetasvpv_4,
                            "cosine of the angle between the 4. SV flight direction and the direction of the 4. SV momentum",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_deltaR_4",
                            sv_deltaR_4,
                            "pseudoangular distance between jet axis and the 4. SV direction",
                            nanoaod::FlatTable::FloatColumn, 10);
  partTable->addColumn<float>("ParT_sv_enratio_4",
                            sv_enratio_4,
                            "ratio of the 4. SV energy ratio to the jet energy",
                            nanoaod::FlatTable::FloatColumn, 10);
  
  iEvent.put(std::move(partTable), nameParT_);  
  
}

// below OK
template< typename T>
void ParTTableProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("nameParT", "Jet");
  desc.add<std::string>("idx_nameParT", "djIdx");
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("slimmedSecondaryVertices"));
  descriptions.addWithDefaultLabel(desc);
}


template< typename T>
float ParTTableProducer<T>::mindistsvpfcand(const reco::TransientTrack track, const reco::VertexCompositePtrCandidateCollection svs_unsorted) {

    float mindist_ = 999.999;
    float out_dist = 0.0;
    for (unsigned int i=0; i<svs_unsorted.size(); ++i) {
      if (!track.isValid()) {continue;}
      reco::Vertex::CovarianceMatrix csv; svs_unsorted.at(i).fillVertexCovariance(csv);
      reco::Vertex vertex(svs_unsorted.at(i).vertex(), csv);
      if (!vertex.isValid()) {continue;}

      GlobalVector direction(svs_unsorted.at(i).px(),svs_unsorted.at(i).py(),svs_unsorted.at(i).pz());


      AnalyticalImpactPointExtrapolator extrapolator(track.field());
      TrajectoryStateOnSurface tsos =  extrapolator.extrapolate(track.impactPointState(), RecoVertex::convertPos(vertex.position()));

      VertexDistance3D dist;

      if (!tsos.isValid()) {continue;}
      GlobalPoint refPoint = tsos.globalPosition();
      GlobalError refPointErr = tsos.cartesianError().position();
      GlobalPoint vertexPosition = RecoVertex::convertPos(vertex.position());
      GlobalError vertexPositionErr = RecoVertex::convertError(vertex.error());

      std::pair<bool, Measurement1D> result(true, dist.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
      if (!result.first) {continue;}

      GlobalPoint impactPoint = tsos.globalPosition();
      GlobalVector IPVec(impactPoint.x() - vertex.x(), impactPoint.y() - vertex.y(), impactPoint.z() - vertex.z());
      double prod = IPVec.dot(direction);
      double sign = (prod >= 0) ? 1. : -1.;

      if(result.second.value() < mindist_){ 
        out_dist = sign * result.second.value();
        mindist_ = result.second.value();
      } 
    }
    return out_dist;
}

typedef ParTTableProducer<pat::Jet> PatJetParTTableProducer;
//typedef ParTTableProducer<reco::GenJet> GenJetParTTableProducer;

DEFINE_FWK_MODULE(PatJetParTTableProducer);
//DEFINE_FWK_MODULE(GenJetParTTableProducer);
