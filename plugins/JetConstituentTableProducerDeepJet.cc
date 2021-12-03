#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "RecoBTag/FeatureTools/interface/sorting_modules.h"

#include <string> 

using namespace btagbtvdeep;

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

template<typename T>
class JetConstituentTableProducerDeepJet : public edm::stream::EDProducer<> {
public:
  explicit JetConstituentTableProducerDeepJet(const edm::ParameterSet &);
  ~JetConstituentTableProducerDeepJet() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  typedef reco::VertexCollection VertexCollection;
  //=====
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;

  const std::string nameDeepJet_;
  const std::string idx_nameDeepJet_;
  const bool readBtag_;
  const double jet_radius_;
  const bool add_DeepJet_noclip_;

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
  
};

//
// constructors and destructor
//
template< typename T>
JetConstituentTableProducerDeepJet<T>::JetConstituentTableProducerDeepJet(const edm::ParameterSet &iConfig)
    : nameDeepJet_(iConfig.getParameter<std::string>("nameDeepJet")),
      idx_nameDeepJet_(iConfig.getParameter<std::string>("idx_nameDeepJet")),
      readBtag_(iConfig.getParameter<bool>("readBtag")),
      jet_radius_(iConfig.getParameter<double>("jet_radius")),
      add_DeepJet_noclip_(iConfig.getParameter<bool>("add_DeepJet_noclip")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))){
  produces<nanoaod::FlatTable>(nameDeepJet_);
}

template< typename T>
JetConstituentTableProducerDeepJet<T>::~JetConstituentTableProducerDeepJet() {}

template< typename T>
void JetConstituentTableProducerDeepJet<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  //auto outCands = std::make_unique<std::vector<reco::CandidatePtr>>();
  //auto outSVs = std::make_unique<std::vector<const reco::VertexCompositePtrCandidate *>> ();
  std::vector<int> jetIdx_dj;
  
    
  auto jets = iEvent.getHandle(jet_token_);
  iEvent.getByToken(vtx_token_, vtxs_);
  iEvent.getByToken(cand_token_, cands_);
  iEvent.getByToken(sv_token_, svs_);

  if(readBtag_){
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);
  }

  
  std::vector<int> jet_N_CPFCands(jets->size());
  std::vector<int> jet_N_NPFCands(jets->size());
  std::vector<int> jet_N_SVs(jets->size()); 
    
    
  // should default to 0 if less than 25 cpf with information
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackEtaRel_0to24(25, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPtRel_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPPar_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackDeltaR_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPParRatio_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip2dVal_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip2dSig_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip3dVal_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip3dSig_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackJetDistVal_0to24(25, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> Cpfcan_ptrel_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_drminsv_0to24(25, std::vector<float>(jets->size()));
  std::vector<std::vector<int>> Cpfcan_VTX_ass_0to24(25, std::vector<int>(jets->size())); 
  std::vector<std::vector<float>> Cpfcan_puppiw_0to24(25, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> Cpfcan_chi2_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<int>> Cpfcan_quality_0to24(25, std::vector<int>(jets->size()));
  // no clip versions
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackEtaRel_0to24_noclip;
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPtRel_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPPar_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackDeltaR_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPParRatio_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip2dVal_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip2dSig_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip3dVal_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip3dSig_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackJetDistVal_0to24_noclip;
  std::vector<std::vector<float>> Cpfcan_ptrel_0to24_noclip; 
  std::vector<std::vector<float>> Cpfcan_drminsv_0to24_noclip;
  std::vector<std::vector<float>> Cpfcan_chi2_0to24_noclip; 
    
    
  // should default to 0 if less than 25 npf with information
  std::vector<std::vector<float>> Npfcan_ptrel_0to24(25, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> Npfcan_deltaR_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<int>> Npfcan_isGamma_0to24(25, std::vector<int>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_HadFrac_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_drminsv_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_puppiw_0to24(25, std::vector<float>(jets->size())); 
  // no clip versions
  std::vector<std::vector<float>> Npfcan_ptrel_0to24_noclip;
  std::vector<std::vector<float>> Npfcan_deltaR_0to24_noclip; 
  std::vector<std::vector<float>> Npfcan_drminsv_0to24_noclip; 
    

  // should default to 0 if less than four SVs with information
  std::vector<std::vector<float>> sv_mass_0to3(4, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> sv_pt_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<int>> sv_ntracks_0to3(4, std::vector<int>(jets->size())); 
  std::vector<std::vector<float>> sv_chi2_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_normchi2_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_dxy_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_dxysig_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_d3d_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_d3dsig_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_costhetasvpv_0to3(4, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> sv_ptrel_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_phirel_0to3(4, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> sv_deltaR_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_enratio_0to3(4, std::vector<float>(jets->size()));
  // no clip versions
  std::vector<std::vector<float>> sv_normchi2_0to3_noclip;
  std::vector<std::vector<float>> sv_deltaR_0to3_noclip; 
  std::vector<std::vector<float>> sv_dxysig_0to3_noclip; 
  std::vector<std::vector<float>> sv_d3dsig_0to3_noclip;
  
  //std::cout << "Successfully initialized all vectors." << std::endl;    
    
  //std::cout << "In producer: add_DeepJet_noclip_ ? " << add_DeepJet_noclip_ << std::endl;  
  
  if (add_DeepJet_noclip_) {
    Cpfcan_BtagPf_trackEtaRel_0to24_noclip.resize(25, std::vector<float>(jets->size()));
    Cpfcan_BtagPf_trackPtRel_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_BtagPf_trackPPar_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_BtagPf_trackDeltaR_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_BtagPf_trackPParRatio_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_BtagPf_trackSip2dVal_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_BtagPf_trackSip2dSig_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_BtagPf_trackSip3dVal_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_BtagPf_trackSip3dSig_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_BtagPf_trackJetDistVal_0to24_noclip.resize(25, std::vector<float>(jets->size()));
    Cpfcan_ptrel_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Cpfcan_drminsv_0to24_noclip.resize(25, std::vector<float>(jets->size()));
    Cpfcan_chi2_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
      
    Npfcan_ptrel_0to24_noclip.resize(25, std::vector<float>(jets->size()));
    Npfcan_deltaR_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
    Npfcan_drminsv_0to24_noclip.resize(25, std::vector<float>(jets->size())); 
      
    sv_normchi2_0to3_noclip.resize(4, std::vector<float>(jets->size()));
    sv_deltaR_0to3_noclip.resize(4, std::vector<float>(jets->size())); 
    sv_dxysig_0to3_noclip.resize(4, std::vector<float>(jets->size())); 
    sv_d3dsig_0to3_noclip.resize(4, std::vector<float>(jets->size()));
      
    //std::cout << "Successfully resized noclip vectors." << std::endl;  
  }
  
  //std::cout << "Start jet loop." << std::endl;  
    
  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet = jets->at(i_jet);
    math::XYZVector jet_dir = jet.momentum().Unit();
    GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
    VertexDistance3D vdist;

    pv_ = &vtxs_->at(0);
    
    //////////////////////
    // Secondary Vertices
    std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
    std::vector<const reco::VertexCompositePtrCandidate *> allSVs;  
      
    
    jetIdx_dj.push_back(i_jet);
      
    jet_N_SVs[i_jet] = 0;
    for (const auto &sv : *svs_) {
      // Factor in cuts in NanoAOD for indexing
      Measurement1D dl= vdist.distance(vtxs_->front(), VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
      if(dl.significance() > 3){
        allSVs.push_back(&sv);
      }
      if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
        jetSVs.push_back(&sv);
        jet_N_SVs[i_jet]++;
      }
    }
      
      
    // sort by dxy significance
    std::sort(jetSVs.begin(),
              jetSVs.end(),
              [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) {
                return sv_vertex_comparator(*sva, *svb, *pv_);
              });
    
    // counter to get flat info per jet for SVs
    unsigned i_sv_in_jet = 0;  

    
    for (const auto &sv : jetSVs) {

      if (readBtag_ && !vtxs_->empty()) {
          
        if (i_sv_in_jet < 4) {   
          sv_mass_0to3[i_sv_in_jet][i_jet] = sv->mass();
          sv_pt_0to3[i_sv_in_jet][i_jet] = sv->pt();
          sv_ntracks_0to3[i_sv_in_jet][i_jet] = sv->numberOfDaughters();
          sv_chi2_0to3[i_sv_in_jet][i_jet] = sv->vertexChi2();
          sv_normchi2_0to3[i_sv_in_jet][i_jet] = catch_infs_and_bound(sv->vertexChi2() / sv->vertexNdof(), 1000, -1000, 1000);
          const auto& dxy_meas = vertexDxy(*sv, *pv_);
          sv_dxy_0to3[i_sv_in_jet][i_jet] = dxy_meas.value();
          sv_dxysig_0to3[i_sv_in_jet][i_jet] = catch_infs_and_bound(dxy_meas.value() / dxy_meas.error(), 0, -1, 800);
          const auto& d3d_meas = vertexD3d(*sv, *pv_);
          sv_d3d_0to3[i_sv_in_jet][i_jet] = d3d_meas.value();
          sv_d3dsig_0to3[i_sv_in_jet][i_jet] = catch_infs_and_bound(d3d_meas.value() / d3d_meas.error(), 0, -1, 800);
          sv_costhetasvpv_0to3[i_sv_in_jet][i_jet] = vertexDdotP(*sv, *pv_);
          // Jet related
          sv_ptrel_0to3[i_sv_in_jet][i_jet] = sv->pt() / jet.pt();
          sv_phirel_0to3[i_sv_in_jet][i_jet] = reco::deltaPhi(*sv, jet);
          sv_deltaR_0to3[i_sv_in_jet][i_jet] = catch_infs_and_bound(std::fabs(reco::deltaR(*sv, jet_dir)) - 0.5, 0, -2, 0);
          sv_enratio_0to3[i_sv_in_jet][i_jet] = sv->energy() / jet.energy();
            
          // no clip
          if (add_DeepJet_noclip_) {
            sv_normchi2_0to3_noclip[i_sv_in_jet][i_jet] = sv->vertexChi2() / sv->vertexNdof();
            sv_dxysig_0to3_noclip[i_sv_in_jet][i_jet] = dxy_meas.value() / dxy_meas.error();
            sv_d3dsig_0to3_noclip[i_sv_in_jet][i_jet] = d3d_meas.value() / d3d_meas.error();
            sv_deltaR_0to3_noclip[i_sv_in_jet][i_jet] = std::fabs(reco::deltaR(*sv, jet_dir));
          }
        } else {
                continue;
        }
          
      } 
      i_sv_in_jet++;
    }

    //std::cout << "Successfully filled SV info, start sorting PF cands now." << std::endl;    
      
    // PF Cands    
    std::vector<reco::CandidatePtr> const & daughters = jet.daughterPtrVector();

    const auto& svs_unsorted = *svs_;  
      
      
      
    std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted, n_sorted;
      
      
      
    // first time looping over all pf candidates
    //     to fill sorted indices and get a connection back to the old indices
    
    jet_N_CPFCands[i_jet] = 0;  
    jet_N_NPFCands[i_jet] = 0;  
      
    for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
      auto cand = jet.daughter(i);
      //if ( cand.isNull() ) continue;
      //auto const *packedCand = dynamic_cast <pat::PackedCandidate const *>(cand.get());
      //if ( packedCand == nullptr ) continue;
      if (cand) {
/*
        // not used currently as it was not part of the JetConstituentTableProducer at first, only for DeepFlavourTagInfoProducer
        // candidates under 950MeV (configurable) are not considered
        // might change if we use also white-listing
        if (cand->pt() < min_candidate_pt_)
          continue;
*/
        if (cand->charge() != 0) {
          //auto& trackinfo = trackinfos.emplace(i, track_builder).first->second;
          //trackinfo.buildTrackInfo(cand, jet_dir, jet_ref_track_dir, *pv_); // *pv_ is an atlernative to vtxs_->at(0)
          btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
            // similar to https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h
          trkinfo.buildTrackInfo(cand, jet_dir, jet_ref_track_dir, vtxs_->at(0));
          c_sorted.emplace_back(
              i, trkinfo.getTrackSip2dSig(), -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet.pt());
          jet_N_CPFCands[i_jet]++;
        } else {
          n_sorted.emplace_back(i, -1, -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet.pt());
          jet_N_NPFCands[i_jet]++;
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
  
      
    //std::cout << "Start looping over PF cands to fill info." << std::endl;    
      
    int i_pf_in_jet = 0;
    for (const auto &cand : daughters) {
      
      auto candPtrs = cands_->ptrs();
      auto candInNewList = std::find( candPtrs.begin(), candPtrs.end(), cand );
      if ( candInNewList == candPtrs.end() ) {
        //std::cout << "Cannot find candidate : " << cand.id() << ", " << cand.key() << ", pt = " << cand->pt() << std::endl;
        continue;
      }
/*
      outCands->push_back(cand);
      jetIdx_pf.push_back(i_jet);
      pfcandIdx.push_back(candInNewList - candPtrs.begin());
      //cand_pt.push_back(cand->pt());
*/
      if (readBtag_ && !vtxs_->empty()) {
        if ( cand.isNull() ) continue;
        auto const *packedCand = dynamic_cast <pat::PackedCandidate const *>(cand.get());
        if ( packedCand == nullptr ) continue;
        float drminpfcandsv = btagbtvdeep::mindrsvpfcand(svs_unsorted, &(*cand));
        if ( packedCand->charge() != 0 ) {
            
            // is charged candidate
            auto entry_sorted = c_sortedindices.at(i_pf_in_jet);
            //std::cout << "Current candidate is " << i_pf_in_jet << " and entry_sorted = c_sortedindices.at(i_pf_in_jet) = " << entry_sorted << std::endl; 
            // need only the first 25 cpfs for DeepJet
            if (entry_sorted > 24) {
                continue;
            }
            
            Cpfcan_puppiw_0to24[entry_sorted][i_jet] = packedCand->puppiWeight();
            Cpfcan_VTX_ass_0to24[entry_sorted][i_jet] = packedCand->pvAssociationQuality();
            Cpfcan_drminsv_0to24[entry_sorted][i_jet] = catch_infs_and_bound(drminpfcandsv, 0, -1. * jet_radius_, 0, -1. * jet_radius_);
            Cpfcan_ptrel_0to24[entry_sorted][i_jet] = catch_infs_and_bound(packedCand->pt() / jet.pt(), 0, -1, 0, -1);
            
            // no clip
            if (add_DeepJet_noclip_) {
              Cpfcan_drminsv_0to24_noclip[entry_sorted][i_jet] = drminpfcandsv;
              Cpfcan_ptrel_0to24_noclip[entry_sorted][i_jet] = packedCand->pt() / jet.pt();
            }
            
            if ( packedCand && packedCand->hasTrackDetails()){
              const auto& pseudo_track = packedCand->pseudoTrack();
              Cpfcan_chi2_0to24[entry_sorted][i_jet] = catch_infs_and_bound(pseudo_track.normalizedChi2(), 300, -1, 300);
                
              // this returns the quality enum not a mask.
              Cpfcan_quality_0to24[entry_sorted][i_jet] = pseudo_track.qualityMask();  
                
                
              btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
                // similar to https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h
              trkinfo.buildTrackInfo(&(*packedCand), jet_dir, jet_ref_track_dir, vtxs_->at(0));
              
              
              // with clip
              Cpfcan_BtagPf_trackEtaRel_0to24[entry_sorted][i_jet]     = catch_infs_and_bound(trkinfo.getTrackEtaRel(), 0, -5, 15);
              Cpfcan_BtagPf_trackPtRel_0to24[entry_sorted][i_jet]      = catch_infs_and_bound(trkinfo.getTrackPtRel(), 0, -1, 4);
              Cpfcan_BtagPf_trackPPar_0to24[entry_sorted][i_jet]       = catch_infs_and_bound(trkinfo.getTrackPPar(), 0, -1e5, 1e5);
              Cpfcan_BtagPf_trackDeltaR_0to24[entry_sorted][i_jet]     = catch_infs_and_bound(trkinfo.getTrackDeltaR(), 0, -5, 5);
              Cpfcan_BtagPf_trackPParRatio_0to24[entry_sorted][i_jet]  = catch_infs_and_bound(trkinfo.getTrackPParRatio(), 0, -10, 100);
              Cpfcan_BtagPf_trackSip2dVal_0to24[entry_sorted][i_jet]   = catch_infs_and_bound(trkinfo.getTrackSip2dVal(), 0, -1, 70);
              Cpfcan_BtagPf_trackSip2dSig_0to24[entry_sorted][i_jet]   = catch_infs_and_bound(trkinfo.getTrackSip2dSig(), 0, -1, 4e4);
              Cpfcan_BtagPf_trackSip3dVal_0to24[entry_sorted][i_jet]   = catch_infs_and_bound(trkinfo.getTrackSip3dVal(), 0, -1, 1e5);
              Cpfcan_BtagPf_trackSip3dSig_0to24[entry_sorted][i_jet]   = catch_infs_and_bound(trkinfo.getTrackSip3dSig(), 0, -1, 4e4);
              Cpfcan_BtagPf_trackJetDistVal_0to24[entry_sorted][i_jet] = catch_infs_and_bound(trkinfo.getTrackJetDistVal(), 0, -20, 1);  
                
              // no clip was default inside JetConstituentTableProducer, but DeepJet seems to use the clipped versions
              if (add_DeepJet_noclip_) {
                Cpfcan_BtagPf_trackEtaRel_0to24_noclip[entry_sorted][i_jet]     = trkinfo.getTrackEtaRel();
                Cpfcan_BtagPf_trackPtRel_0to24_noclip[entry_sorted][i_jet]      = trkinfo.getTrackPtRel();
                Cpfcan_BtagPf_trackPPar_0to24_noclip[entry_sorted][i_jet]       = trkinfo.getTrackPPar();
                Cpfcan_BtagPf_trackDeltaR_0to24_noclip[entry_sorted][i_jet]     = trkinfo.getTrackDeltaR();
                Cpfcan_BtagPf_trackPParRatio_0to24_noclip[entry_sorted][i_jet]  = trkinfo.getTrackPParRatio();
                Cpfcan_BtagPf_trackSip2dVal_0to24_noclip[entry_sorted][i_jet]   = trkinfo.getTrackSip2dVal();
                Cpfcan_BtagPf_trackSip2dSig_0to24_noclip[entry_sorted][i_jet]   = trkinfo.getTrackSip2dSig();
                Cpfcan_BtagPf_trackSip3dVal_0to24_noclip[entry_sorted][i_jet]   = trkinfo.getTrackSip3dVal();
                Cpfcan_BtagPf_trackSip3dSig_0to24_noclip[entry_sorted][i_jet]   = trkinfo.getTrackSip3dSig();
                Cpfcan_BtagPf_trackJetDistVal_0to24_noclip[entry_sorted][i_jet] = trkinfo.getTrackJetDistVal();
                
                Cpfcan_chi2_0to24_noclip[entry_sorted][i_jet] = pseudo_track.normalizedChi2();
              }
                
              //c_pf_features.btagPf_trackPtRatio    = catch_infs_and_bound(track_info.getTrackPtRatio(), 0, -1, 10);
                
              //Cpfcan_ptrel.push_back(catch_infs_and_bound(cand->pt() / jet.pt(), 0, -1, 0, -1));
                
              // note: in JetConstituentTableProducer, there is an additional 0.4 inside the mindrsvpfcand call below
              // but this is not present for other producers I looked at
              //Cpfcan_drminsv.push_back(catch_infs_and_bound(mindrsvpfcand(svs_unsorted, &(*cand), 0.4), 0, -1. * jet_radius_, 0, -1. * jet_radius_));
            } else {
                    // default negative chi2 and loose track if notTrackDetails
                    Cpfcan_chi2_0to24[entry_sorted][i_jet] = catch_infs_and_bound(-1, 300, -1, 300);
                    Cpfcan_quality_0to24[entry_sorted][i_jet] = (1 << reco::TrackBase::loose);
                    
                    
                    // no clip
                    if (add_DeepJet_noclip_) {
                      Cpfcan_chi2_0to24_noclip[entry_sorted][i_jet] = -1;
                    }
                
                    // vector defaults to 0 anyway (initialization), no need to push_back zeros here for other variables

           }
        } else {
            
            
            // is neutral candidate
            auto entry_sorted = n_sortedindices.at(i_pf_in_jet);
            //std::cout << "Current candidate is " << i_pf_in_jet << " and entry_sorted = n_sortedindices.at(i_pf_in_jet) = " << entry_sorted << std::endl; 
            // need only the first 25 cpfs for DeepJet
            if (entry_sorted > 24) {
                continue;
            }
            
            Npfcan_puppiw_0to24[entry_sorted][i_jet] = packedCand->puppiWeight();
            Npfcan_HadFrac_0to24[entry_sorted][i_jet] = packedCand->hcalFraction();
            if (std::abs(packedCand->pdgId()) == 22)
              Npfcan_isGamma_0to24[entry_sorted][i_jet] = 1;
            // catch
            Npfcan_deltaR_0to24[entry_sorted][i_jet] = catch_infs_and_bound(reco::deltaR(*packedCand, jet), 0, -0.6, 0, -0.6);
            Npfcan_drminsv_0to24[entry_sorted][i_jet] = catch_infs_and_bound(drminpfcandsv, 0, -1. * jet_radius_, 0, -1. * jet_radius_);
            Npfcan_ptrel_0to24[entry_sorted][i_jet] = catch_infs_and_bound(packedCand->pt() / jet.pt(), 0, -1, 0, -1);
            // for all cases where catch_infs_and_bound appears, also store the raw quantities without any modifications (no clip)
            if (add_DeepJet_noclip_) {
              Npfcan_deltaR_0to24_noclip[entry_sorted][i_jet] = reco::deltaR(*packedCand, jet);
              Npfcan_drminsv_0to24_noclip[entry_sorted][i_jet] = drminpfcandsv;
              Npfcan_ptrel_0to24_noclip[entry_sorted][i_jet] = packedCand->pt() / jet.pt();
            }
        }
      }
      i_pf_in_jet++;
    }  // end jet loop
  }

    
  //std::cout << "Successfully made it through the jet loop, now starting filling the columns into the table." << std::endl; 
    
  // DeepJetInputs table
  auto djTable = std::make_unique<nanoaod::FlatTable>(jetIdx_dj.size(), nameDeepJet_, false, true);
  djTable->addColumn<int>("DeepJet_jetIdx", jetIdx_dj, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
    
    
  djTable->addColumn<int>("DeepJet_nCpfcand", jet_N_CPFCands, "Number of charged PF candidates in the jet", nanoaod::FlatTable::IntColumn);
  djTable->addColumn<int>("DeepJet_nNpfcand", jet_N_NPFCands, "Number of neutral PF candidates in the jet", nanoaod::FlatTable::IntColumn);
  djTable->addColumn<int>("DeepJet_nsv", jet_N_SVs, "Number of secondary vertices in the jet", nanoaod::FlatTable::IntColumn);
    
    
  std::string input_name;
  std::string description;  
  for (unsigned int p = 0; p < 25; p++) {
      auto s = std::to_string(p);
      
      
      // ============================================================== Cpfs ===================================================================
      input_name = "DeepJet_Cpfcan_puppiw_" + s;
      description = "charged candidate PUPPI weight of the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_puppiw_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_VTX_ass_" + s;
      description = "integer flag that indicates whether the track was used in the primary vertex fit for the " + s + ". cpf";
      djTable->addColumn<int>(input_name, Cpfcan_VTX_ass_0to24[p], description, nanoaod::FlatTable::IntColumn, 10);
      input_name = "DeepJet_Cpfcan_drminsv_" + s;
      description = "track pseudoangular distance from the closest secondary vertex of the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_drminsv_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_ptrel_" + s;
      description = "fraction of the jet momentum carried by the track for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_ptrel_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_quality_" + s;
      description = "integer flag which indicates the quality of the fitted track, based on number of detector hits used for the reconstruction as well as the overall chi2 of the charged track fit for the " + s + ". cpf";
      djTable->addColumn<int>(input_name, Cpfcan_quality_0to24[p], description, nanoaod::FlatTable::IntColumn, 10);
      input_name = "DeepJet_Cpfcan_chi2_" + s;
      description = "chi2 of the charged track fit for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_chi2_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      
      input_name = "DeepJet_Cpfcan_BtagPf_trackDeltaR_" + s;
      description = "track pseudoangular distance from the jet axis for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackDeltaR_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackEtaRel_" + s;
      description = "track pseudorapidity, relative to the jet axis for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackEtaRel_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackJetDistVal_" + s;
      description = "minimum track approach distance to jet axis for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackJetDistVal_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackPPar_" + s;
      description = "dot product of the jet and track momentum for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPPar_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackPParRatio_" + s;
      description = "dot product of the jet and track momentum divided by the magnitude of the jet momentum for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPParRatio_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackPtRel_" + s;
      description = "track transverse momentum, relative to the jet axis for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPtRel_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackSip2dSig_" + s;
      description = "track 2D signed impact parameter significance for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip2dSig_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackSip3dSig_" + s;
      description = "track 3D signed impact parameter significance for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip3dSig_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackSip2dVal_" + s;
      description = "track 2D signed impact parameter for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip2dVal_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Cpfcan_BtagPf_trackSip3dVal_" + s;
      description = "track 3D signed impact parameter for the " + s + ". cpf";
      djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip3dVal_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      
      // ============================================================== Npfs ===================================================================
      input_name = "DeepJet_Npfcan_puppiw_" + s;
      description = "neutral candidate PUPPI weight for the " + s + ". npf";
      djTable->addColumn<float>(input_name, Npfcan_puppiw_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Npfcan_deltaR_" + s;
      description = "pseudoangular distance between the neutral candidate and the jet axis for the " + s + ". npf";
      djTable->addColumn<float>(input_name, Npfcan_deltaR_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Npfcan_drminsv_" + s;
      description = "pseudoangular distance between the neutral candidate and the closest secondary vertex for the " + s + ". npf";
      djTable->addColumn<float>(input_name, Npfcan_drminsv_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Npfcan_HadFrac_" + s;
      description = "fraction of the neutral candidate energy deposited in the hadronic calorimeter for the " + s + ". npf";
      djTable->addColumn<float>(input_name, Npfcan_HadFrac_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Npfcan_ptrel_" + s;
      description = "fraction of the jet momentum carried by the neutral candidate for the " + s + ". npf";
      djTable->addColumn<float>(input_name, Npfcan_ptrel_0to24[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_Npfcan_isGamma_" + s;
      description = "integer flag indicating whether the neutral candidate is a photon for the " + s + ". npf";
      djTable->addColumn<int>(input_name, Npfcan_isGamma_0to24[p], description, nanoaod::FlatTable::IntColumn, 10);
      
      
      if (add_DeepJet_noclip_) {
          // cpf
          input_name = "DeepJet_Cpfcan_drminsv_noclip_" + s;
          description = "track pseudoangular distance from the closest secondary vertex of the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_drminsv_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_ptrel_noclip_" + s;
          description = "fraction of the jet momentum carried by the track for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_ptrel_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_chi2_noclip_" + s;
          description = "chi2 of the charged track fit for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_chi2_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          
          // cpf from track_builder
          input_name = "DeepJet_Cpfcan_BtagPf_trackDeltaR_noclip_" + s;
          description = "track pseudoangular distance from the jet axis for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackDeltaR_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackEtaRel_noclip_" + s;
          description = "track pseudorapidity, relative to the jet axis for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackEtaRel_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackJetDistVal_noclip_" + s;
          description = "minimum track approach distance to jet axis for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackJetDistVal_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackPPar_noclip_" + s;
          description = "dot product of the jet and track momentum for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPPar_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackPParRatio_noclip_" + s;
          description = "dot product of the jet and track momentum divided by the magnitude of the jet momentum for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPParRatio_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackPtRel_noclip_" + s;
          description = "track transverse momentum, relative to the jet axis for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPtRel_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackSip2dSig_noclip_" + s;
          description = "track 2D signed impact parameter significance for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip2dSig_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackSip3dSig_noclip_" + s;
          description = "track 3D signed impact parameter significance for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip3dSig_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackSip2dVal_noclip_" + s;
          description = "track 2D signed impact parameter for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip2dVal_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Cpfcan_BtagPf_trackSip3dVal_noclip_" + s;
          description = "track 3D signed impact parameter for the " + s + ". cpf (no clip)";
          djTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip3dVal_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          
          // npf
          input_name = "DeepJet_Npfcan_deltaR_noclip_" + s;
          description = "pseudoangular distance between the neutral candidate and the jet axis for the " + s + ". npf (no clip)";
          djTable->addColumn<float>(input_name, Npfcan_deltaR_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Npfcan_drminsv_noclip_" + s;
          description = "pseudoangular distance between the neutral candidate and the closest secondary vertex for the " + s + ". npf (no clip)";
          djTable->addColumn<float>(input_name, Npfcan_drminsv_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_Npfcan_ptrel_noclip_" + s;
          description = "fraction of the jet momentum carried by the neutral candidate for the " + s + ". npf (no clip)";
          djTable->addColumn<float>(input_name, Npfcan_ptrel_0to24_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
      }
      
  }  
    
  //std::cout << "Successfully filled the PFCands columns into the table, now start with SVs." << std::endl; 
    
  // ============================================================== SVs ===================================================================
  for (unsigned int p = 0; p < 4; p++) {
      auto s = std::to_string(p);
      
      input_name = "DeepJet_sv_mass_" + s;
      description = "SV mass of the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_mass_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_pt_" + s;
      description = "SV pt of the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_pt_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_ntracks_" + s;
      description = "Number of tracks asociated to the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_ntracks_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_chi2_" + s;
      description = "chi2 of the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_chi2_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_normchi2_" + s;
      description = "chi2/dof of the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_normchi2_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_dxy_" + s;
      description = "2D impact parameter (flight distance) value of the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_dxy_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_dxysig_" + s;
      description = "2D impact parameter (flight distance) significance of the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_dxysig_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_d3d_" + s;
      description = "3D impact parameter (flight distance) value of the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_d3d_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_d3dsig_" + s;
      description = "3D impact parameter (flight distance) significance of the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_d3dsig_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_costhetasvpv_" + s;
      description = "cosine of the angle between the " + s + ". SV flight direction and the direction of the " + s + ". SV momentum";
      djTable->addColumn<float>(input_name, sv_costhetasvpv_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJetExtra_sv_phirel_" + s;
      description = "DeltaPhi(sv, jet) for the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_phirel_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJetExtra_sv_ptrel_" + s;
      description = "pT relative to parent jet for the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_ptrel_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_deltaR_" + s;
      description = "pseudoangular distance between jet axis and the " + s + ". SV direction";
      djTable->addColumn<float>(input_name, sv_deltaR_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_enratio_" + s;
      description = "ratio of the " + s + ". SV energy ratio to the jet energy";
      djTable->addColumn<float>(input_name, sv_enratio_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      
      if (add_DeepJet_noclip_) {
          input_name = "DeepJet_sv_normchi2_noclip_" + s;
          description = "chi2/dof of the " + s + ". SV (no clip)";
          djTable->addColumn<float>(input_name, sv_normchi2_0to3_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_sv_dxysig_noclip_" + s;
          description = "2D impact parameter (flight distance) significance of the " + s + ". SV (no clip)";
          djTable->addColumn<float>(input_name, sv_dxysig_0to3_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_sv_d3dsig_noclip_" + s;
          description = "3D impact parameter (flight distance) significance of the " + s + ". SV (no clip)";
          djTable->addColumn<float>(input_name, sv_d3dsig_0to3_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "DeepJet_sv_deltaR_noclip_" + s;
          description = "pseudoangular distance between jet axis and the " + s + ". SV direction (no clip)";
          djTable->addColumn<float>(input_name, sv_deltaR_0to3_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
      }
  }
  
    
  iEvent.put(std::move(djTable), nameDeepJet_);  
  
}

template< typename T>
void JetConstituentTableProducerDeepJet<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("nameDeepJet", "Jet");
  desc.add<std::string>("idx_nameDeepJet", "djIdx");
  desc.add<double>("jet_radius", true);
  desc.add<bool>("readBtag", true);
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK8"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("slimmedSecondaryVertices"));
  desc.add<bool>("add_DeepJet_noclip", false);
  descriptions.addWithDefaultLabel(desc);
}

typedef JetConstituentTableProducerDeepJet<pat::Jet> PatJetConstituentTableProducerDeepJet;
//typedef JetConstituentTableProducerDeepJet<reco::GenJet> GenJetConstituentTableProducerDeepJet;

DEFINE_FWK_MODULE(PatJetConstituentTableProducerDeepJet);
//DEFINE_FWK_MODULE(GenJetConstituentTableProducerDeepJet);
