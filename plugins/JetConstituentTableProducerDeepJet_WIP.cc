#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"

#include "RecoBTag/FeatureTools/interface/JetConverter.h"
#include "RecoBTag/FeatureTools/interface/ShallowTagInfoConverter.h"
#include "RecoBTag/FeatureTools/interface/SecondaryVertexConverter.h"
#include "RecoBTag/FeatureTools/interface/NeutralCandidateConverter.h"
#include "RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
using namespace btagbtvdeep;

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

template<typename T>
class JetConstituentTableProducerDeepJet_WIP : public edm::stream::EDProducer<> {
public:
  explicit JetConstituentTableProducerDeepJet_WIP(const edm::ParameterSet &);
  ~JetConstituentTableProducerDeepJet_WIP() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  typedef reco::VertexCollection VertexCollection;
  //=====
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
    
  typedef edm::View<reco::ShallowTagInfo> ShallowTagInfoCollection;

  //const std::string name_;
  //const std::string name_;
  //const std::string nameSV_;
  const std::string nameDeepJet_;
  //const std::string idx_name_;
  //const std::string idx_nameSV_;
  const std::string idx_nameDeepJet_;
  const bool readBtag_;
  const double jet_radius_;
    
  // from DeepFlavourTagInfoProducer
  const double min_candidate_pt_;
  const bool flip_;

  edm::EDGetTokenT<edm::View<T>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<reco::CandidateView> cand_token_;
  edm::EDGetTokenT<SVCollection> sv_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<reco::CandidateView> cands_;
  edm::Handle<SVCollection> svs_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;

  const reco::Vertex *pv_ = nullptr;
    
    
  // from DeepFlavourTagInfoProducer  
  edm::EDGetTokenT<ShallowTagInfoCollection> shallow_tag_info_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> pvasq_value_map_token_;
  edm::EDGetTokenT<edm::Association<VertexCollection>> pvas_token_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> candidateToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;  
    
  bool use_puppi_value_map_;
  bool use_pvasq_value_map_;

  bool fallback_puppi_weight_;
  bool fallback_vertex_association_;
  
};

//
// constructors and destructor
//
template< typename T>
JetConstituentTableProducerDeepJet_WIP<T>::JetConstituentTableProducerDeepJet_WIP(const edm::ParameterSet &iConfig)
    : //name_(iConfig.getParameter<std::string>("name")),
      //nameSV_(iConfig.getParameter<std::string>("nameSV")),
      nameDeepJet_(iConfig.getParameter<std::string>("nameDeepJet")),
      //idx_name_(iConfig.getParameter<std::string>("idx_name")),
      //idx_nameSV_(iConfig.getParameter<std::string>("idx_nameSV")),
      idx_nameDeepJet_(iConfig.getParameter<std::string>("idx_nameDeepJet")),
      readBtag_(iConfig.getParameter<bool>("readBtag")),
      jet_radius_(iConfig.getParameter<double>("jet_radius")),
      min_candidate_pt_(iConfig.getParameter<double>("min_candidate_pt")),
      flip_(iConfig.getParameter<bool>("flip")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
      shallow_tag_info_token_(
          consumes<ShallowTagInfoCollection>(iConfig.getParameter<edm::InputTag>("shallow_tag_infos"))),
      candidateToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("candidates"))),
      track_builder_token_(
          esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      use_puppi_value_map_(false),
      use_pvasq_value_map_(false),
      fallback_puppi_weight_(iConfig.getParameter<bool>("fallback_puppi_weight")),
      fallback_vertex_association_(iConfig.getParameter<bool>("fallback_vertex_association")){
  //produces<nanoaod::FlatTable>(name_);
  //produces<nanoaod::FlatTable>(name_);
  //produces<nanoaod::FlatTable>(nameSV_);
  produces<nanoaod::FlatTable>(nameDeepJet_);
  //produces<std::vector<reco::CandidatePtr>>();
          
          
  // from DeepFlavourTagInfoProducer
  const auto& puppi_value_map_tag = iConfig.getParameter<edm::InputTag>("puppi_value_map");
  if (!puppi_value_map_tag.label().empty()) {
    puppi_value_map_token_ = consumes<edm::ValueMap<float>>(puppi_value_map_tag);
    use_puppi_value_map_ = true;
  }
          
  const auto& pvas_tag = iConfig.getParameter<edm::InputTag>("vertex_associator");
  if (!pvas_tag.label().empty()) {
    pvasq_value_map_token_ = consumes<edm::ValueMap<int>>(pvas_tag);
    pvas_token_ = consumes<edm::Association<VertexCollection>>(pvas_tag);
    use_pvasq_value_map_ = true;
  }
          
          
}

template< typename T>
JetConstituentTableProducerDeepJet_WIP<T>::~JetConstituentTableProducerDeepJet_WIP() {}

template< typename T>
void JetConstituentTableProducerDeepJet_WIP<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  //auto outCands = std::make_unique<std::vector<reco::CandidatePtr>>();
  //auto outSVs = std::make_unique<std::vector<const reco::VertexCompositePtrCandidate *>> ();
  std::vector<int> jetIdx_dj;
  //std::vector<int> jetIdx_pf, jetIdx_sv, pfcandIdx, svIdx;
/* 
  // PF Cands
  std::vector<float> btagEtaRel, btagPtRatio, btagPParRatio, btagSip3dVal, btagSip3dSig, btagJetDistVal, cand_pt;
  std::vector<float> Cpfcan_BtagPf_trackEtaRel, Cpfcan_BtagPf_trackPtRel, Cpfcan_BtagPf_trackPPar, Cpfcan_BtagPf_trackDeltaR,
    Cpfcan_BtagPf_trackPParRatio, Cpfcan_BtagPf_trackSip2dVal, Cpfcan_BtagPf_trackSip2dSig, Cpfcan_BtagPf_trackSip3dVal,
    Cpfcan_BtagPf_trackSip3dSig, Cpfcan_BtagPf_trackJetDistVal, Cpfcan_ptrel, Cpfcan_drminsv,
    Cpfcan_VTX_ass, Cpfcan_puppiw, Cpfcan_chi2, Cpfcan_quality;
  std::vector<float> Npfcan_ptrel, Npfcan_deltaR, Npfcan_isGamma, Npfcan_HadFrac, Npfcan_drminsv, Npfcan_puppiw;
*/
/*
  // Secondary vertices
  std::vector<float> sv_mass, sv_pt, sv_ntracks, sv_chi2, sv_normchi2, sv_dxy, sv_dxysig, sv_d3d, sv_d3dsig, sv_costhetasvpv;
  std::vector<float> sv_ptrel, sv_phirel, sv_deltaR, sv_enratio;
*/
/*
  // SV but flat
  std::vector<float> sv_mass_0, sv_pt_0, sv_ntracks_0, sv_chi2_0, sv_normchi2_0, sv_dxy_0, sv_dxysig_0, sv_d3d_0, sv_d3dsig_0, sv_costhetasvpv_0;
  std::vector<float> sv_ptrel_0, sv_phirel_0, sv_deltaR_0, sv_enratio_0;
  std::vector<float> sv_mass_1, sv_pt_1, sv_ntracks_1, sv_chi2_1, sv_normchi2_1, sv_dxy_1, sv_dxysig_1, sv_d3d_1, sv_d3dsig_1, sv_costhetasvpv_1;
  std::vector<float> sv_ptrel_1, sv_phirel_1, sv_deltaR_1, sv_enratio_1;
  std::vector<float> sv_mass_2, sv_pt_2, sv_ntracks_2, sv_chi2_2, sv_normchi2_2, sv_dxy_2, sv_dxysig_2, sv_d3d_2, sv_d3dsig_2, sv_costhetasvpv_2;
  std::vector<float> sv_ptrel_2, sv_phirel_2, sv_deltaR_2, sv_enratio_2;
  std::vector<float> sv_mass_3, sv_pt_3, sv_ntracks_3, sv_chi2_3, sv_normchi2_3, sv_dxy_3, sv_dxysig_3, sv_d3d_3, sv_d3dsig_3, sv_costhetasvpv_3;
  std::vector<float> sv_ptrel_3, sv_phirel_3, sv_deltaR_3, sv_enratio_3;
*/
    
  

    
  auto jets = iEvent.getHandle(jet_token_);
  iEvent.getByToken(vtx_token_, vtxs_);
  iEvent.getByToken(cand_token_, cands_);
  iEvent.getByToken(sv_token_, svs_);

  if(readBtag_){
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);
  }
    
  // from DeepFlavourTagInfoProducer
/*    
  edm::Handle<ShallowTagInfoCollection> shallow_tag_infos;
  iEvent.getByToken(shallow_tag_info_token_, shallow_tag_infos);
  double negative_cut = 0;  //used only with flip_
  if (flip_) {              //FIXME: Check if can do even less often than once per event
    const edm::Provenance* prov = shallow_tag_infos.provenance();
    const edm::ParameterSet& psetFromProvenance = edm::parameterSet(prov->stable(), iEvent.processHistory());
    negative_cut = ((psetFromProvenance.getParameter<edm::ParameterSet>("computer"))
                        .getParameter<edm::ParameterSet>("trackSelection"))
                       .getParameter<double>("sip3dSigMax");
  }
*/
  edm::Handle<edm::ValueMap<float>> puppi_value_map;
  if (use_puppi_value_map_) {
    iEvent.getByToken(puppi_value_map_token_, puppi_value_map);
  }

  edm::Handle<edm::ValueMap<int>> pvasq_value_map;
  edm::Handle<edm::Association<VertexCollection>> pvas;
  if (use_pvasq_value_map_) {
    iEvent.getByToken(pvasq_value_map_token_, pvasq_value_map);
    iEvent.getByToken(pvas_token_, pvas);
  }

  edm::ESHandle<TransientTrackBuilder> track_builder = iSetup.getHandle(track_builder_token_);
      
    
    
  
  
    
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
  
  // should default to 0 if less than 25 npf with information
  std::vector<std::vector<float>> Npfcan_ptrel_0to24(25, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> Npfcan_deltaR_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<int>> Npfcan_isGamma_0to24(25, std::vector<int>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_HadFrac_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_drminsv_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_puppiw_0to24(25, std::vector<float>(jets->size())); 

  //auto nJets = jets->size();
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
  
    
  std::cout << "Just after creating the vectors:" << std::endl;
  std::cout << "Size of Cpfcan_puppiw_0to24 is " << Cpfcan_puppiw_0to24.size() << std::endl;
  std::cout << "Size of the first column of Cpfcan_puppiw_0to24 is " << Cpfcan_puppiw_0to24[0].size() << std::endl;  
    
  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet = jets->at(i_jet);

      
      
    // from DeepFlavourTagInfoProducer  
    const auto* pf_jet = dynamic_cast<const reco::PFJet*>(&jet);
    const auto* pat_jet = dynamic_cast<const pat::Jet*>(&jet);
    //edm::RefToBase<reco::Jet> jet_ref(jets, i_jet);  
      
      
    math::XYZVector jet_dir = jet.momentum().Unit();
    GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
    VertexDistance3D vdist;

    pv_ = &vtxs_->at(0);
      
      
      
    //jetIdx_dj.push_back(i_jet);
      
      
      
    //////////////////////
    // Secondary Vertices
    std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
    //std::vector<const reco::VertexCompositePtrCandidate *> allSVs;  
      
      
    for (const auto &sv : *svs_) {
      // Factor in cuts in NanoAOD for indexing
      Measurement1D dl= vdist.distance(vtxs_->front(), VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
/*
      if(dl.significance() > 3){
        allSVs.push_back(&sv);
      }
*/
      if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
        jetSVs.push_back(&sv);
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
      // I'd like to try something like this: https://stackoverflow.com/questions/50870374/better-way-to-map-string-fields-to-variables
      // have a map from key to variables inside the SV loop such that I do not need to write every single constituent (more relevant for PFcands later)
    
    for (const auto &sv : jetSVs) {
/*
      // auto svPtrs = svs_->ptrs();
      auto svInNewList = std::find(allSVs.begin(), allSVs.end(), sv );
      if (svInNewList == allSVs.end()) {
        // continue;
        svIdx.push_back(-1);
      } else{
        svIdx.push_back(svInNewList - allSVs.begin());
      }
      outSVs->push_back(sv);
      jetIdx_sv.push_back(i_jet);
*/
      if (readBtag_ && !vtxs_->empty()) {
/*
        // Jet independent
        sv_mass.push_back(sv->mass());
        sv_pt.push_back(sv->pt());

        sv_ntracks.push_back(sv->numberOfDaughters());
        sv_chi2.push_back(sv->vertexChi2());
        sv_normchi2.push_back(catch_infs_and_bound(sv->vertexChi2() / sv->vertexNdof(), 1000, -1000, 1000));
        const auto& dxy_meas = vertexDxy(*sv, *pv_);
        sv_dxy.push_back(dxy_meas.value());
        sv_dxysig.push_back(catch_infs_and_bound(dxy_meas.value() / dxy_meas.error(), 0, -1, 800));
        const auto& d3d_meas = vertexD3d(*sv, *pv_);
        sv_d3d.push_back(d3d_meas.value());
        sv_d3dsig.push_back(catch_infs_and_bound(d3d_meas.value() / d3d_meas.error(), 0, -1, 800));
        sv_costhetasvpv.push_back(vertexDdotP(*sv, *pv_));
        // Jet related
        sv_ptrel.push_back(sv->pt() / jet.pt());
        sv_phirel.push_back(reco::deltaPhi(*sv, jet));
        sv_deltaR.push_back(catch_infs_and_bound(std::fabs(reco::deltaR(*sv, jet_dir)) - 0.5, 0, -2, 0));
        sv_enratio.push_back(sv->energy() / jet.energy());
*/          
          
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
        } else {
                continue;
        }
          
      } 
      i_sv_in_jet++;
    }

      
      
    // PF Cands    
/*
    std::vector<reco::CandidatePtr> const & daughters = jet.daughterPtrVector();

    const auto& svs_unsorted = *svs_;  
      
     
    std::vector<const eco::CandidatePtr *> cPFs;
    std::vector<const eco::CandidatePtr *> nPFs;  
*/    
      
      
    // ========================================================================================================================
    // this is all copied from DeepFlavourTagInfoProducer
      
    std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted, n_sorted;
    // to cache the TrackInfo
    std::map<unsigned int, btagbtvdeep::TrackInfoBuilder> trackinfos;

    // unsorted reference to sv
    const auto& svs_unsorted = *svs_;
    // fill collection, from DeepTNtuples plus some styling
    for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
      auto cand = jet.daughter(i);
      if (cand) {
        // candidates under 950MeV (configurable) are not considered
        // might change if we use also white-listing
        if (cand->pt() < min_candidate_pt_)
          continue;
        if (cand->charge() != 0) {
          auto& trackinfo = trackinfos.emplace(i, track_builder).first->second;
          trackinfo.buildTrackInfo(cand, jet_dir, jet_ref_track_dir, *pv_); // *pv_ is an atlernative to vtxs_->at(0)
          c_sorted.emplace_back(
              i, trackinfo.getTrackSip2dSig(), -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet.pt());
        } else {
          n_sorted.emplace_back(i, -1, -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet.pt());
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
      
    // ------------------------------------------------------------------------------------------------------------------------
      
    // now: cpfs and npfs with sorted indices are available, one can go from the normal indices to the ones in the sorted collections
      
      
    // need to find out if it's a charged candidate or not, store both types separately to later on perform sorting on the two collections  
/*    for (const auto &dau : daughters) {
      // Factor in cuts in NanoAOD for indexing
      Measurement1D dl= vdist.distance(vtxs_->front(), VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
      if(dl.significance() > 3){
        allSVs.push_back(&sv);
      }
      if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
        jetSVs.push_back(&sv);
      }
    }
*/      
      
    // sort by something t.b.c for PF candidates, and ond only for neutral or charged ones speparately
    // some example for pt sorting:
    // std::sort(daughters.begin(), daughters.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); });
    // stolen from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#Example_code_accessing_all_high
    
/*
    // example for sorting as done for SVs
    std::sort(jetSVs.begin(),
              jetSVs.end(),
              [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) {
                return sv_vertex_comparator(*sva, *svb, *pv_);
              });
*/      
      
      
      
      
      
    if (readBtag_ && !vtxs_->empty()) {  
      for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
          // get pointer and check that is correct
          auto cand = dynamic_cast<const reco::Candidate*>(jet.daughter(i));
          if (!cand)
            continue;
          // candidates under 950MeV are not considered
          // might change if we use also white-listing
          if (cand->pt() < 0.95)
            continue;

          auto packed_cand = dynamic_cast<const pat::PackedCandidate*>(cand);
          auto reco_cand = dynamic_cast<const reco::PFCandidate*>(cand);

          // need some edm::Ptr or edm::Ref if reco candidates
          reco::PFCandidatePtr reco_ptr;
          if (pf_jet) {
            reco_ptr = pf_jet->getPFConstituent(i);
          } else if (pat_jet && reco_cand) {
            reco_ptr = pat_jet->getPFConstituent(i);
          }

          // get PUPPI weight from value map
          float puppiw = 1.0;  // fallback value
          if (reco_cand && use_puppi_value_map_) {
            puppiw = (*puppi_value_map)[reco_ptr];
          } else if (reco_cand && !fallback_puppi_weight_) {
            throw edm::Exception(edm::errors::InvalidReference, "PUPPI value map missing")
                << "use fallback_puppi_weight option to use " << puppiw << "as default";
          }


          float drminpfcandsv = btagbtvdeep::mindrsvpfcand(svs_unsorted, cand);

          if (cand->charge() != 0) {
            // is charged candidate
            auto entry = c_sortedindices.at(i);
            std::cout << "Current candidate is " << i << " and entry = c_sortedindices.at(i) = " << entry << std::endl; 
            // need only the first 25 cpfs for DeepJet
            if (entry > 24) {
                continue;
            }
            // get cached track info
            auto& trackinfo = trackinfos.at(i);
/*              
            if (flip_ && (trackinfo.getTrackSip3dSig() > negative_cut)) {
              continue;
            }
*/            
/*
            // get_ref to vector element
            auto& c_pf_features = features.c_pf_features.at(entry);
*/
            std::cout << "Prior to filling with entries for this candidate:" << std::endl;
            std::cout << "Size of Cpfcan_puppiw_0to24 is " << Cpfcan_puppiw_0to24.size() << std::endl;
            std::cout << "Size of the first column of Cpfcan_puppiw_0to24 is " << Cpfcan_puppiw_0to24[0].size() << std::endl;
            // fill feature structure
            if (packed_cand) {
                Cpfcan_puppiw_0to24[entry][i_jet] = packed_cand->puppiWeight();
            
            std::cout << "After filling with entries for this candidate:" << std::endl;
            std::cout << "Size of Cpfcan_puppiw_0to24 is " << Cpfcan_puppiw_0to24.size() << std::endl;
            std::cout << "Size of the first column of Cpfcan_puppiw_0to24 is " << Cpfcan_puppiw_0to24[0].size() << std::endl;
/*
              btagbtvdeep::packedCandidateToFeatures(
                  packed_cand, jet, trackinfo, drminpfcandsv, static_cast<float>(jet_radius_), c_pf_features, flip_);
*/
            } else if (reco_cand) {
              // get vertex association quality
              int pv_ass_quality = 0;  // fallback value
              if (use_pvasq_value_map_) {
                pv_ass_quality = (*pvasq_value_map)[reco_ptr];
              } else if (!fallback_vertex_association_) {
                throw edm::Exception(edm::errors::InvalidReference, "vertex association missing")
                    << "use fallback_vertex_association option to use" << pv_ass_quality
                    << "as default quality and closest dz PV as criteria";
              }
              // getting the PV as PackedCandidatesProducer
              // but using not the slimmed but original vertices
              auto ctrack = reco_cand->bestTrack();
              int pvi = -1;
              float dist = 1e99;
              for (size_t ii = 0; ii < vtxs_->size(); ii++) {
                float dz = (ctrack) ? std::abs(ctrack->dz(((*vtxs_)[ii]).position())) : 0;
                if (dz < dist) {
                  pvi = ii;
                  dist = dz;
                }
              }
              auto PV = reco::VertexRef(vtxs_, pvi);
              if (use_pvasq_value_map_) {
                const reco::VertexRef& PV_orig = (*pvas)[reco_ptr];
                if (PV_orig.isNonnull())
                  PV = reco::VertexRef(vtxs_, PV_orig.key());
              }
                
              Cpfcan_puppiw_0to24[entry][i_jet] = puppiw;  
/*                
              btagbtvdeep::recoCandidateToFeatures(reco_cand,
                                                   jet,
                                                   trackinfo,
                                                   drminpfcandsv,
                                                   static_cast<float>(jet_radius_),
                                                   puppiw,
                                                   pv_ass_quality,
                                                   PV,
                                                   c_pf_features,
                                                   flip_);
*/
            }
          } else {
            // is neutral candidate
            auto entry = n_sortedindices.at(i);
            // need only the first 25 npfs for DeepJet
            if (entry > 24) {
                continue;
            }
/*
            // get_ref to vector element
            auto& n_pf_features = features.n_pf_features.at(entry);
*/
            // fill feature structure
            if (packed_cand) {
/*
                btagbtvdeep::packedCandidateToFeatures(
                  packed_cand, jet, drminpfcandsv, static_cast<float>(jet_radius_), n_pf_features);
*/
            } else if (reco_cand) {
/*
                btagbtvdeep::recoCandidateToFeatures(
                  reco_cand, jet, drminpfcandsv, static_cast<float>(jet_radius_), puppiw, n_pf_features);
*/
            }
          }
        }

        //output_tag_infos->emplace_back(features, jet_ref);
    }
      
      
      
      
      
      
/*      
    for (const auto &cand : daughters) {
      auto candPtrs = cands_->ptrs();
      auto candInNewList = std::find( candPtrs.begin(), candPtrs.end(), cand );
      if ( candInNewList == candPtrs.end() ) {
        //std::cout << "Cannot find candidate : " << cand.id() << ", " << cand.key() << ", pt = " << cand->pt() << std::endl;
        continue;
      }

      outCands->push_back(cand);
      jetIdx_pf.push_back(i_jet);
      pfcandIdx.push_back(candInNewList - candPtrs.begin());
      //cand_pt.push_back(cand->pt());

      if (readBtag_ && !vtxs_->empty()) {
        if ( cand.isNull() ) continue;
        auto const *packedCand = dynamic_cast <pat::PackedCandidate const *>(cand.get());
        if ( packedCand == nullptr ) continue;
        if ( packedCand && packedCand->hasTrackDetails()){
          btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
            // similar to https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h
          trkinfo.buildTrackInfo(&(*packedCand), jet_dir, jet_ref_track_dir, vtxs_->at(0));
          Cpfcan_BtagPf_trackEtaRel.push_back(trkinfo.getTrackEtaRel());
          Cpfcan_BtagPf_trackPtRel.push_back(trkinfo.getTrackPtRel());
          Cpfcan_BtagPf_trackPPar.push_back(trkinfo.getTrackPPar());
          Cpfcan_BtagPf_trackDeltaR.push_back(trkinfo.getTrackDeltaR());
          Cpfcan_BtagPf_trackPParRatio.push_back(trkinfo.getTrackPParRatio());
          Cpfcan_BtagPf_trackSip2dVal.push_back(trkinfo.getTrackSip2dVal());
          Cpfcan_BtagPf_trackSip2dSig.push_back(trkinfo.getTrackSip2dSig());
          Cpfcan_BtagPf_trackSip3dVal.push_back(trkinfo.getTrackSip3dVal());
          Cpfcan_BtagPf_trackSip3dSig.push_back(trkinfo.getTrackSip3dSig());
          Cpfcan_BtagPf_trackJetDistVal.push_back(trkinfo.getTrackJetDistVal());
          Cpfcan_ptrel.push_back(catch_infs_and_bound(cand->pt() / jet.pt(), 0, -1, 0, -1));
          Cpfcan_drminsv.push_back(catch_infs_and_bound(mindrsvpfcand(svs_unsorted, &(*cand), 0.4), 0, -1. * jet_radius_, 0, -1. * jet_radius_));
        } else {
                Cpfcan_BtagPf_trackEtaRel.push_back(0);
                Cpfcan_BtagPf_trackPtRel.push_back(0);
                Cpfcan_BtagPf_trackPPar.push_back(0);
                Cpfcan_BtagPf_trackDeltaR.push_back(0);
                Cpfcan_BtagPf_trackPParRatio.push_back(0);
                Cpfcan_BtagPf_trackSip2dVal.push_back(0);
                Cpfcan_BtagPf_trackSip2dSig.push_back(0);
                Cpfcan_BtagPf_trackSip3dVal.push_back(0);
                Cpfcan_BtagPf_trackSip3dSig.push_back(0);
                Cpfcan_BtagPf_trackJetDistVal.push_back(0);
                Cpfcan_ptrel.push_back(0);
                Cpfcan_drminsv.push_back(0);

        }
      }
    }  */// end jet loop
  }
/*
  auto candTable = std::make_unique<nanoaod::FlatTable>(outCands->size(), name_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTableProducer
  candTable->addColumn<int>(idx_name_, pfcandIdx, "Index in the candidate list", nanoaod::FlatTable::IntColumn);
  candTable->addColumn<int>("jetIdx", jetIdx_pf, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
  if (readBtag_) {
    //candTable->addColumn<float>("pt", cand_pt, "pt", nanoaod::FlatTable::FloatColumn, 10);  // to check matchind down the line
    candTable->addColumn<float>("Cpfcan_BtagPf_trackEtaRel", Cpfcan_BtagPf_trackEtaRel, "btagEtaRel", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackPtRel", Cpfcan_BtagPf_trackPtRel, "btagPtRel", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackPPar", Cpfcan_BtagPf_trackPPar, "btagPPar", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackDeltaR", Cpfcan_BtagPf_trackDeltaR, "btagDeltaR", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackPParRatio", Cpfcan_BtagPf_trackPParRatio, "btagPParRatio", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackSip2dVal", Cpfcan_BtagPf_trackSip2dVal, "btagSip2dVal", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackSip2dSig", Cpfcan_BtagPf_trackSip2dSig, "btagSip2dSig", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackSip3dVal", Cpfcan_BtagPf_trackSip3dVal, "btagSip3dVal", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackSip3dSig", Cpfcan_BtagPf_trackSip3dSig, "btagSip3dSig", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_BtagPf_trackJetDistVal", Cpfcan_BtagPf_trackJetDistVal, "btagJetDistVal", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_ptrel", Cpfcan_ptrel, "Cpfcan_ptrel", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("Cpfcan_drminsv", Cpfcan_drminsv, "Cpfcan_drminsv", nanoaod::FlatTable::FloatColumn, 10);
  }
  iEvent.put(std::move(candTable), name_);
*/
/*   // SV table
  auto svTable = std::make_unique<nanoaod::FlatTable>(outSVs->size(), nameSV_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTnameableProducer
  svTable->addColumn<int>("jetIdx", jetIdx_sv, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
  svTable->addColumn<int>(idx_nameSV_, svIdx, "Index in the SV list", nanoaod::FlatTable::IntColumn);
  if (readBtag_) {
    svTable->addColumn<float>("mass", sv_mass, "SV mass", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("pt", sv_pt, "SV pt", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("ntracks", sv_ntracks, "Number of tracks associated to SV", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("chi2", sv_chi2, "chi2", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("normchi2", sv_normchi2, "chi2/ndof", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("dxy", sv_dxy, "", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("dxysig", sv_dxysig, "", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("d3d", sv_d3d, "", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("d3dsig", sv_d3dsig, "", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("costhetasvpv", sv_costhetasvpv, "", nanoaod::FlatTable::FloatColumn, 10);
    // Jet related
    svTable->addColumn<float>("phirel", sv_phirel, "DeltaPhi(sv, jet)", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("ptrel", sv_ptrel, "pT relative to parent jet", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("deltaR", sv_deltaR, "dR from parent jet", nanoaod::FlatTable::FloatColumn, 10);
    svTable->addColumn<float>("enratio", sv_enratio, "energy relative to parent jet", nanoaod::FlatTable::FloatColumn, 10);
  }
  iEvent.put(std::move(svTable), nameSV_);
*/  
    
  // DeepJetInputs table
  //auto djTable = std::make_unique<nanoaod::FlatTable>(jets->size(), nameDeepJet_, false);
  auto djTable = std::make_unique<nanoaod::FlatTable>(jetIdx_dj.size(), nameDeepJet_, false, true);
  //djTable->addColumn<int>("DeepJet_jetIdx", jetIdx_dj, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
  
  std::cout << "Just before adding the column:" << std::endl;
  std::cout << "Size of Cpfcan_puppiw_0to24 is " << Cpfcan_puppiw_0to24.size() << std::endl;
  std::cout << "Size of the first column of Cpfcan_puppiw_0to24 is " << Cpfcan_puppiw_0to24[0].size() << std::endl;  
  std::cout << "Size of sv_mass_0to3 is " << sv_mass_0to3.size() << std::endl;
  std::cout << "Size of the first column of sv_mass_0to3 is " << sv_mass_0to3[0].size() << std::endl; 
  
  for (unsigned int k = 0; k < Cpfcan_puppiw_0to24[0].size(); k++) {
      std::cout << "Cpfcan_puppiw_0to24[0][" << k << "] = " << Cpfcan_puppiw_0to24[0][k] << std::endl;
  }
    
  //djTable->addColumn<float>("DeepJet_Cpfcan_puppiw_0", Cpfcan_puppiw_0to24[0], "Cpfcan_puppiw of the 1. cpf", nanoaod::FlatTable::FloatColumn, 10),  
    
    
  djTable->addColumn<float>("DeepJet_sv_mass_0", sv_mass_0to3[0], "SV mass of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_mass_1", sv_mass_0to3[1], "SV mass of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_mass_2", sv_mass_0to3[2], "SV mass of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_mass_3", sv_mass_0to3[3], "SV mass of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_pt_0", sv_pt_0to3[0], "SV pt of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_pt_1", sv_pt_0to3[1], "SV pt of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_pt_2", sv_pt_0to3[2], "SV pt of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_pt_3", sv_pt_0to3[3], "SV pt of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_ntracks_0", sv_ntracks_0to3[0], "Number of tracks asociated to the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_ntracks_1", sv_ntracks_0to3[1], "Number of tracks asociated to the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_ntracks_2", sv_ntracks_0to3[2], "Number of tracks asociated to the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_ntracks_3", sv_ntracks_0to3[3], "Number of tracks asociated to the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_chi2_0", sv_chi2_0to3[0], "chi2 of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_chi2_1", sv_chi2_0to3[1], "chi2 of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_chi2_2", sv_chi2_0to3[2], "chi2 of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_chi2_3", sv_chi2_0to3[3], "chi2 of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_normchi2_0", sv_normchi2_0to3[0], "chi2/dof of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_normchi2_1", sv_normchi2_0to3[1], "chi2/dof of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_normchi2_2", sv_normchi2_0to3[2], "chi2/dof of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_normchi2_3", sv_normchi2_0to3[3], "chi2/dof of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_dxy_0", sv_dxy_0to3[0], "dxy of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_dxy_1", sv_dxy_0to3[1], "dxy of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_dxy_2", sv_dxy_0to3[2], "dxy of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_dxy_3", sv_dxy_0to3[3], "dxy of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_dxysig_0", sv_dxysig_0to3[0], "dxysig of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_dxysig_1", sv_dxysig_0to3[1], "dxysig of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_dxysig_2", sv_dxysig_0to3[2], "dxysig of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_dxysig_3", sv_dxysig_0to3[3], "dxysig of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_d3d_0", sv_d3d_0to3[0], "d3d of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_d3d_1", sv_d3d_0to3[1], "d3d of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_d3d_2", sv_d3d_0to3[2], "d3d of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_d3d_3", sv_d3d_0to3[3], "d3d of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_d3dsig_0", sv_d3dsig_0to3[0], "d3dsig of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_d3dsig_1", sv_d3dsig_0to3[1], "d3dsig of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_d3dsig_2", sv_d3dsig_0to3[2], "d3dsig of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_d3dsig_3", sv_d3dsig_0to3[3], "d3dsig of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_costhetasvpv_0", sv_costhetasvpv_0to3[0], "costhetasvpv of the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_costhetasvpv_1", sv_costhetasvpv_0to3[1], "costhetasvpv of the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_costhetasvpv_2", sv_costhetasvpv_0to3[2], "costhetasvpv of the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_costhetasvpv_3", sv_costhetasvpv_0to3[3], "costhetasvpv of the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_phirel_0", sv_phirel_0to3[0], "DeltaPhi(sv, jet) for the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_phirel_1", sv_phirel_0to3[1], "DeltaPhi(sv, jet) for the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_phirel_2", sv_phirel_0to3[2], "DeltaPhi(sv, jet) for the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_phirel_3", sv_phirel_0to3[3], "DeltaPhi(sv, jet) for the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_ptrel_0", sv_ptrel_0to3[0], "pT relative to parent jet for the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_ptrel_1", sv_ptrel_0to3[1], "pT relative to parent jet for the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_ptrel_2", sv_ptrel_0to3[2], "pT relative to parent jet for the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_ptrel_3", sv_ptrel_0to3[3], "pT relative to parent jet for the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_deltaR_0", sv_deltaR_0to3[0], "dR from parent jet for the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_deltaR_1", sv_deltaR_0to3[1], "dR from parent jet for the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_deltaR_2", sv_deltaR_0to3[2], "dR from parent jet for the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_deltaR_3", sv_deltaR_0to3[3], "dR from parent jet for the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
    
  djTable->addColumn<float>("DeepJet_sv_enratio_0", sv_enratio_0to3[0], "energy relative to parent jet for the first SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_enratio_1", sv_enratio_0to3[1], "energy relative to parent jet for the second SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_enratio_2", sv_enratio_0to3[2], "energy relative to parent jet for the third SV", nanoaod::FlatTable::FloatColumn, 10);
  djTable->addColumn<float>("DeepJet_sv_enratio_3", sv_enratio_0to3[3], "energy relative to parent jet for the fourth SV", nanoaod::FlatTable::FloatColumn, 10);
  iEvent.put(std::move(djTable), nameDeepJet_);
  
    
    
  //iEvent.put(std::move(outCands));
}

template< typename T>
void JetConstituentTableProducerDeepJet_WIP<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  //desc.add<std::string>("name", "JetPFCands");
  //desc.add<std::string>("nameSV", "JetSV");
  desc.add<std::string>("nameDeepJet", "Jet");
  //desc.add<std::string>("idx_name", "candIdx");
  //desc.add<std::string>("idx_nameSV", "svIdx");
  desc.add<std::string>("idx_nameDeepJet", "djIdx");
  desc.add<double>("jet_radius", true);
  desc.add<bool>("readBtag", true);
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK8"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("slimmedSecondaryVertices"));
  desc.add<edm::InputTag>("shallow_tag_infos", edm::InputTag("pfDeepCSVTagInfos"));
  desc.add<edm::InputTag>("puppi_value_map", edm::InputTag("puppi"));
  desc.add<edm::InputTag>("vertex_associator", edm::InputTag("primaryVertexAssociation", "original"));
  desc.add<double>("min_candidate_pt", 0.95);
  desc.add<bool>("flip", false);
  desc.add<bool>("fallback_puppi_weight", false);
  desc.add<bool>("fallback_vertex_association", false);
  descriptions.addWithDefaultLabel(desc);
}

typedef JetConstituentTableProducerDeepJet_WIP<pat::Jet> PatJetConstituentTableProducerDeepJet_WIP;
//typedef JetConstituentTableProducer<reco::GenJet> GenJetConstituentTableProducer;

DEFINE_FWK_MODULE(PatJetConstituentTableProducerDeepJet_WIP);
//DEFINE_FWK_MODULE(GenJetConstituentTableProducer);
