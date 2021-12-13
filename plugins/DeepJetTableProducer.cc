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

// anstein: add tag info and a way to go back to the jet reference
#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/DeepFlavourTagInfo.h"

template<typename T>
class DeepJetTableProducer : public edm::stream::EDProducer<> {
public:
  explicit DeepJetTableProducer(const edm::ParameterSet &);
  ~DeepJetTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;
    
  const std::string nameDeepJet_;
  const std::string idx_nameDeepJet_;

  edm::EDGetTokenT<edm::View<T>> jet_token_;
    
  // anstein: from ONNX producer
  typedef std::vector<reco::DeepFlavourTagInfo> TagInfoCollection;
  const edm::EDGetTokenT<TagInfoCollection> tag_info_src_;
  // anstein: could be used later to manually change number of features and constituents per feature
  /*
  constexpr static unsigned n_features_global_ = 15;
  constexpr static unsigned n_cpf_ = 25;
  constexpr static unsigned n_features_cpf_ = 16;
  constexpr static unsigned n_npf_ = 25;
  constexpr static unsigned n_features_npf_ = 6;
  constexpr static unsigned n_sv_ = 4;
  constexpr static unsigned n_features_sv_ = 12;
  */
};

//
// constructors and destructor
//
template< typename T>
DeepJetTableProducer<T>::DeepJetTableProducer(const edm::ParameterSet &iConfig)
    : nameDeepJet_(iConfig.getParameter<std::string>("nameDeepJet")),
      idx_nameDeepJet_(iConfig.getParameter<std::string>("idx_nameDeepJet")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      tag_info_src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("tagInfo_src"))){
  produces<nanoaod::FlatTable>(nameDeepJet_);
}

template< typename T>
DeepJetTableProducer<T>::~DeepJetTableProducer() {}

template< typename T>
void DeepJetTableProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  // anstein: only necessary to explicitly check correct matching of jets
  //  std::vector<int> jetIdx_dj;
  
  auto jets = iEvent.getHandle(jet_token_);

  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(tag_info_src_, tag_infos);
    
  
  std::vector<int> jet_N_CPFCands(jets->size());
  std::vector<int> jet_N_NPFCands(jets->size());
  std::vector<int> jet_N_PVs(jets->size()); 
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
    
    
  // should default to 0 if less than 25 npf with information
  std::vector<std::vector<float>> Npfcan_ptrel_0to24(25, std::vector<float>(jets->size()));
  std::vector<std::vector<float>> Npfcan_deltaR_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<int>> Npfcan_isGamma_0to24(25, std::vector<int>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_HadFrac_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_drminsv_0to24(25, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> Npfcan_puppiw_0to24(25, std::vector<float>(jets->size()));
    
    
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
  /*
  // anstein: not relevant for this version of the tagger
  std::vector<std::vector<float>> sv_ptrel_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_phirel_0to3(4, std::vector<float>(jets->size()));
  */
  std::vector<std::vector<float>> sv_deltaR_0to3(4, std::vector<float>(jets->size())); 
  std::vector<std::vector<float>> sv_enratio_0to3(4, std::vector<float>(jets->size())); 
  
  if (!tag_infos->empty()) {  

      for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
          // anstein: new version of the jet loop which reads tag info instead of contituent info

          const auto& taginfo = (*tag_infos)[i_jet];
          const auto& features = taginfo.features();
          // jet variables
          //const auto& jet_features = features.jet_features;
          // anstein: already part of jet table
          /*
          *ptr = jet_features.pt;
          *(++ptr) = jet_features.eta;
          */
          // number of elements in different collections
          jet_N_CPFCands[i_jet] = features.c_pf_features.size();
          jet_N_NPFCands[i_jet] = features.n_pf_features.size();
          jet_N_SVs[i_jet] = features.sv_features.size();
          jet_N_PVs[i_jet] = features.npv;
          // variables from ShallowTagInfo
          // anstein: already included via DeepCSV
          /*
          const auto& tag_info_features = features.tag_info_features;
          *(++ptr) = tag_info_features.trackSumJetEtRatio;
          *(++ptr) = tag_info_features.trackSumJetDeltaR;
          *(++ptr) = tag_info_features.vertexCategory;
          *(++ptr) = tag_info_features.trackSip2dValAboveCharm;
          *(++ptr) = tag_info_features.trackSip2dSigAboveCharm;
          *(++ptr) = tag_info_features.trackSip3dValAboveCharm;
          *(++ptr) = tag_info_features.trackSip3dSigAboveCharm;
          *(++ptr) = tag_info_features.jetNSelectedTracks;
          *(++ptr) = tag_info_features.jetNTracksEtaRel;
          */
          
          // c_pf candidates
          auto max_c_pf_n = std::min(features.c_pf_features.size(), (std::size_t)25);
          // anstein: I don't need the offset (similar for npf and sv)
          //offset = i_jet * input_sizes_[kChargedCandidates];
          for (std::size_t c_pf_n = 0; c_pf_n < max_c_pf_n; c_pf_n++) {
            const auto& c_pf_features = features.c_pf_features.at(c_pf_n);
            Cpfcan_BtagPf_trackEtaRel_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackEtaRel;
            Cpfcan_BtagPf_trackPtRel_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackPtRel;
            Cpfcan_BtagPf_trackPPar_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackPPar;
            Cpfcan_BtagPf_trackDeltaR_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackDeltaR;
            Cpfcan_BtagPf_trackPParRatio_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackPParRatio;
            Cpfcan_BtagPf_trackSip2dVal_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackSip2dVal;
            Cpfcan_BtagPf_trackSip2dSig_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackSip2dSig;
            Cpfcan_BtagPf_trackSip3dVal_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackSip3dVal;
            Cpfcan_BtagPf_trackSip3dSig_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackSip3dSig;
            Cpfcan_BtagPf_trackJetDistVal_0to24[c_pf_n][i_jet] = c_pf_features.btagPf_trackJetDistVal;
            Cpfcan_ptrel_0to24[c_pf_n][i_jet] = c_pf_features.ptrel;
            Cpfcan_drminsv_0to24[c_pf_n][i_jet] = c_pf_features.drminsv;
            Cpfcan_VTX_ass_0to24[c_pf_n][i_jet] = c_pf_features.vtx_ass;
            Cpfcan_puppiw_0to24[c_pf_n][i_jet] = c_pf_features.puppiw;
            Cpfcan_chi2_0to24[c_pf_n][i_jet] = c_pf_features.chi2;
            Cpfcan_quality_0to24[c_pf_n][i_jet] = c_pf_features.quality;
          }
          
          // n_pf candidates
          auto max_n_pf_n = std::min(features.n_pf_features.size(), (std::size_t)25);
          for (std::size_t n_pf_n = 0; n_pf_n < max_n_pf_n; n_pf_n++) {
            const auto& n_pf_features = features.n_pf_features.at(n_pf_n);
            Npfcan_ptrel_0to24[n_pf_n][i_jet] = n_pf_features.ptrel;
            Npfcan_deltaR_0to24[n_pf_n][i_jet] = n_pf_features.deltaR;
            Npfcan_isGamma_0to24[n_pf_n][i_jet] = n_pf_features.isGamma;
            Npfcan_HadFrac_0to24[n_pf_n][i_jet] = n_pf_features.hadFrac;
            Npfcan_drminsv_0to24[n_pf_n][i_jet] = n_pf_features.drminsv;
            Npfcan_puppiw_0to24[n_pf_n][i_jet] = n_pf_features.puppiw;
          }

          // sv candidates
          auto max_sv_n = std::min(features.sv_features.size(), (std::size_t)4);
          for (std::size_t sv_n = 0; sv_n < max_sv_n; sv_n++) {
            const auto& sv_features = features.sv_features.at(sv_n);
            sv_pt_0to3[sv_n][i_jet] = sv_features.pt;
            sv_deltaR_0to3[sv_n][i_jet] = sv_features.deltaR;
            sv_mass_0to3[sv_n][i_jet] = sv_features.mass;
            sv_ntracks_0to3[sv_n][i_jet] = sv_features.ntracks;
            sv_chi2_0to3[sv_n][i_jet] = sv_features.chi2;
            sv_normchi2_0to3[sv_n][i_jet] = sv_features.normchi2;
            sv_dxy_0to3[sv_n][i_jet] = sv_features.dxy;
            sv_dxysig_0to3[sv_n][i_jet] = sv_features.dxysig;
            sv_d3d_0to3[sv_n][i_jet] = sv_features.d3d;
            sv_d3dsig_0to3[sv_n][i_jet] = sv_features.d3dsig;
            sv_costhetasvpv_0to3[sv_n][i_jet] = sv_features.costhetasvpv;
            sv_enratio_0to3[sv_n][i_jet] = sv_features.enratio;
          }
      }
  }
    
  // DeepJetInputs table
  auto djTable = std::make_unique<nanoaod::FlatTable>(jet_N_CPFCands.size(), nameDeepJet_, false, true);
  //djTable->addColumn<int>("DeepJet_jetIdx", jetIdx_dj, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
    
    
  djTable->addColumn<int>("DeepJet_nCpfcand", jet_N_CPFCands, "Number of charged PF candidates in the jet", nanoaod::FlatTable::IntColumn);
  djTable->addColumn<int>("DeepJet_nNpfcand", jet_N_NPFCands, "Number of neutral PF candidates in the jet", nanoaod::FlatTable::IntColumn);
  djTable->addColumn<int>("DeepJet_nsv", jet_N_SVs, "Number of secondary vertices in the jet", nanoaod::FlatTable::IntColumn);
  djTable->addColumn<int>("DeepJet_npv", jet_N_PVs, "Number of primary vertices", nanoaod::FlatTable::IntColumn);
    
    
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
      
  }  

    
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
      /*
      // anstein: only relevant if also included in the tag info, not yet, maybe in future versions of the tagger
      input_name = "DeepJetExtra_sv_phirel_" + s;
      description = "DeltaPhi(sv, jet) for the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_phirel_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJetExtra_sv_ptrel_" + s;
      description = "pT relative to parent jet for the " + s + ". SV";
      djTable->addColumn<float>(input_name, sv_ptrel_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      */
      input_name = "DeepJet_sv_deltaR_" + s;
      description = "pseudoangular distance between jet axis and the " + s + ". SV direction";
      djTable->addColumn<float>(input_name, sv_deltaR_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "DeepJet_sv_enratio_" + s;
      description = "ratio of the " + s + ". SV energy ratio to the jet energy";
      djTable->addColumn<float>(input_name, sv_enratio_0to3[p], description, nanoaod::FlatTable::FloatColumn, 10);

  }
  
    
  iEvent.put(std::move(djTable), nameDeepJet_);  
  
}

template< typename T>
void DeepJetTableProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("nameDeepJet", "Jet");
  desc.add<std::string>("idx_nameDeepJet", "djIdx");
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("tagInfo_src", edm::InputTag("pfDeepFlavourTagInfosWithDeepInfo"));
  descriptions.addWithDefaultLabel(desc);
}

typedef DeepJetTableProducer<pat::Jet> PatJetDeepJetTableProducer;
//typedef DeepJetTableProducer<reco::GenJet> GenJetDeepJetTableProducer;

DEFINE_FWK_MODULE(PatJetDeepJetTableProducer);
//DEFINE_FWK_MODULE(GenJetDeepJetTableProducer);
