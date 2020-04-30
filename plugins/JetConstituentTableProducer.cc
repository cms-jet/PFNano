#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

template<typename T>
class JetConstituentTableProducer : public edm::stream::EDProducer<> {
public:
  explicit JetConstituentTableProducer(const edm::ParameterSet &);
  ~JetConstituentTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  typedef reco::VertexCollection VertexCollection;

  const std::string name_;
  const bool readBtag_;

  edm::EDGetTokenT<edm::View<T>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<reco::CandidateView> cand_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<reco::CandidateView> cands_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;
};

//
// constructors and destructor
//
template< typename T>
JetConstituentTableProducer<T>::JetConstituentTableProducer(const edm::ParameterSet &iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      readBtag_(iConfig.getParameter<bool>("readBtag")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))) {
  produces<nanoaod::FlatTable>(name_);
  produces<std::vector<reco::CandidatePtr>>();
}

template< typename T>
JetConstituentTableProducer<T>::~JetConstituentTableProducer() {}

template< typename T>
void JetConstituentTableProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  auto outCands = std::make_unique<std::vector<reco::CandidatePtr>>();
  std::vector<int> jetIdx, candIdx;
  std::vector<float> btagEtaRel, btagPtRatio, btagPParRatio, btagSip3dVal, btagSip3dSig, btagJetDistVal;
  auto jets = iEvent.getHandle(jet_token_);
  iEvent.getByToken(cand_token_, cands_);

  if(readBtag_){
    iEvent.getByToken(vtx_token_, vtxs_);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);
  }

  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet = jets->at(i_jet);
    math::XYZVector jet_dir = jet.momentum().Unit();
    GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
    
    std::vector<reco::CandidatePtr> const & daughters = jet.daughterPtrVector();

    for (const auto &cand : daughters) {
      auto candPtrs = cands_->ptrs();
      auto candInNewList = std::find( candPtrs.begin(), candPtrs.end(), cand );
      if ( candInNewList == candPtrs.end() ) {
	std::cout << "Cannot find candidate : " << cand.id() << ", " << cand.key() << ", pt = " << cand->pt() << std::endl;
	continue;
      }
      outCands->push_back(cand);
      jetIdx.push_back(i_jet);
      candIdx.push_back( candInNewList - candPtrs.begin() );
      if (readBtag_ && !vtxs_->empty()) {
	if ( cand.isNull() ) continue;
	auto const *packedCand = dynamic_cast <pat::PackedCandidate const *>(cand.get());
	if ( packedCand == nullptr ) continue;
	if ( packedCand && packedCand->hasTrackDetails()){
	  btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
	  trkinfo.buildTrackInfo(&(*packedCand), jet_dir, jet_ref_track_dir, vtxs_->at(0));
	  btagEtaRel.push_back(trkinfo.getTrackEtaRel());
	  btagPtRatio.push_back(trkinfo.getTrackPtRatio());
	  btagPParRatio.push_back(trkinfo.getTrackPParRatio());
	  btagSip3dVal.push_back(trkinfo.getTrackSip3dVal());
	  btagSip3dSig.push_back(trkinfo.getTrackSip3dSig());
	  btagJetDistVal.push_back(trkinfo.getTrackJetDistVal());
	} else {
          btagEtaRel.push_back(0);
          btagPtRatio.push_back(0);
          btagPParRatio.push_back(0);
          btagSip3dVal.push_back(0);
          btagSip3dSig.push_back(0);
          btagJetDistVal.push_back(0);
        }
      }
    }  // end jet loop
  }

  auto candTable = std::make_unique<nanoaod::FlatTable>(outCands->size(), name_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTableProducer
  candTable->addColumn<int>("candIdx", candIdx, "Index in the candidate list", nanoaod::FlatTable::IntColumn);
  candTable->addColumn<int>("jetIdx", jetIdx, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
  if (readBtag_) {
    candTable->addColumn<float>("btagEtaRel", btagEtaRel, "btagEtaRel", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagPtRatio", btagPtRatio, "btagPtRatio", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagPParRatio", btagPParRatio, "btagPParRatio", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagSip3dVal", btagSip3dVal, "btagSip3dVal", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagSip3dSig", btagSip3dSig, "btagSip3dSig", nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>("btagJetDistVal", btagJetDistVal, "btagJetDistVal", nanoaod::FlatTable::FloatColumn, 10);
  }
  iEvent.put(std::move(candTable), name_);
  iEvent.put(std::move(outCands));
}

template< typename T>
void JetConstituentTableProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("name", "JetPFCands");
  desc.add<bool>("readBtag", true);
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK8"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  descriptions.addWithDefaultLabel(desc);
}

typedef JetConstituentTableProducer<pat::Jet> PatJetConstituentTableProducer;
typedef JetConstituentTableProducer<reco::GenJet> GenJetConstituentTableProducer;

DEFINE_FWK_MODULE(PatJetConstituentTableProducer);
DEFINE_FWK_MODULE(GenJetConstituentTableProducer);
