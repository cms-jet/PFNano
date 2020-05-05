#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// SV
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
using namespace btagbtvdeep;

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class JetConstituentTableProducer : public edm::stream::EDProducer<> {
public:
  explicit JetConstituentTableProducer(const edm::ParameterSet &);
  ~JetConstituentTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  typedef edm::Ptr<pat::PackedCandidate> CandidatePtr;
  typedef edm::View<pat::PackedCandidate> CandidateView;
  typedef reco::VertexCollection VertexCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;

  const std::string name_;
  const StringCutObjectSelector<pat::Jet> jetCut_;
  const double jet_radius_;

  edm::EDGetTokenT<edm::View<pat::Jet>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<CandidateView> pfcand_token_;
  edm::EDGetTokenT<SVCollection> sv_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<CandidateView> pfcands_;
  edm::Handle<SVCollection> svs_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;

  const reco::Vertex *pv_ = nullptr;


};

//
// constructors and destructor
//
JetConstituentTableProducer::JetConstituentTableProducer(const edm::ParameterSet &iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      jetCut_(iConfig.getParameter<std::string>("cut"), true),
      jet_radius_(iConfig.getParameter<double>("jet_radius")),
      jet_token_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      pfcand_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pf_candidates"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))){
  produces<nanoaod::FlatTable>(name_);
  produces<std::vector<CandidatePtr>>();
}

JetConstituentTableProducer::~JetConstituentTableProducer() {}


void JetConstituentTableProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  auto outCands = std::make_unique<std::vector<CandidatePtr>>();
  std::vector<int> jetIdx;
  std::vector<float> btagEtaRel, btagPtRatio, btagPParRatio, btagSip3dVal, btagSip3dSig, btagJetDistVal;

  iEvent.getByToken(vtx_token_, vtxs_);
  if (!vtxs_->empty()) {
    pv_ = &vtxs_->at(0);
    auto jets = iEvent.getHandle(jet_token_);
    iEvent.getByToken(pfcand_token_, pfcands_);
    // edm::Handle<SVCollection> svs;
    iEvent.getByToken(sv_token_, svs_);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);

    for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
      const auto &jet = jets->at(i_jet);

      std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " << std::endl;
      //reco::SecondaryVertexTagInfo const* svinfo = jet.tagInfoSecondaryVertex("pfInclusiveSecondaryVertexFinderTagInfos");
      //reco::SecondaryVertexTagInfo const* svinfo = jet.tagInfoSecondaryVertex();
      //reco::CandSecondaryVertexTagInfo const* svinfo = jet.tagInfoSecondaryVertex();

      //auto svinfo = jet.tagInfoSecondaryVertex("svTagInfos");
      // std::cout << "A" << std::endl;
      // if (svinfo != nullptr) {
      //   int nSVs = svinfo->nVertices();
      //   //int nSVs = svinfo.nVertices();
      //   std::cout << "XX " << nSVs << std::endl;
      // }
      // else {
      //   std::cout << "FUCK"  << std::endl;
      // }

      std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
      for (const auto &sv : *svs_) {
        if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
          std::cout << "FUCK" << sv.pt()  << std::endl;
          jetSVs.push_back(&sv);
        }
      }
      // sort by dxy significance
      std::sort(jetSVs.begin(),
            jetSVs.end(),
            [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) {
              return sv_vertex_comparator(*sva, *svb, *pv_);
            });

      if (!jetCut_(jet))
        continue;

      math::XYZVector jet_dir = jet.momentum().Unit();
      GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());

      std::vector<CandidatePtr> daughters;
      for (const auto &cand : jet.daughterPtrVector()) {
        const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
        assert(packed_cand != nullptr);
        // remove particles w/ extremely low puppi weights (needed for 2017 MiniAOD)
        if (packed_cand->puppiWeight() < 0.01)
          continue;
        // get the original reco/packed candidate not scaled by the puppi weight
        daughters.push_back(pfcands_->ptrAt(cand.key()));
      }
      // sort by original pt (not Puppi-weighted)
      std::sort(daughters.begin(), daughters.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });

      for (const auto &cand : daughters) {
        outCands->push_back(cand);
        jetIdx.push_back(i_jet);
        if (cand->hasTrackDetails()){
          btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
          trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, vtxs_->at(0));
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
  candTable->addColumn<int>("jetIdx", jetIdx, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
  candTable->addColumn<float>("btagEtaRel", btagEtaRel, "btagEtaRel", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("btagPtRatio", btagPtRatio, "btagPtRatio", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("btagPParRatio", btagPParRatio, "btagPParRatio", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("btagSip3dVal", btagSip3dVal, "btagSip3dVal", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("btagSip3dSig", btagSip3dSig, "btagSip3dSig", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("btagJetDistVal", btagJetDistVal, "btagJetDistVal", nanoaod::FlatTable::FloatColumn, 10);

  iEvent.put(std::move(candTable), name_);
  iEvent.put(std::move(outCands));
}

void JetConstituentTableProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("slimmedJetsAK8"));
  desc.add<std::string>("name", "AK8PFCands");
  desc.add<std::string>("cut", "pt()>170");
  desc.add<double>("jet_radius", 0.4);
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("pf_candidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("slimmedSecondaryVertices"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(JetConstituentTableProducer);
