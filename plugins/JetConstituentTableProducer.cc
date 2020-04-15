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

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <unordered_map>
#include <unordered_set>

class JetConstituentTableProducer : public edm::stream::EDProducer<> {
public:
  explicit JetConstituentTableProducer(const edm::ParameterSet &);
  ~JetConstituentTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  typedef edm::Ptr<pat::PackedCandidate> CandidatePtr;
  typedef CandidatePtr::key_type Key;
  typedef edm::View<pat::PackedCandidate> CandidateView;
  typedef reco::VertexCollection VertexCollection;

  void produce(edm::Event &, const edm::EventSetup &) override;

  template <typename T>
  T getval(const std::unordered_map<Key, T> &m, Key key, T fallback = 0) {
    auto it = m.find(key);
    return it == m.end() ? fallback : it->second;
  }

  std::vector<CandidatePtr> getDaughters(const pat::Jet &jet, bool isPuppi);

  const std::string name_;
  const bool check_indices_;

  unsigned ncols_ = 0;
  std::vector<std::string> jet_names_;
  std::vector<bool> jet_ispuppi_;
  std::vector<StringCutObjectSelector<pat::Jet>> jet_cuts_;
  std::vector<edm::EDGetTokenT<edm::View<pat::Jet>>> jet_tokens_;
  std::vector<std::string> npf_valuemap_names_;
  std::vector<std::string> xref_table_names_;

  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<CandidateView> pfcand_token_;

  std::vector<edm::Handle<edm::View<pat::Jet>>> jets_;
  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<CandidateView> pfcands_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;
};

//
// constructors and destructor
//
JetConstituentTableProducer::JetConstituentTableProducer(const edm::ParameterSet &iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      check_indices_(iConfig.getParameter<bool>("check_indices")),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      pfcand_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pf_candidates"))) {
  const auto &pset = iConfig.getParameterSet("jets");
  jet_names_ = pset.getParameterNames();
  ncols_ = jet_names_.size();
  for (const auto &jetname : jet_names_) {
    const auto &p = pset.getParameterSet(jetname);
    jet_ispuppi_.emplace_back(p.getParameter<bool>("isPuppi"));
    jet_cuts_.emplace_back(p.getParameter<std::string>("cut"), true);
    jet_tokens_.emplace_back(consumes<edm::View<pat::Jet>>(p.getParameter<edm::InputTag>("src")));
  }

  produces<nanoaod::FlatTable>(name_);
  produces<std::vector<CandidatePtr>>("constituents");
  for (const auto &jetname : jet_names_) {
    produces<edm::ValueMap<int>>(npf_valuemap_names_.emplace_back(jetname + "Npfcand"));
    produces<nanoaod::FlatTable>(xref_table_names_.emplace_back(jetname + "To" + name_));
  }
}

JetConstituentTableProducer::~JetConstituentTableProducer() {}

void JetConstituentTableProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  auto outCands = std::make_unique<std::vector<CandidatePtr>>();
  std::vector<std::vector<int>> nPFCands(ncols_);
  std::vector<std::vector<int>> dauIdx(ncols_);  // for each vector: size = sum(nPFcand) for all jets in the event
  std::vector<std::vector<float>> btagEtaRel(ncols_), btagPtRatio(ncols_), btagPParRatio(ncols_), btagSip3dVal(ncols_),
      btagSip3dSig(ncols_), btagJetDistVal(ncols_);  // size = outCands.size()

  for (unsigned ic = 0; ic < ncols_; ++ic) {
    jets_.emplace_back();
    iEvent.getByToken(jet_tokens_[ic], jets_[ic]);
    // nPFCands must have same size as the jet collection in all cases
    nPFCands[ic] = std::vector<int>(jets_[ic]->size(), 0);
  }

  iEvent.getByToken(vtx_token_, vtxs_);
  if (!vtxs_->empty()) {
    iEvent.getByToken(pfcand_token_, pfcands_);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);

    std::vector<Key> kvec;                   // stores cand keys in an ordered vector
    std::unordered_set<Key> kset;            // used to check for duplicates -- not adding to kvec if already in
    std::unordered_map<Key, unsigned> kmap;  // map from cand's key to its index in the kvec

    // vmaps below: one for each jet collection
    std::vector<std::vector<Key>> vvec_dauKeys(ncols_);
    std::vector<std::unordered_map<Key, float>> vmap_btagEtaRel(ncols_), vmap_btagPtRatio(ncols_),
        vmap_btagPParRatio(ncols_), vmap_btagSip3dVal(ncols_), vmap_btagSip3dSig(ncols_), vmap_btagJetDistVal(ncols_);

    for (unsigned ic = 0; ic < ncols_; ++ic) {
      auto &v_dauKeys = vvec_dauKeys[ic];
      auto &m_btagEtaRel = vmap_btagEtaRel[ic];
      auto &m_btagPtRatio = vmap_btagPtRatio[ic];
      auto &m_btagPParRatio = vmap_btagPParRatio[ic];
      auto &m_btagSip3dVal = vmap_btagSip3dVal[ic];
      auto &m_btagSip3dSig = vmap_btagSip3dSig[ic];
      auto &m_btagJetDistVal = vmap_btagJetDistVal[ic];

      for (unsigned i_jet = 0; i_jet < jets_[ic]->size(); ++i_jet) {
        const auto &jet = jets_[ic]->at(i_jet);
        if (!jet_cuts_[ic](jet)) {
          continue;
        }

        math::XYZVector jet_dir = jet.momentum().Unit();
        GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());

        auto daughters = getDaughters(jet, jet_ispuppi_[ic]);
        nPFCands[ic][i_jet] = daughters.size();
        for (const auto &cand : daughters) {
          auto result = kset.insert(cand.key());
          if (result.second) {
            // if cand.key is not in kset
            kmap[cand.key()] = kvec.size();  // new val will be added to the end of kvec
            kvec.push_back(cand.key());
          }
          v_dauKeys.push_back(cand.key());  // add for the current jet collection
          if (cand->hasTrackDetails()) {
            btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
            trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, vtxs_->at(0));
            m_btagEtaRel[cand.key()] = trkinfo.getTrackEtaRel();
            m_btagPtRatio[cand.key()] = trkinfo.getTrackPtRatio();
            m_btagPParRatio[cand.key()] = trkinfo.getTrackPParRatio();
            m_btagSip3dVal[cand.key()] = trkinfo.getTrackSip3dVal();
            m_btagSip3dSig[cand.key()] = trkinfo.getTrackSip3dSig();
            m_btagJetDistVal[cand.key()] = trkinfo.getTrackJetDistVal();
          } else {
            m_btagEtaRel[cand.key()] = 0;
            m_btagPtRatio[cand.key()] = 0;
            m_btagPParRatio[cand.key()] = 0;
            m_btagSip3dVal[cand.key()] = 0;
            m_btagSip3dSig[cand.key()] = 0;
            m_btagJetDistVal[cand.key()] = 0;
          }
        }
      }  // end jet loop
    }

    // use pfcand key to determine the index in the output pfcand collection/table
    for (unsigned ic = 0; ic < ncols_; ++ic) {
      for (const auto &key : vvec_dauKeys[ic]) {
        dauIdx[ic].push_back(kmap.at(key));
      }
    }

    for (const auto &key : kvec) {
      outCands->push_back(pfcands_->ptrAt(key));
      // fill variables
      for (unsigned ic = 0; ic < ncols_; ++ic) {
        // fill 0 if this pfcand does not appear in this jet collection
        btagEtaRel[ic].push_back(getval(vmap_btagEtaRel[ic], key));
        btagPtRatio[ic].push_back(getval(vmap_btagPtRatio[ic], key));
        btagPParRatio[ic].push_back(getval(vmap_btagPParRatio[ic], key));
        btagSip3dVal[ic].push_back(getval(vmap_btagSip3dVal[ic], key));
        btagSip3dSig[ic].push_back(getval(vmap_btagSip3dSig[ic], key));
        btagJetDistVal[ic].push_back(getval(vmap_btagJetDistVal[ic], key));
      }
    }

    // for testing the implementation -- turn if off for production
    if (check_indices_) {
      for (unsigned ic = 0; ic < ncols_; ++ic) {
        unsigned idau_global = 0;
        for (unsigned i_jet = 0; i_jet < jets_[ic]->size(); ++i_jet) {
          const auto &jet = jets_[ic]->at(i_jet);
          if (nPFCands[ic][i_jet] == 0)
            continue;
          auto daughters = getDaughters(jet, jet_ispuppi_[ic]);
          for (unsigned idau = 0; idau < daughters.size(); ++idau) {
            if (daughters[idau].key() != outCands->at(dauIdx[ic].at(idau_global)).key()) {
              throw cms::Exception("RuntimeError") << "Inconsistent index for jet collection " << jet_names_[ic];
            }
            ++idau_global;
          }
        }
      }
    }  // check_indices_
  }

  auto candTable = std::make_unique<nanoaod::FlatTable>(outCands->size(), name_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTableProducer
  for (unsigned ic = 0; ic < ncols_; ++ic) {
    auto suffix = "_" + jet_names_[ic];
    auto doc = "(" + jet_names_[ic] + ")";
    candTable->addColumn<float>(
        "btagEtaRel" + suffix, btagEtaRel[ic], "btagEtaRel" + doc, nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>(
        "btagPtRatio" + suffix, btagPtRatio[ic], "btagPtRatio" + doc, nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>(
        "btagPParRatio" + suffix, btagPParRatio[ic], "btagPParRatio" + doc, nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>(
        "btagSip3dVal" + suffix, btagSip3dVal[ic], "btagSip3dVal" + doc, nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>(
        "btagSip3dSig" + suffix, btagSip3dSig[ic], "btagSip3dSig" + doc, nanoaod::FlatTable::FloatColumn, 10);
    candTable->addColumn<float>(
        "btagJetDistVal" + suffix, btagJetDistVal[ic], "btagJetDistVal" + doc, nanoaod::FlatTable::FloatColumn, 10);
  }

  iEvent.put(std::move(candTable), name_);
  iEvent.put(std::move(outCands), "constituents");
  for (unsigned ic = 0; ic < ncols_; ++ic) {
    auto npfMap = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler filler(*npfMap);
    filler.insert(jets_[ic], nPFCands[ic].begin(), nPFCands[ic].end());
    filler.fill();
    iEvent.put(std::move(npfMap), npf_valuemap_names_[ic]);

    auto xrefTable = std::make_unique<nanoaod::FlatTable>(dauIdx[ic].size(), xref_table_names_[ic], false);
    xrefTable->addColumn<int>("candIdx",
                              dauIdx[ic],
                              "Indices of the jet constitutes in the PFCand table. Use nPFCands in the jet table to "
                              "separate these indices for each jet.",
                              nanoaod::FlatTable::IntColumn);
    iEvent.put(std::move(xrefTable), xref_table_names_[ic]);
  }
}

std::vector<JetConstituentTableProducer::CandidatePtr> JetConstituentTableProducer::getDaughters(const pat::Jet &jet,
                                                                                                 bool isPuppi) {
  std::vector<CandidatePtr> daughters;
  for (const auto &cand : jet.daughterPtrVector()) {
    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    assert(packed_cand != nullptr);
    // remove particles w/ extremely low puppi weights (needed for 2017 MiniAOD)
    if (isPuppi && packed_cand->puppiWeight() < 0.01)
      continue;
    // get the original reco/packed candidate not scaled by the puppi weight
    daughters.push_back(pfcands_->ptrAt(cand.key()));
  }
  // sort by original pt (not Puppi weighted)
  std::sort(daughters.begin(), daughters.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });
  return daughters;
}

void JetConstituentTableProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("name", "JetPFCands");
  desc.add<bool>("check_indices", false);
  edm::ParameterSetDescription jets;
  jets.setAllowAnything();
  desc.add<edm::ParameterSetDescription>("jets", jets);
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("pf_candidates", edm::InputTag("packedPFCandidates"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(JetConstituentTableProducer);
