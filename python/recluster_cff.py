import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################################################################################
######### For AK8 PUPPI jets
finalJetsAK8Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                            src = cms.InputTag("updatedJetsAK8"),
                                            cut = cms.string("pt > 170.0")
                                            )
genJetsAK8Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                            src = cms.InputTag("slimmedGenJetsAK8"),
                                            cut = cms.string("pt > 100.0")
                                            )



##################### Tables for final output and docs ##########################
finalJetsAK8ConstituentsTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("finalJetsAK8Constituents", "constituents"),
    cut = cms.string(""), #we should not filter after pruning
    name= cms.string("PFCandsAK8"),
    doc = cms.string("interesting gen particles from AK8 jets"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the AK8 constituents
    variables = cms.PSet(CandVars,
                            puppiWeight = Var("puppiWeight()", float, doc="Puppi weight",precision=10),
                            puppiWeightNoLep = Var("puppiWeightNoLep()", float, doc="Puppi weight removing leptons",precision=10),
                            vtxChi2 = Var("?hasTrackDetails()?vertexChi2():-1", float, doc="vertex chi2",precision=10),
                            trkChi2 = Var("?hasTrackDetails()?pseudoTrack().normalizedChi2():-1", float, doc="normalized trk chi2", precision=10),
                            dz = Var("?hasTrackDetails()?dz():-1", float, doc="pf dz", precision=10),
                            dzErr = Var("?hasTrackDetails()?dzError():-1", float, doc="pf dz err", precision=10),
                            d0 = Var("?hasTrackDetails()?dxy():-1", float, doc="pf d0", precision=10),
                            d0Err = Var("?hasTrackDetails()?dxyError():-1", float, doc="pf d0 err", precision=10),
                            pvAssocQuality = Var("pvAssociationQuality()", int, doc="primary vertex association quality"),
                            lostInnerHits = Var("lostInnerHits()", int, doc="lost inner hits"),
                            trkQuality = Var("?hasTrackDetails()?pseudoTrack().qualityMask():0", int, doc="track quality mask"),
    )
)

##################### Tables for final output and docs ##########################
genJetsAK8ParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("genJetsAK8Constituents", "constituents"),
    cut = cms.string(""), #we should not filter after pruning
    name= cms.string("GenPartAK8"),
    doc = cms.string("interesting gen particles from AK8 jets"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the AK8 constituents
    variables = cms.PSet(CandVars
    )
)

##################################################################################
######### For AK4 CHS jets
finalJetsAK4Constituents = finalJetsAK8Constituents.clone( src = 'updatedJets', cut = 'pt>10.0' )
finalJetsAK4ConstituentsTable = finalJetsAK8ConstituentsTable.clone(
                                                                src = cms.InputTag("finalJetsAK4Constituents", "constituents"),
                                                                name= cms.string("PFCandsAK4"),
                                                                doc = cms.string("interesting gen particles from AK4 jets"),
                                                                )
genJetsAK4Constituents = genJetsAK8Constituents.clone(
                                            src = cms.InputTag("slimmedGenJets"),
                                            cut = cms.string("pt > 10.0")
                                            )
genJetsAK4ParticleTable = genJetsAK8ParticleTable.clone(
                                                    src = cms.InputTag("genJetsAK4Constituents", "constituents"),
                                                    name= cms.string("GenPartAK4"),
                                                    doc = cms.string("interesting gen particles from AK4 jets"),
                                                    )

jetReclusterSequence = cms.Sequence(finalJetsAK4Constituents+finalJetsAK8Constituents)
jetReclusterMCSequence = cms.Sequence(genJetsAK4Constituents+genJetsAK8Constituents)
jetReclusterTable = cms.Sequence(finalJetsAK4ConstituentsTable+finalJetsAK8ConstituentsTable)
jetReclusterMCTable = cms.Sequence(genJetsAK4ParticleTable+genJetsAK8ParticleTable)

