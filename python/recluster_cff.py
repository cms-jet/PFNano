import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *





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

jetReclusterSequence = cms.Sequence(finalJetsAK8Constituents)
jetReclusterMCSequence = cms.Sequence(genJetsAK8Constituents)
jetReclusterTable = cms.Sequence(finalJetsAK8ConstituentsTable)
jetReclusterMCTable = cms.Sequence(genJetsAK8ParticleTable)

