import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addPFCands(process, runOnMC=False, onlyAK4=False, onlyAK8=False):
    process.customizedPFCandsTask = cms.Task( )
    process.schedule.associate(process.customizedPFCandsTask)

    process.finalJetsAK8Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                            src = cms.InputTag("finalJetsAK8"),
                                            cut = cms.string("pt > 170.0")
                                            )
    process.finalJetsAK4Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                            src = cms.InputTag("finalJets"),
                                            cut = cms.string("pt > 30.0")
                                            )
    if onlyAK4:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
    elif onlyAK8:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK8Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)
    else:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"), cms.InputTag("finalJetsAK8Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
        process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)
    process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger",
                                                  src = candList)
    process.customConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                        src = cms.InputTag("finalJetsConstituents"),
                                                        cut = cms.string(""), #we should not filter after pruning
                                                        name = cms.string("JetPFCands"),
                                                        doc = cms.string("interesting particles from AK4 and AK8 jets"),
                                                        singleton = cms.bool(False), # the number of entries is variable
                                                        extension = cms.bool(False), # this is the extension table for the AK8 constituents
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
    
    process.customizedPFCandsTask.add(process.finalJetsConstituents)
    process.customizedPFCandsTask.add(process.customConstituentsExtTable)

    if runOnMC:

        process.genJetsAK8Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                    src = cms.InputTag("slimmedGenJetsAK8"),
                                                    cut = cms.string("pt > 100.0")
                                                    )

      
        process.genJetsAK4Constituents = process.genJetsAK8Constituents.clone(
                                                    src = cms.InputTag("slimmedGenJets"),
                                                    cut = cms.string("pt > 20.0")
                                                    )
        if onlyAK4:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK4Constituents", "constituents"))
            process.customizedPFCandsTask.add(process.genJetsAK4Constituents)
        elif onlyAK8:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK8Constituents", "constituents"))
            process.customizedPFCandsTask.add(process.genJetsAK8Constituents)
        else:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK4Constituents", "constituents"), cms.InputTag("genJetsAK8Constituents", "constituents"))
            process.customizedPFCandsTask.add(process.genJetsAK4Constituents)
            process.customizedPFCandsTask.add(process.genJetsAK8Constituents)
        process.genJetsConstituents = cms.EDProducer("PackedGenParticlePtrMerger",
                                                    src = genCandList
                                                    )
        process.genJetsParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                         src = cms.InputTag("genJetsConstituents"),
                                                         cut = cms.string(""), #we should not filter after pruning
                                                         name= cms.string("GenJetCands"),
                                                         doc = cms.string("interesting gen particles from AK4 and AK8 jets"),
                                                         singleton = cms.bool(False), # the number of entries is variable
                                                         extension = cms.bool(False), # this is the main table for the AK8 constituents
                                                         variables = cms.PSet(CandVars
                                                                          )
                                                     )
        process.customizedPFCandsTask.add(process.genJetsConstituents)
        process.customizedPFCandsTask.add(process.genJetsParticleTable)
      
    return process
