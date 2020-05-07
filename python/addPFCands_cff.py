import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addPFCands(process, runOnMC=False, onlyAK4=False, onlyAK8=False):
    process.customizedPFCandsTask = cms.Task( )
    process.schedule.associate(process.customizedPFCandsTask)

    process.customAK8ConstituentsTable = cms.EDProducer("JetConstituentTableProducer",
                                                        src=cms.InputTag("finalJetsAK8"),
                                                        jet_radius=cms.double(0.8),
                                                        cut=cms.string("pt()>170"),
                                                        namePF=cms.string("FatJetPFCands"),
                                                        nameSV=cms.string("FatJetSV"))

    process.customAK8ConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                       src = cms.InputTag("customAK8ConstituentsTable"),
                                       cut = cms.string(""), #we should not filter after pruning
                                       name = cms.string("FatJetPFCands"),
                                       doc = cms.string("interesting particles from AK8 jets"),
                                       singleton = cms.bool(False), # the number of entries is variable
                                       extension = cms.bool(True), # this is the extension table for the AK8 constituents
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

    process.customAK4ConstituentsTable = process.customAK8ConstituentsTable.clone(
        src='finalJets', jet_radius=cms.double(0.4), cut='pt()>20', namePF='JetPFCands', nameSV="JetSV")
    process.customAK4ConstituentsExtTable = process.customAK8ConstituentsExtTable.clone(
        src='customAK4ConstituentsTable',
        name='JetPFCands',
        doc='interesting particles from AK4 jets')

    if not onlyAK4:
        process.customizedPFCandsTask.add(process.customAK8ConstituentsTable)
        process.customizedPFCandsTask.add(process.customAK8ConstituentsExtTable)
    if not onlyAK8:
        process.customizedPFCandsTask.add(process.customAK4ConstituentsTable)
        process.customizedPFCandsTask.add(process.customAK4ConstituentsExtTable)

    if runOnMC:

        process.genJetsAK8Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                    src = cms.InputTag("slimmedGenJetsAK8"),
                                                    cut = cms.string("pt > 100.0")
                                                    )

        process.genJetsAK8ParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src = cms.InputTag("genJetsAK8Constituents", "constituents"),
            cut = cms.string(""), #we should not filter after pruning
            name= cms.string("GenJetAK8Cands"),
            doc = cms.string("interesting gen particles from AK8 jets"),
            singleton = cms.bool(False), # the number of entries is variable
            extension = cms.bool(False), # this is the main table for the AK8 constituents
            variables = cms.PSet(CandVars
            )
        )

        process.genJetsAK4Constituents = process.genJetsAK8Constituents.clone(
                                                    src = cms.InputTag("slimmedGenJets"),
                                                    cut = cms.string("pt > 20.0")
                                                    )
        process.genJetsAK4ParticleTable = process.genJetsAK8ParticleTable.clone(
                                                            src = cms.InputTag("genJetsAK4Constituents", "constituents"),
                                                            name= cms.string("GenJetCands"),
                                                            doc = cms.string("interesting gen particles from AK4 jets"),
                                                            )

        if not onlyAK8:
            process.customizedPFCandsTask.add(process.genJetsAK4Constituents)
            process.customizedPFCandsTask.add(process.genJetsAK4ParticleTable)
        if not onlyAK4:
            process.customizedPFCandsTask.add(process.genJetsAK8Constituents)
            process.customizedPFCandsTask.add(process.genJetsAK8ParticleTable)

    return process
