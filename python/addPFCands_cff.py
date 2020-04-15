import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *


def addPFCands(process, runOnMC=False, onlyAK4=False, onlyAK8=False):
    process.customizedPFCandsTask = cms.Task()
    process.schedule.associate(process.customizedPFCandsTask)

    addAK4 = not onlyAK8
    addAK8 = not onlyAK4
    jets = cms.PSet()
    if addAK4:
        jets.Jet = cms.PSet(
            src = cms.InputTag("linkedObjects","jets"),
            isPuppi = cms.bool(False),
            cut = cms.string('pt > 20'),
            )
    if addAK8:
        jets.FatJet = cms.PSet(
            src = cms.InputTag('updatedJetsAK8WithUserData'),
            isPuppi = cms.bool(True),
            cut = cms.string('pt > 170'),
            )

    process.jetConstituentsTable = cms.EDProducer("JetConstituentTableProducer",
                                                  jets = jets,
                                                  name = cms.string("PFCands"),
                                                  check_indices = cms.bool(False),  # for debugging
                                                  )

    process.jetConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("jetConstituentsTable", "constituents"),
        cut = cms.string(""), #we should not filter after pruning
        name = cms.string("PFCands"),
        doc = cms.string("interesting particles from jets"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(True), # this is the extension table for the AK8 constituents
        variables = cms.PSet(CandVars,
            puppiWeight = Var("puppiWeight()", float, doc="Puppi weight", precision=10),
            puppiWeightNoLep = Var("puppiWeightNoLep()", float, doc="Puppi weight removing leptons", precision=10),
            vtxChi2 = Var("?hasTrackDetails()?vertexChi2():-1", float, doc="vertex chi2", precision=10),
            trkChi2 = Var("?hasTrackDetails()?pseudoTrack().normalizedChi2():-1", float, doc="normalized trk chi2", precision=10),
            dz = Var("dz()", float, doc="pf dz", precision=10),
            dzErr = Var("?hasTrackDetails()?dzError():-1", float, doc="pf dz err", precision=10),
            d0 = Var("dxy()", float, doc="pf d0", precision=10),
            d0Err = Var("?hasTrackDetails()?dxyError():-1", float, doc="pf d0 err", precision=10),
            pvAssocQuality = Var("pvAssociationQuality()", int, doc="primary vertex association quality"),
            lostInnerHits = Var("lostInnerHits()", int, doc="lost inner hits"),
            trkQuality = Var("?hasTrackDetails()?pseudoTrack().qualityMask():0", int, doc="track quality mask"),
            )
        )

    process.customizedPFCandsTask.add(process.jetConstituentsTable, process.jetConstituentsExtTable)

    if addAK4:
        ext = getattr(process.jetTable, 'externalVariables', cms.PSet())
        ext.nPFCand = ExtVar(cms.InputTag("jetConstituentsTable", "JetNpfcand"), int, doc="number of PFCands stored in the PFCands table")
        process.jetTable.externalVariables = ext
    if addAK8:
        ext = getattr(process.fatJetTable, 'externalVariables', cms.PSet())
        ext.nPFCand = ExtVar(cms.InputTag("jetConstituentsTable", "FatJetNpfcand"), int, doc="number of PFCands stored in the PFCands table")
        process.fatJetTable.externalVariables = ext

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

        if addAK4:
            process.customizedPFCandsTask.add(process.genJetsAK4Constituents)
            process.customizedPFCandsTask.add(process.genJetsAK4ParticleTable)
        if addAK8:
            process.customizedPFCandsTask.add(process.genJetsAK8Constituents)
            process.customizedPFCandsTask.add(process.genJetsAK8ParticleTable)

    return process
