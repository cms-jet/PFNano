import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addPFCands(process, runOnMC=False, saveAll=False, addAK4=False, addAK8=False, addAK15=False, saveAllGen=False):
    '''
        Add PFCands and/or jet-PFCand association tables to NanoAOD

        To control which PFCands are saved:
            - If saveAll is True, save all PFcands
            - Otherwise, save the PFCands corresponding to the addAK* options.

        To control which jet-PFCand tables are made, use the addAK* options.

    '''
    if not (saveAll or addAK4 or addAK8 or addAK15):
        raise ValueError("In call to addPFCands(), at least one save* or add* option must be True.")

    if saveAll and (addAK4 or addAK8 or addAK15):
        raise ValueError("In call to addPFCands(), you can't specify saveAll=True and addAK*=True at the same time.")

    process.customizedPFCandsTask = cms.Task()
    process.schedule.associate(process.customizedPFCandsTask)

    if saveAll:
        candInput = cms.InputTag("packedPFCandidates")
    else:
        candList = cms.VInputTag()
        if addAK4:
            process.finalJetsAK4Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                                    src = cms.InputTag("finalJets"),
                                                    cut = cms.string("")
                                                    )
            candList.append(cms.InputTag("finalJetsAK4Constituents", "constituents"))
            process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)

        if addAK8:
            process.finalJetsAK8Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                                    src = cms.InputTag("finalJetsAK8"),
                                                    cut = cms.string("")
                                                    )
            candList.append(cms.InputTag("finalJetsAK8Constituents", "constituents"))
            process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)

        if addAK15:
            process.finalJetsAK15Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                                    src = cms.InputTag("finalJetsAK15"),
                                                    cut = cms.string("")
                                                    )
            candList.append(cms.InputTag("finalJetsAK15Constituents", "constituents"))
            process.customizedPFCandsTask.add(process.finalJetsAK15Constituents)
        process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger", 
                                                        src = candList, 
                                                        skipNulls = cms.bool(True), 
                                                        warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituents")

    # Make constituent table producers
    process.customConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                        src = candInput,
                                                        cut = cms.string(""), #we should not filter after pruning
                                                        name = cms.string("PFCands"),
                                                        doc = cms.string("interesting particles from various jet collections"),
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

    if addAK4:
        process.customAK4ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                            #candidates = cms.InputTag("packedPFCandidates"),
                                                            candidates = candInput,
                                                            jets = cms.InputTag("finalJets"),
                                                            jet_radius = cms.double(0.4),
                                                            name = cms.string("JetPFCands"),
                                                            idx_name = cms.string("pFCandsIdx"),
                                                            nameSV = cms.string("JetSVs"),
                                                            idx_nameSV = cms.string("sVIdx"),
                                                            )
    if addAK8:
        process.customAK8ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                            candidates = candInput,
                                                            jets = cms.InputTag("finalJetsAK8"),
                                                            jet_radius = cms.double(0.8),
                                                            name = cms.string("FatJetPFCands"),
                                                            idx_name = cms.string("pFCandsIdx"),
                                                            nameSV = cms.string("FatJetSVs"),
                                                            idx_nameSV = cms.string("sVIdx"),
                                                            )

    if addAK15:
        process.customAK15ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                            candidates = candInput,
                                                            jets       = cms.InputTag("finalJetsAK15"),
                                                            jet_radius = cms.double(1.5),
                                                            name       = cms.string("FatJetAK15PFCands"),
                                                            idx_name = cms.string("pFCandsIdx"),
                                                            nameSV     = cms.string("JetSVsAK15"),
                                                            idx_nameSV = cms.string("sVIdx"),
                                                            )

    # Add constituents tables to customizedPFCandsTask
    if not saveAll:
        process.customizedPFCandsTask.add(process.finalJetsConstituents)
    process.customizedPFCandsTask.add(process.customConstituentsExtTable)
    if addAK4:
        process.customizedPFCandsTask.add(process.customAK4ConstituentsTable)
    if addAK8:
        process.customizedPFCandsTask.add(process.customAK8ConstituentsTable)
    if addAK15:
        process.customizedPFCandsTask.add(process.customAK15ConstituentsTable)

    if runOnMC:
        if addAK4:
            process.genJetsAK4Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                        src = cms.InputTag("slimmedGenJets"),
                                                        cut = cms.string("pt > 20")
                                                        )
        if addAK8:
            process.genJetsAK8Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                        src = cms.InputTag("slimmedGenJetsAK8"),
                                                        cut = cms.string("pt > 100.")
                                                        )
        if addAK15:
            process.genJetsAK15Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                        src = cms.InputTag("ak15GenJetsNoNu"), # "slimmedGenJetsAK15"
                                                        cut = cms.string("pt > 100")
                                                        )

        if saveAll or saveAllGen:
            genCandInput = cms.InputTag("packedGenParticles")
        else:
            genCandList = cms.VInputTag()
            if addAK4:
                genCandList.append(cms.InputTag("genJetsAK4Constituents", "constituents"))
            if addAK8:
                genCandList.append(cms.InputTag("genJetsAK8Constituents", "constituents"))
            if addAK15:
                genCandList.append(cms.InputTag("genJetsAK15Constituents", "constituents"))

            # genJetsConstituents = merged set of PFCands
            process.genJetsConstituents = cms.EDProducer("PackedGenParticlePtrMerger", 
                                                            src        = genCandList, 
                                                            skipNulls  = cms.bool(True), 
                                                            warnOnSkip = cms.bool(True))
            genCandInput =  cms.InputTag("genJetsConstituents")


        process.genJetsParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                         src = genCandInput,
                                                         cut = cms.string(""), #we should not filter after pruning
                                                         name= cms.string("GenCands"),
                                                         doc = cms.string("interesting gen particles from various jet collections"),
                                                         singleton = cms.bool(False), # the number of entries is variable
                                                         extension = cms.bool(False), # this is the main table for the AK8 constituents
                                                         variables = cms.PSet(CandVars)
                                                     )
        if addAK8:
            process.genAK8ConstituentsTable = cms.EDProducer("GenJetConstituentTableProducer",
                                                             candidates = genCandInput,
                                                             jets = cms.InputTag("genJetsAK8Constituents"), # Note: The name has "Constituents" in it, but these are the jets
                                                             name = cms.string("GenFatJetCands"),
                                                             nameSV = cms.string("GenFatJetSVs"),
                                                             idx_name = cms.string("pFCandsIdx"),
                                                             idx_nameSV = cms.string("sVIdx"),
                                                             readBtag = cms.bool(False))
        if addAK4:
            process.genAK4ConstituentsTable = cms.EDProducer("GenJetConstituentTableProducer",
                                                             candidates = genCandInput,
                                                             jets = cms.InputTag("genJetsAK4Constituents"), # Note: The name has "Constituents" in it, but these are the jets
                                                             name = cms.string("GenJetCands"),
                                                             nameSV = cms.string("GenJetSVs"),
                                                             idx_name = cms.string("pFCandsIdx"),
                                                             idx_nameSV = cms.string("sVIdx"),
                                                             readBtag = cms.bool(False))
        if addAK15:
            process.genAK15ConstituentsTable = cms.EDProducer("GenJetConstituentTableProducer",
                                                             candidates = genCandInput,
                                                             jets = cms.InputTag("genJetsAK15Constituents"), # Note: The name has "Constituents" in it, but these are the jets
                                                             name = cms.string("GenJetCandsAK15"),
                                                             nameSV = cms.string("GenJetSVsAK15"),
                                                             idx_name = cms.string("pFCandsIdx"),
                                                             idx_nameSV = cms.string("sVIdx"),
                                                             readBtag = cms.bool(False))

        # Add everything to customizedPFCandsTask
        if addAK4:
            process.customizedPFCandsTask.add(process.genJetsAK4Constituents) #Note: For gen need to add jets to the process to keep pt cuts.
        if addAK8:
            process.customizedPFCandsTask.add(process.genJetsAK8Constituents)
        if addAK15:
            process.customizedPFCandsTask.add(process.genJetsAK15Constituents)
        if not (saveAll or saveAllGen):
            process.customizedPFCandsTask.add(process.genJetsConstituents)
        process.customizedPFCandsTask.add(process.genJetsParticleTable)
        if addAK4:
            process.customizedPFCandsTask.add(process.genAK4ConstituentsTable)
        if addAK8:
            process.customizedPFCandsTask.add(process.genAK8ConstituentsTable)
        if addAK15:
            process.customizedPFCandsTask.add(process.genAK15ConstituentsTable)

    return process
