Where do the AK15 jets come from?

process.finalJetsAK15 = cms.EDFilter("PATJetRefSelector",
    cut = cms.string(''),
    src = cms.InputTag("ak15WithUserData")
)

process.ak15WithUserData = cms.EDProducer("PATJetUserDataEmbedder",
    src = cms.InputTag("selectedUpdatedPatJetsAK15"),
    userFloats = cms.PSet(

    ),
    userInts = cms.PSet(
        tightId = cms.InputTag("tightJetIdAK15Puppi"),
        tightIdLepVeto = cms.InputTag("tightJetIdLepVetoAK15Puppi")
    )
)

process.selectedUpdatedPatJetsAK15 = cms.EDFilter("PATJetSelector",
    cut = cms.string(''),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("updatedPatJetsTransientCorrectedAK15")
)

process.updatedPatJetsTransientCorrectedAK15 = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(
        cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTagsAK15"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15","probXcc"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15","probXbb"), cms.InputTag("pfDeepCSVJetTagsAK15","probbb"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15","probXqq"), 
        cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15","probQCDothers"), cms.InputTag("pfJetProbabilityBJetTagsAK15"), cms.InputTag("pfDeepCSVJetTagsAK15","probb")
    ),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsTransientCorrectedAK15")),
    jetSource = cms.InputTag("updatedPatJetsAK15"),
    printWarning = cms.bool(True),
    sort = cms.bool(True),
    tagInfoSources = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosAK15"), cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosAK15"), cms.InputTag("pfParticleNetTagInfosAK15"), cms.InputTag("pfDeepCSVTagInfosAK15")),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)

process.updatedPatJetsAK15 = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsAK15")),
    jetSource = cms.InputTag("selectedUpdatedPatJetsAK15WithDeepInfo"),
    printWarning = cms.bool(True),
    sort = cms.bool(True),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)

process.selectedUpdatedPatJetsAK15WithDeepInfo = cms.EDFilter("PATJetSelector",
    cut = cms.string(''),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("updatedPatJetsTransientCorrectedAK15WithDeepInfo")
)

process.updatedPatJetsTransientCorrectedAK15WithDeepInfo = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(True),
    discriminatorSources = cms.VInputTag(
        cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probTbqq"), cms.InputTag("pfMassIndependentDeepDoubleCvLV2JetTagsAK15WithDeepInfo","probHcc"), cms.InputTag("pfDeepCSVJetTagsAK15WithDeepInfo","probudsg"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probQCDbb"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probTbc"), 
        cms.InputTag("pfMassIndependentDeepDoubleCvLJetTagsAK15WithDeepInfo","probHcc"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probTbmu"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probTbq"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15WithDeepInfo","probQCDbb"), cms.InputTag("pfMassIndependentDeepDoubleBvLV2JetTagsAK15WithDeepInfo","probHbb"), 
        cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probZbb"), cms.InputTag("pfMassIndependentDeepDoubleCvBV2JetTagsAK15WithDeepInfo","probHcc"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probQCDcc"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probQCDothers"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15WithDeepInfo","probXbb"), 
        cms.InputTag("pfMassIndependentDeepDoubleCvBJetTagsAK15WithDeepInfo","probHcc"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15WithDeepInfo","probXqq"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probZqq"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15WithDeepInfo","probQCDc"), cms.InputTag("pfJetProbabilityBJetTagsAK15WithDeepInfo"), 
        cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probHbb"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probWcq"), cms.InputTag("pfDeepCSVJetTagsAK15WithDeepInfo","probc"), cms.InputTag("pfDeepCSVJetTagsAK15WithDeepInfo","probb"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probHqqqq"), 
        cms.InputTag("pfMassIndependentDeepDoubleBvLJetTagsAK15WithDeepInfo","probHbb"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probWqq"), cms.InputTag("pfDeepCSVJetTagsAK15WithDeepInfo","probbb"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15WithDeepInfo","probQCDothers"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probQCDc"), 
        cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probQCDb"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probZcc"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15WithDeepInfo","probQCDb"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probTbel"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probTbta"), 
        cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15WithDeepInfo","probXcc"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probHcc"), cms.InputTag("pfMassDecorrelatedParticleNetJetTagsAK15WithDeepInfo","probQCDcc"), cms.InputTag("pfParticleNetJetTagsAK15WithDeepInfo","probTbcq"), cms.InputTag("pfMassDecorrelatedParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","XqqvsQCD"), 
        cms.InputTag("pfParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","WvsQCD"), cms.InputTag("pfMassDecorrelatedParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","XbbvsQCD"), cms.InputTag("pfParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","ZbbvsQCD"), cms.InputTag("pfParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","H4qvsQCD"), cms.InputTag("pfParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","HbbvsQCD"), 
        cms.InputTag("pfParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","HccvsQCD"), cms.InputTag("pfParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","TvsQCD"), cms.InputTag("pfMassDecorrelatedParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","XccvsQCD"), cms.InputTag("pfParticleNetDiscriminatorsJetTagsAK15WithDeepInfo","ZvsQCD")
    ),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsTransientCorrectedAK15WithDeepInfo")),
    jetSource = cms.InputTag("updatedPatJetsAK15WithDeepInfo"),
    printWarning = cms.bool(True),
    sort = cms.bool(False),
    tagInfoSources = cms.VInputTag(
        cms.InputTag("pfParticleNetTagInfosAK15WithDeepInfo"), cms.InputTag("pfDeepDoubleXTagInfosAK15WithDeepInfo"), cms.InputTag("pfBoostedDoubleSVAK8TagInfosAK15WithDeepInfo"), cms.InputTag("pfImpactParameterAK8TagInfosAK15WithDeepInfo"), cms.InputTag("pfInclusiveSecondaryVertexFinderAK8TagInfosAK15WithDeepInfo"), 
        cms.InputTag("pfDeepCSVTagInfosAK15WithDeepInfo"), cms.InputTag("pfImpactParameterTagInfosAK15WithDeepInfo"), cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosAK15WithDeepInfo"), cms.InputTag("pfDeepDoubleXTagInfosAK15WithDeepInfo")
    ),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)

process.updatedPatJetsAK15WithDeepInfo = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsAK15WithDeepInfo")),
    jetSource = cms.InputTag("packedPatJetsAK15PFPuppiSoftDrop"),
    printWarning = cms.bool(False),
    sort = cms.bool(False),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)

process.packedPatJetsAK15PFPuppiSoftDrop = cms.EDProducer("JetSubstructurePacker",
    algoLabels = cms.vstring('SoftDrop'),
    algoTags = cms.VInputTag(cms.InputTag("selectedPatJetsAK15PFPuppiSoftDropPacked")),
    distMax = cms.double(1.5),
    fixDaughters = cms.bool(False),
    jetSrc = cms.InputTag("selectedPatJetsAK15PFPuppi")
)

process.selectedPatJetsAK15PFPuppi = cms.EDFilter("PATJetSelector",
    cut = cms.string('pt > 160.0 && abs(rapidity()) < 2.4'),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("patJetsAK15PFPuppi")
)

process.patJetsAK15PFPuppi = cms.EDProducer("PATJetProducer",
    JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationAK15PFPuppi"),
    JetPartonMapSource = cms.InputTag("patJetFlavourAssociationLegacyAK15PFPuppi"),
    addAssociatedTracks = cms.bool(False),
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addEfficiencies = cms.bool(False),
    addGenJetMatch = cms.bool(True),
    addGenPartonMatch = cms.bool(True),
    addJetCharge = cms.bool(False),
    addJetCorrFactors = cms.bool(True),
    addJetFlavourInfo = cms.bool(True),
    addJetID = cms.bool(False),
    addPartonJetMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    addTagInfos = cms.bool(True),
    discriminatorSources = cms.VInputTag(cms.InputTag("pfJetProbabilityBJetTagsAK15PFPuppi"), cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTagsAK15PFPuppi"), cms.InputTag("pfDeepCSVJetTagsAK15PFPuppi","probbb"), cms.InputTag("pfDeepCSVJetTagsAK15PFPuppi","probb")),
    efficiencies = cms.PSet(),
    embedGenJetMatch = cms.bool(True),
    embedGenPartonMatch = cms.bool(True),
    embedPFCandidates = cms.bool(False),
    genJetMatch = cms.InputTag("patJetGenJetMatchAK15PFPuppi"),
    genPartonMatch = cms.InputTag("patJetPartonMatchAK15PFPuppi"),
    getJetMCFlavour = cms.bool(True),
    jetChargeSource = cms.InputTag(""),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsAK15PFPuppi")),
    jetIDMap = cms.InputTag("ak4JetID"),
    jetSource = cms.InputTag("ak15PFJetsPuppi"),
    partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
    resolutions = cms.PSet(),
    tagInfoSources = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosAK15PFPuppi"), cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosAK15PFPuppi"), cms.InputTag("pfDeepCSVTagInfosAK15PFPuppi")),
    trackAssociationSource = cms.InputTag(""),
    useLegacyJetMCFlavour = cms.bool(False),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag(
                "", "ak15PFJetsPuppiSoftDropMass", "NjettinessAK15Puppi:tau1", "NjettinessAK15Puppi:tau2", "NjettinessAK15Puppi:tau3", 
                "ak15PFJetsPuppiSoftDropValueMap:nb1AK15PuppiSoftDropN2", "ak15PFJetsPuppiSoftDropValueMap:nb1AK15PuppiSoftDropN3"
            )
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)
process.ak15PFJetsPuppi = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(5.0),
    Rho_EtaMax = cms.double(4.4),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(True),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.0),
    jetAlgorithm = cms.string('AntiKt'),
    jetPtMin = cms.double(5.0),
    jetType = cms.string('PFJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    rParam = cms.double(1.5),
    src = cms.InputTag("puppi"),
    srcPVs = cms.InputTag(""),
    useDeterministicSeed = cms.bool(True),
    voronoiRfact = cms.double(-0.9)
)


# Where do the AK15 PFcands come from?
process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger",
    skipNulls = cms.bool(True),
    src = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents","constituents"), cms.InputTag("finalJetsAK8Constituents","constituents"), cms.InputTag("finalJetsAK15Constituents","constituents")),
    warnOnSkip = cms.bool(True)
)

process.finalJetsAK15Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
    cut = cms.string(''),
    src = cms.InputTag("finalJetsAK15")
)

# ... so these should come from 'puppi'
