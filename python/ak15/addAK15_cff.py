import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_run2_jme_2016_cff import run2_jme_2016
from Configuration.Eras.Modifier_run2_jme_2017_cff import run2_jme_2017
from Configuration.Eras.Modifier_run2_nanoAOD_106Xv1_cff import run2_nanoAOD_106Xv1

from PhysicsTools.NanoAOD.common_cff import Var,P4Vars
from PhysicsTools.PFNano.addPFCands_cff import addPFCands
from PhysicsTools.PFNano.addBTV import add_BTV

# ---------------------------------------------------------

def setupAK15(process, runOnMC=False, path=None, runParticleNet=False, runParticleNetMD=False):
    # recluster Puppi jets
    bTagDiscriminators = [
        'pfJetProbabilityBJetTags',
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probbb',
    ]
    subjetBTagDiscriminators = [
        'pfJetProbabilityBJetTags',
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probbb',
    ]
    JETCorrLevels = ['L2Relative', 'L3Absolute', 'L2L3Residual']

    #print("Flag 1")
    #for aatt in dir(process):
    #    if "ak4" in aatt.lower():
    #        print(aatt)


    from PhysicsTools.PFNano.ak15.jetToolbox_cff import jetToolbox
    jetToolbox(process, 'ak15', 'dummySeqAK15', 'noOutput',
        PUMethod='Puppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels,
        Cut='pt > 160.0 && abs(rapidity()) < 2.4',
        runOnMC=runOnMC,
        addNsub=True, maxTau=3,
        addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels,
        bTagDiscriminators=bTagDiscriminators, subjetBTagDiscriminators=subjetBTagDiscriminators, 
        addEnergyCorrFunc=True, ecfType = "N", ecfBeta = 1.0, ecfN3 = False,
        addEnergyCorrFuncSubjets=True, ecfSubjetType = "N", ecfSubjetBeta = 1.0, ecfSubjetN3 = False,
    )

    ### EXPERIMENTAL
    ### Try to add jetCorrFactors like in jets_cff.py
    from  PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
    from  PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets
    process.jetCorrFactorsAK15 = patJetCorrFactors.clone(src='packedPatJetsAK15PFPuppiSoftDrop',
        levels = cms.vstring('L1FastJet',
            'L2Relative',
            'L3Absolute',
        'L2L3Residual'),
        payload = cms.string('AK8PFPuppi'),
        primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    )

    process.updatedJetsAK15 = updatedPatJets.clone(
        addBTagInfo=False,
        jetSource='packedPatJetsAK15PFPuppiSoftDrop',
        jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactorsAK15") ),
    )

    #print("\nFlag 2")
    #for aatt in dir(process):
    #    if "ak4" in aatt.lower():
    #        print(aatt)

    if runOnMC:
        process.ak15GenJetsNoNu.jetPtMin = 100
        process.ak15GenJetsNoNuSoftDrop.jetPtMin = 100

    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

    from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsProbs as pfParticleNetJetTagsProbs

    if runParticleNet:
        bTagDiscriminators += pfParticleNetJetTagsProbs
    if runParticleNetMD:
        ak15_flav_names = ["probQCDothers", "probXbb", "probXcc", "probXqq"]  # FIXME
        pfMassDecorrelatedParticleNetJetTagsProbs = ['pfMassDecorrelatedParticleNetJetTags:' + n for n in ak15_flav_names]  # FIXME
        bTagDiscriminators += pfMassDecorrelatedParticleNetJetTagsProbs

    print("*** Adding ParticleNet discriminators for AK15 ***")
    updateJetCollection(
        process,
        jetSource=cms.InputTag('updatedJetsAK15'),
        rParam=1.5,
        jetCorrections=('AK8PFPuppi', cms.vstring(JETCorrLevels), 'None'),
        btagDiscriminators=bTagDiscriminators,
        postfix='AK15',
    )
    print("*** Done adding ParticleNet discriminators for AK15 ***")


    # DEBUG : Print process attributes with AK15 in name
    #for aatt in dir(process):
    #    if "ak15" in aatt.lower() or "ak4" in aatt.lower():
    #        print(aatt)
    #print("\n")
    #for aatt in dir(process):
    #    if "pfParticleNet" in aatt:
    #        print(aatt)


    # configure DeepAK15
    from RecoBTag.ONNXRuntime.pfDeepBoostedJet_cff import pfDeepBoostedJetTags as _pfDeepBoostedJetTags
    # if runParticleNet:
    #     process.pfParticleNetTagInfosAK15.jet_radius = 1.5
    #     from PhysicsTools.PFNano.ak15.pfParticleNetPreprocessParamsAK15_cfi import pfParticleNetPreprocessParamsAK15
    #     process.pfParticleNetJetTagsAK15ParticleNet.preprocessParams = pfParticleNetPreprocessParamsAK15
    #     process.pfParticleNetJetTagsAK15ParticleNet.model_path = 'PhysicsTools/PFNano/data/ParticleNet/ak15/v00/ParticleNet-symbol.json'
    #     process.pfParticleNetJetTagsAK15ParticleNet.param_path = 'PhysicsTools/PFNano/data/ParticleNet/ak15/v00/ParticleNet-0000.params'

    if runParticleNetMD:
        from RecoBTag.ONNXRuntime.pfParticleNet_cff import pfMassDecorrelatedParticleNetJetTags
        process.pfParticleNetTagInfosAK15.jet_radius = 1.5                                                                                                                                                
        process.pfMassDecorrelatedParticleNetJetTagsAK15ParticleNet = pfMassDecorrelatedParticleNetJetTags.clone(
            src = 'pfParticleNetTagInfos',
            preprocess_json = 'PhysicsTools/PFNano/data/ParticleNet-MD/ak15/v00/preprocess.json',
            model_path = 'PhysicsTools/PFNano/data/ParticleNet-MD/ak15/v00/particle-net.onnx',
        )

#         process.pfParticleNetTagInfosAK15.jet_radius = 1.5
#         from PhysicsTools.PFNano.ak15.pfMassDecorrelatedParticleNetPreprocessParamsAK15_cfi import pfMassDecorrelatedParticleNetPreprocessParamsAK15
#         process.pfMassDecorrelatedParticleNetJetTagsAK15 = _pfDeepBoostedJetTags.clone(
#             src = process.pfMassDecorrelatedParticleNetJetTagsAK15.src,
#             flav_names = ak15_flav_names,
#             preprocessParams = pfMassDecorrelatedParticleNetPreprocessParamsAK15,
#             model_path = 'PhysicsTools/PFNano/data/ParticleNet-MD/ak15/ParticleNetMD.onnx',
# #             debugMode=True
#             )

    # src
    srcJets = cms.InputTag('selectedUpdatedPatJetsAK15')

    # jetID
    process.looseJetIdAK15Puppi = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
            version=cms.string('WINTER16'),
            quality=cms.string('LOOSE'),
        ),
        src=srcJets
    )

    process.tightJetIdAK15Puppi = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
            version=cms.string('SUMMER18PUPPI'),
            quality=cms.string('TIGHT'),
        ),
        src=srcJets
    )

    process.tightJetIdLepVetoAK15Puppi = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
            version=cms.string('SUMMER18PUPPI'),
            quality=cms.string('TIGHTLEPVETO'),
        ),
        src=srcJets
    )

    run2_jme_2016.toModify(process.tightJetIdAK15Puppi.filterParams, version="WINTER16")
    run2_jme_2016.toModify(process.tightJetIdLepVetoAK15Puppi.filterParams, version="WINTER16")
    run2_jme_2017.toModify(process.tightJetIdAK15Puppi.filterParams, version="WINTER17PUPPI")
    run2_jme_2017.toModify(process.tightJetIdLepVetoAK15Puppi.filterParams, version="WINTER17PUPPI")

    process.ak15WithUserData = cms.EDProducer("PATJetUserDataEmbedder",
        src=srcJets,
        userFloats=cms.PSet(),
        userInts=cms.PSet(
            tightId=cms.InputTag("tightJetIdAK15Puppi"),
            tightIdLepVeto=cms.InputTag("tightJetIdLepVetoAK15Puppi"),
        ),
    )
    run2_jme_2016.toModify(process.ak15WithUserData.userInts,
        looseId=cms.InputTag("looseJetIdAK15Puppi"),
        tightIdLepVeto=None,
    )

    # Final AK15 jets
    process.finalJetsAK15 = cms.EDFilter("PATJetRefSelector",
        src = cms.InputTag("ak15WithUserData"),
        cut = cms.string("") # pt > 170
    )


    process.ak15Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("finalJetsAK15"), # ak15WithUserData
        name=cms.string("FatJetAK15"), # AK15Puppi
        cut=cms.string(""),
        doc=cms.string("ak15 puppi jets"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            jetId           = Var("userInt('tightId')*2+4*userInt('tightIdLepVeto')", int, doc="Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"),
            area            = Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
            rawFactor       = Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            nPFConstituents = Var("numberOfDaughters()", int, doc="Number of PF candidate constituents"),
            tau1            = Var("userFloat('NjettinessAK15Puppi:tau1')", float, doc="Nsubjettiness (1 axis)", precision=10),
            tau2            = Var("userFloat('NjettinessAK15Puppi:tau2')", float, doc="Nsubjettiness (2 axis)", precision=10),
            tau3            = Var("userFloat('NjettinessAK15Puppi:tau3')", float, doc="Nsubjettiness (3 axis)", precision=10),
            #tau4            = Var("userFloat('NjettinessAK15Puppi:tau4')", float, doc="Nsubjettiness (4 axis)", precision=10),
            n2b1            = Var("?hasUserFloat('nb1AK15PuppiSoftDrop:ecfN2')?userFloat('nb1AK15PuppiSoftDrop:ecfN2'):-99999.", float, doc="N2 with beta=1 (for jets with raw pT>250 GeV)", precision=10),
            n3b1            = Var("?hasUserFloat('nb1AK15PuppiSoftDrop:ecfN3')?userFloat('nb1AK15PuppiSoftDrop:ecfN3'):-99999.", float, doc="N3 with beta=1 (for jets with raw pT>250 GeV)", precision=10), 
            msoftdrop       = Var("groomedMass()", float, doc="Corrected soft drop mass with PUPPI", precision=10),
            btagCSVV2       = Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc="pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            btagDeepB       = Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')", float, doc="DeepCSV b+bb tag discriminator", precision=10),
            btagJP          = Var("bDiscriminator('pfJetProbabilityBJetTags')", float, doc="pfJetProbabilityBJetTags b-tag discriminator (aka JP)", precision=10),
            nBHadrons       = Var("jetFlavourInfo().getbHadrons().size()", int, doc="number of b-hadrons"),
            nCHadrons       = Var("jetFlavourInfo().getcHadrons().size()", int, doc="number of c-hadrons"),
            subJetIdx1      = Var("?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].key():-1", int, doc="index of first subjet"),
            subJetIdx2      = Var("?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].key():-1", int, doc="index of second subjet"),
        )
    )
    run2_jme_2016.toModify(process.ak15Table.variables, jetId=Var("userInt('tightId')*2+userInt('looseId')", int, doc="Jet ID flags bit1 is loose, bit2 is tight"))
    process.ak15Table.variables.pt.precision = 10

    # add nominal taggers
    if runParticleNet:
        for prob in pfParticleNetJetTagsProbs:
            name = 'ParticleNet_' + prob.split(':')[1]
            setattr(process.ak15Table.variables, name, Var("bDiscriminator('%s')" % prob, float, doc=prob, precision=-1))

    # add mass-decorelated taggers
    if runParticleNetMD:
        for prob in pfMassDecorrelatedParticleNetJetTagsProbs:
            name = 'ParticleNetMD_' + prob.split(':')[1]
            name = name.replace('QCDothers', 'QCD')  # FIXME
            setattr(process.ak15Table.variables, name, Var("bDiscriminator('%s')" % prob, float, doc=prob, precision=-1))

    process.ak15SubJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("selectedPatJetsAK15PFPuppiSoftDropPacked", "SubJets"),
        cut=cms.string(""),
        name=cms.string("FatJetAK15SubJet"), # AK15PuppiSubJet
        doc=cms.string("ak15 puppi subjets"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            area=Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
            rawFactor=Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            btagDeepB=Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')", float, doc="DeepCSV b+bb tag discriminator", precision=10),
            btagCSVV2=Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            btagJP=Var("bDiscriminator('pfJetProbabilityBJetTags')", float, doc="pfJetProbabilityBJetTags b-tag discriminator (aka JP)", precision=10),
            nBHadrons=Var("jetFlavourInfo().getbHadrons().size()", int, doc="number of b-hadrons"),
            nCHadrons=Var("jetFlavourInfo().getcHadrons().size()", int, doc="number of c-hadrons"),
        )
    )
    process.ak15SubJetTable.variables.pt.precision = 10

    process.ak15Task = cms.Task(
        process.jetCorrFactorsAK15,
        process.updatedJetsAK15,
        process.tightJetIdAK15Puppi,
        process.tightJetIdLepVetoAK15Puppi,
        process.ak15WithUserData,
        process.finalJetsAK15,
        process.ak15Table,
        process.ak15SubJetTable,
    )

    if runOnMC:
        process.genJetAK15Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src=cms.InputTag("ak15GenJetsNoNu"),
            cut=cms.string("pt > 100."),
            name=cms.string("GenJetAK15"),
            doc=cms.string("AK15 GenJets made with visible genparticles"),
            singleton=cms.bool(False),  # the number of entries is variable
            extension=cms.bool(False),  # this is the main table for the genjets
            variables=cms.PSet(P4Vars,
            )
        )
        process.genJetAK15Table.variables.pt.precision = 10

        process.genSubJetAK15Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src=cms.InputTag("ak15GenJetsNoNuSoftDrop", "SubJets"),
            cut=cms.string(""),
            name=cms.string("GenSubJetAK15"),
            doc=cms.string("AK15 Gen-SubJets made with visible genparticles"),
            singleton=cms.bool(False),  # the number of entries is variable
            extension=cms.bool(False),  # this is the main table for the genjets
            variables=cms.PSet(P4Vars,
            )
        )
        process.genSubJetAK15Table.variables.pt.precision = 10

        process.genJetAK15SoftDropTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                         src=cms.InputTag("ak15GenJetsNoNuSoftDrop"),
                                                         cut=cms.string("pt > 100."),
                                                         name=cms.string("SoftDropGenJetAK15"),
                                                         doc=cms.string("AK15 GenJets made with visible genparticles and SD"),
                                                         singleton=cms.bool(False),
                                                         extension=cms.bool(False),
                                                         variables=cms.PSet(P4Vars,
                                                                        ))
        process.genJetAK15SoftDropTable.variables.pt.precision = 10

        process.ak15Task.add(process.genJetAK15Table)
        process.ak15Task.add(process.genSubJetAK15Table)
        process.ak15Task.add(process.genJetAK15SoftDropTable)

        ###### hack to avoid circular dependency ######
        process.jetMC.remove(process.patJetPartons)
        process.ak15Task.add(process.patJetPartons)
        ###############################################

    _ak15Task_2016 = process.ak15Task.copy()
    _ak15Task_2016.replace(process.tightJetIdLepVetoAK15Puppi, process.looseJetIdAK15Puppi)
    run2_jme_2016.toReplaceWith(process.ak15Task, _ak15Task_2016)

    if path is None:
        process.schedule.associate(process.ak15Task)
    else:
        getattr(process, path).associate(process.ak15Task)



#
# Functions for cmsDriver.py --customise
#

def setupPFNanoAK15_data(process):
    if hasattr(process, "patJetsPuppi"):
        process.patJetsPuppi.addGenPartonMatch = cms.bool(False)
        process.patJetsPuppi.addGenJetMatch = cms.bool(False)
        process.patJetsPuppi.JetPartonMapSource = cms.InputTag("")
        process.patJetsPuppi.JetFlavourInfoSource = cms.InputTag("")
    setupAK15(process, runOnMC=False, runParticleNet=True, runParticleNetMD=True)
    addPFCands(process, runOnMC=False, saveAll=False, addAK4=False, addAK8=True, addAK15=True)
    add_BTV(process, runOnMC=False, addAK4=True, addAK8=True, addAK15=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def setupPFNanoAK15_mc(process):
    setupAK15(process, runOnMC=True, runParticleNet=True, runParticleNetMD=True)
    addPFCands(process, runOnMC=False, saveAll=False, addAK4=False, addAK8=True, addAK15=True, saveAllGen=False) # runOnMC=False because we don't need GenCands associated to GenJets
    add_BTV(process, runOnMC=True, addAK4=True, addAK8=True, addAK15=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process
