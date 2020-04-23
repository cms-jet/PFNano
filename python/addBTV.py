import FWCore.ParameterSet.Config as cms
# from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.NanoAOD.jets_cff import jetTable, fatJetTable, subJetTable
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection


bTagInfos = [
    'pfImpactParameterTagInfos', 'pfSecondaryVertexTagInfos',
    'pfInclusiveSecondaryVertexFinderTagInfos', 'pfSecondaryVertexNegativeTagInfos',
    'pfInclusiveSecondaryVertexFinderNegativeTagInfos', 'softPFMuonsTagInfos',
    'softPFElectronsTagInfos', 'pfInclusiveSecondaryVertexFinderCvsLTagInfos',
    'pfInclusiveSecondaryVertexFinderNegativeCvsLTagInfos', 'pfDeepFlavourTagInfos'
]


def update_jets_AK4(process):
    # Based on ``nanoAOD_addDeepInfo``
    # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
    _btagDiscriminators = [
        'pfJetProbabilityBJetTags'
    ]
    updateJetCollection(
        process,
        jetSource=cms.InputTag('slimmedJets'),
        jetCorrections=('AK4PFchs',
                        cms.vstring(
                            ['L1FastJet', 'L2Relative', 'L3Absolute',
                             'L2L3Residual']), 'None'),
        btagDiscriminators=_btagDiscriminators,
        postfix='WithDeepInfo',
    )
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.jetCorrFactorsNano.src = "selectedUpdatedPatJetsWithDeepInfo"
    process.updatedJets.jetSource = "selectedUpdatedPatJetsWithDeepInfo"
    return process


def update_jets_AK8(process):
    # Based on ``nanoAOD_addDeepInfoAK8``
    # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
    _btagDiscriminators = [
        'pfJetProbabilityBJetTags'
        ]
    updateJetCollection(
        process,
        jetSource=cms.InputTag('slimmedJetsAK8'),
        pvSource=cms.InputTag('offlineSlimmedPrimaryVertices'),
        svSource=cms.InputTag('slimmedSecondaryVertices'),
        rParam=0.8,
        jetCorrections=('AK8PFPuppi',
                        cms.vstring([
                            'L1FastJet', 'L2Relative', 'L3Absolute',
                            'L2L3Residual'
                        ]), 'None'),
        btagDiscriminators=_btagDiscriminators,
        postfix='AK8WithDeepInfo',
        printWarning=False)
    process.jetCorrFactorsAK8.src = "selectedUpdatedPatJetsAK8WithDeepInfo"
    process.updatedJetsAK8.jetSource = "selectedUpdatedPatJetsAK8WithDeepInfo"
    return process


def update_jets_AK8_subjet(process):
    # Based on ``nanoAOD_addDeepInfoAK8``
    # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
    # and https://github.com/alefisico/RecoBTag-PerformanceMeasurements/blob/10_2_X_boostedCommissioning/test/runBTagAnalyzer_cfg.py
    _btagDiscriminators = [
        'pfJetProbabilityBJetTags'
        ]
    updateJetCollection(
        process,
        labelName='SoftDropSubjetsPF',
        jetSource=cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropSubjets"),
        jetCorrections=('AK4PFPuppi',
                        ['L2Relative', 'L3Absolute'], 'None'),
        pfCandidates=cms.InputTag('puppi'),
        pvSource=cms.InputTag('offlineSlimmedPrimaryVertices'),
        svSource=cms.InputTag('slimmedSecondaryVertices'),
        muSource=cms.InputTag('slimmedMuons'),
        elSource=cms.InputTag('slimmedElectrons'),
        btagInfos=bTagInfos,
        btagDiscriminators=list(_btagDiscriminators),
        explicitJTA=True,  # needed for subjet b tagging
        svClustering=False,  # needed for subjet b tagging (IMPORTANT: Needs to be set to False to disable ghost-association which does not work with slimmed jets)
        fatJets=cms.InputTag('slimmedJetsAK8'),  # needed for subjet b tagging
        rParam=0.8,  # needed for subjet b tagging
        #algo='SoftDropPuppi',  # has to be defined but is not used since svClustering=False
        postfix='AK8SubjetsWithDeepInfo')

    # process.jetCorrFactorsAK8.src = "selectedUpdatedPatJetsAK8SubjetsWithDeepInfo"
    # process.updatedJetsAK8.jetSource = "selectedUpdatedPatJetsAK8SubjetsWithDeepInfo"
    process.jetCorrFactorsAK8.src = "selectedUpdatedPatJetsAK8WithDeepInfo"
    process.updatedJetsAK8.jetSource = "selectedUpdatedPatJetsAK8WithDeepInfo"
    return process


def add_BTV(process, runOnMC=False, onlyAK4=False, onlyAK8=False):
    addAK4 = not onlyAK8
    addAK8 = not onlyAK4

    if addAK4:
        process = update_jets_AK4(process)
    if addAK8:
        process = update_jets_AK8(process)
        process = update_jets_AK8_subjet(process)

    process.customizeJetTask = cms.Task()
    process.schedule.associate(process.customizeJetTask)

    CommonVars = cms.PSet(
        Proba=Var("bDiscriminator('pfJetProbabilityBJetTags')",
                  float,
                  doc="Jet Probability (Usage:BTV)",
                  precision=10),
        nBHadrons=Var("jetFlavourInfo().getbHadrons().size()",
                      int,
                      doc="number of b-hadrons"),
        nCHadrons=Var("jetFlavourInfo().getcHadrons().size()",
                      int,
                      doc="number of c-hadrons"),
    )

    # AK4
    process.customJetExtTable = cms.EDProducer(
        "SimpleCandidateFlatTableProducer",
        src=jetTable.src,
        cut=jetTable.cut,
        name=jetTable.name,
        doc=jetTable.doc,
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(True),  # this is the extension table for Jets
        variables=cms.PSet(
            CommonVars,
        ))

    # AK8
    process.customFatJetExtTable = cms.EDProducer(
        "SimpleCandidateFlatTableProducer",
        src=fatJetTable.src,
        cut=fatJetTable.cut,
        name=fatJetTable.name,
        doc=fatJetTable.doc,
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(True),  # this is the extension table for FatJets
        variables=cms.PSet(
            CommonVars,
            # Proba=Var("bDiscriminator('pfJetProbabilityBJetTags')",
            #           float,
            #           doc="Jet Probability (Usage:BTV)",
            #           precision=10),
            # nBHadrons=Var("jetFlavourInfo().getbHadrons().size()",
            #               int,
            #               doc="number of b-hadrons"),
            # nCHadrons=Var("jetFlavourInfo().getcHadrons().size()",
            #               int,
            #               doc="number of c-hadrons"),
        ))

    # Subjets
    process.customSubJetExtTable = cms.EDProducer(
        "SimpleCandidateFlatTableProducer",
        src=subJetTable.src,
        cut=subJetTable.cut,
        name=subJetTable.name,
        doc=subJetTable.doc,
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(True),  # this is the extension table for FatJets
        variables=cms.PSet(
            CommonVars,
            # Proba=Var("bDiscriminator('pfJetProbabilityBJetTags')",
            #           float,
            #           doc="Jet Probability (Usage:BTV)",
            #           precision=10),
            # nBHadrons=Var("jetFlavourInfo().getbHadrons().size()",
            #               int,
            #               doc="number of b-hadrons"),
            # nCHadrons=Var("jetFlavourInfo().getcHadrons().size()",
            #               int,
            #               doc="number of c-hadrons"),
        ))

    if addAK4:
        process.customizeJetTask.add(process.customJetExtTable)
    if addAK8:
        process.customizeJetTask.add(process.customFatJetExtTable)
        process.customizeJetTask.add(process.customSubJetExtTable)

    return process
