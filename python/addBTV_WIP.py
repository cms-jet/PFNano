import FWCore.ParameterSet.Config as cms
# from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.NanoAOD.jets_cff import jetTable, fatJetTable, subJetTable
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from PhysicsTools.PatAlgos.tools.helpers import addToProcessAndTask, getPatAlgosToolsTask


def update_jets_AK4(process):
    # Based on ``nanoAOD_addDeepInfo``
    # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
    _btagDiscriminators = [
        'pfJetProbabilityBJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probbb',
        'pfDeepCSVJetTags:probudsg',
        # start adding DeepFlavour (DeepJet) here, also part of nanoAOD_addDeepInfo
        # DeepJet flav_names as found in https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/plugins/DeepFlavourONNXJetTagsProducer.cc#L86
        # and https://twiki.cern.ch/twiki/bin/view/CMS/DeepJet
        'pfDeepFlavourJetTags:probb',
        'pfDeepFlavourJetTags:probbb',
        'pfDeepFlavourJetTags:problepb',
        'pfDeepFlavourJetTags:probuds',
        'pfDeepFlavourJetTags:probg',
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

    process.updatedPatJetsTransientCorrectedWithDeepInfo.tagInfoSources.append(cms.InputTag("pfDeepCSVTagInfosWithDeepInfo"))
    # DeepJet here, "append" sounds like one could have that line once for DeepCSV and again for DeepJet
    process.updatedPatJetsTransientCorrectedWithDeepInfo.tagInfoSources.append(cms.InputTag("pfDeepFlavourTagInfosWithDeepInfo"))
    
    
    # Try to also get constituents information
    process.finalJetsAK4Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                            src = cms.InputTag("finalJets"),
                                            cut = cms.string("")
                                            )
    
    candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"))
    process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
    # candidates here
    process.updatedPatJetsTransientCorrectedWithDeepInfo.tagInfoSources.append(candList)
    
    
    process.updatedPatJetsTransientCorrectedWithDeepInfo.addTagInfos = cms.bool(True)
    
    
    
    
    
    
    return process


def update_jets_AK8(process):
    # Based on ``nanoAOD_addDeepInfoAK8``
    # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
    # Care needs to be taken to make sure no discriminators from stock Nano are excluded -> would results in unfilled vars
    _btagDiscriminators = [
        'pfJetProbabilityBJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probbb',
        'pfDeepCSVJetTags:probudsg',
        'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
        'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
        'pfMassIndependentDeepDoubleCvBJetTags:probHcc',
        'pfMassIndependentDeepDoubleBvLV2JetTags:probHbb',
        'pfMassIndependentDeepDoubleCvLV2JetTags:probHcc',
        'pfMassIndependentDeepDoubleCvBV2JetTags:probHcc',
        ]
    from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
    _btagDiscriminators += pfParticleNetJetTagsAll
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
        # this should work but doesn't seem to enable the tag info with addTagInfos
        # btagInfos=['pfDeepDoubleXTagInfos'],
        printWarning=False)
    process.jetCorrFactorsAK8.src = "selectedUpdatedPatJetsAK8WithDeepInfo"
    process.updatedJetsAK8.jetSource = "selectedUpdatedPatJetsAK8WithDeepInfo"
    # add DeepDoubleX taginfos
    process.updatedPatJetsTransientCorrectedAK8WithDeepInfo.tagInfoSources.append(cms.InputTag("pfDeepDoubleXTagInfosAK8WithDeepInfo"))
    process.updatedPatJetsTransientCorrectedAK8WithDeepInfo.addTagInfos = cms.bool(True)
    return process


def update_jets_AK8_subjet(process):
    # Based on ``nanoAOD_addDeepInfoAK8``
    # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
    # and https://github.com/alefisico/RecoBTag-PerformanceMeasurements/blob/10_2_X_boostedCommissioning/test/runBTagAnalyzer_cfg.py
    _btagDiscriminators = [
        'pfJetProbabilityBJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probbb',
        'pfDeepCSVJetTags:probudsg',
        ]
    updateJetCollection(
        process,
        labelName='SoftDropSubjetsPF',
        jetSource=cms.InputTag("slimmedJetsAK8PFPuppiSoftDropPacked", "SubJets"),
        jetCorrections=('AK4PFPuppi',
                        ['L2Relative', 'L3Absolute'], 'None'),
        btagDiscriminators=list(_btagDiscriminators),
        explicitJTA=True,  # needed for subjet b tagging
        svClustering=False,  # needed for subjet b tagging (IMPORTANT: Needs to be set to False to disable ghost-association which does not work with slimmed jets)
        fatJets=cms.InputTag('slimmedJetsAK8'),  # needed for subjet b tagging
        rParam=0.8,  # needed for subjet b tagging
        sortByPt=False, # Don't change order (would mess with subJetIdx for FatJets)
        postfix='AK8SubjetsWithDeepInfo')

    process.subJetTable.src = 'selectedUpdatedPatJetsSoftDropSubjetsPFAK8SubjetsWithDeepInfo' 
    

    return process

def get_DDX_vars():
    # retreive 27 jet-level features used in double-b and deep double-x taggers
    # defined in arXiv:1712.07158
    DDXVars = cms.PSet(
        DDX_jetNTracks = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.jetNTracks", int, doc="number of tracks associated with the jet"),
        DDX_jetNSecondaryVertices = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.jetNSecondaryVertices", int, doc="number of SVs associated with the jet"),
        DDX_tau1_trackEtaRel_0 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_trackEtaRel_0", float, doc="1st smallest track pseudorapidity, relative to the jet axis, associated to the 1st N-subjettiness axis", precision=10),
        DDX_tau1_trackEtaRel_1 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_trackEtaRel_1", float, doc="2nd smallest track pseudorapidity, relative to the jet axis, associated to the 1st N-subjettiness axis", precision=10),
        DDX_tau1_trackEtaRel_2 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_trackEtaRel_2", float, doc="3rd smallest track pseudorapidity, relative to the jet axis, associated to the 1st N-subjettiness axis", precision=10),
        DDX_tau2_trackEtaRel_0 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau2_trackEtaRel_0", float, doc="1st smallest track pseudorapidity, relative to the jet axis, associated to the 2nd N-subjettiness axis", precision=10),
        DDX_tau2_trackEtaRel_1 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau2_trackEtaRel_1", float, doc="2nd smallest track pseudorapidity, relative to the jet axis, associated to the 2nd N-subjettiness axis", precision=10),
        DDX_tau2_trackEtaRel_3 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau2_trackEtaRel_2", float, doc="3rd smallest track pseudorapidity, relative to the jet axis, associated to the 2nd N-subjettiness axis", precision=10),
        DDX_tau1_flightDistance2dSig = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_flightDistance2dSig", float, doc="transverse distance significance between primary and secondary vertex associated to the 1st N-subjettiness axis", precision=10),
        DDX_tau2_flightDistance2dSig = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau2_flightDistance2dSig", float, doc="transverse distance significance between primary and secondary vertex associated to the 2nd N-subjettiness axis", precision=10),
        DDX_tau1_vertexDeltaR = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_vertexDeltaR", float, doc="deltaR between the 1st N-subjettiness axis and secondary vertex direction", precision=10),
        DDX_tau1_vertexEnergyRatio = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_vertexEnergyRatio", float, doc="ratio of energy at secondary vertex over total energy associated to the 1st N-subjettiness axis", precision=10),
        DDX_tau2_vertexEnergyRatio = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau2_vertexEnergyRatio", float, doc="ratio of energy at secondary vertex over total energy associated to the 2nd N-subjettiness axis", precision=10),
        DDX_tau1_vertexMass = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_vertexMass", float, doc="mass of track sum at secondary vertex associated to the 1st N-subjettiness axis", precision=10),
        DDX_tau2_vertexMass = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau2_vertexMass", float, doc="mass of track sum at secondary vertex associated to the 2nd N-subjettiness axis", precision=10),
        DDX_trackSip2dSigAboveBottom_0 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.trackSip2dSigAboveBottom_0", float, doc="track 2D signed impact parameter significance of 1st track lifting mass above bottom", precision=10),
        DDX_trackSip2dSigAboveBottom_1 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.trackSip2dSigAboveBottom_1", float, doc="track 2D signed impact parameter significance of 2nd track lifting mass above bottom", precision=10),
        DDX_trackSip2dSigAboveCharm = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.trackSip2dSigAboveCharm", float, doc="track 2D signed impact parameter significance of 1st track lifting mass above charm", precision=10),
        DDX_trackSip3dSig_0 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.trackSip3dSig_0", float, doc="1st largest track 3D signed impact parameter significance", precision=10),
        DDX_tau1_trackSip3dSig_0 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_trackSip3dSig_0", float, doc="1st largest track 3D signed impact parameter significance associated to the 1st N-subjettiness axis", precision=10),
        DDX_tau1_trackSip3dSig_1 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau1_trackSip3dSig_1", float, doc="2nd largest track 3D signed impact parameter significance associated to the 1st N-subjettiness axis", precision=10),
        DDX_trackSip3dSig_1 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.trackSip3dSig_1", float, doc="2nd largest track 3D signed impact parameter significance", precision=10),
        DDX_tau2_trackSip3dSig_0 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau2_trackSip3dSig_0", float, doc="1st largest track 3D signed impact parameter significance associated to the 2nd N-subjettiness axis", precision=10),
        DDX_tau2_trackSip3dSig_1 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.tau2_trackSip3dSig_1", float, doc="2nd largest track 3D signed impact parameter significance associated to the 2nd N-subjettiness axis", precision=10),
        DDX_trackSip3dSig_2 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.trackSip3dSig_2", float, doc="3rd largest track 3D signed impact parameter significance", precision=10),
        DDX_trackSip3dSig_3 = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.trackSip3dSig_3", float, doc="4th largest track 3D signed impact parameter significance", precision=10),
        DDX_z_ratio = Var("tagInfo(\'pfDeepDoubleX\').features().tag_info_features.z_ratio", float, doc="z = deltaR(SV0,SV1)*pT(SV1)/m(SV0,SV1), defined in Eq. 7 of arXiv:1712.07158", precision=10)
    )
    return DDXVars

def get_DeepCSV_vars():
    DeepCSVVars = cms.PSet(
        # Tagger inputs also include jet pt and eta
        # Track based
        DeepCSV_trackPtRel_0 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRel\')[0]:-999", float, doc="track transverse momentum, relative to the jet axis", precision=10),
        DeepCSV_trackPtRel_1 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRel\')[1]:-999", float, doc="track transverse momentum, relative to the jet axis", precision=10),
        DeepCSV_trackPtRel_2 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRel\')[2]:-999", float, doc="track transverse momentum, relative to the jet axis", precision=10),
        DeepCSV_trackPtRel_3 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRel\')[3]:-999", float, doc="track transverse momentum, relative to the jet axis", precision=10),
        DeepCSV_trackPtRel_4 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRel\')[4]:-999", float, doc="track transverse momentum, relative to the jet axis", precision=10),
        DeepCSV_trackPtRel_5 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRel\')[5]:-999", float, doc="track transverse momentum, relative to the jet axis", precision=10),
        DeepCSV_trackJetDistVal_0 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackJetDistVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackJetDistVal\')[0]:-999", float, doc="minimum track approach distance to jet axis", precision=10),
        DeepCSV_trackJetDistVal_1 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackJetDistVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackJetDistVal\')[1]:-999", float, doc="minimum track approach distance to jet axis", precision=10),
        DeepCSV_trackJetDistVal_2 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackJetDistVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackJetDistVal\')[2]:-999", float, doc="minimum track approach distance to jet axis", precision=10),
        DeepCSV_trackJetDistVal_3 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackJetDistVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackJetDistVal\')[3]:-999", float, doc="minimum track approach distance to jet axis", precision=10),
        DeepCSV_trackJetDistVal_4 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackJetDistVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackJetDistVal\')[4]:-999", float, doc="minimum track approach distance to jet axis", precision=10),
        DeepCSV_trackJetDistVal_5 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackJetDistVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackJetDistVal\')[5]:-999", float, doc="minimum track approach distance to jet axis", precision=10),
        DeepCSV_trackDeltaR_0 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDeltaR\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDeltaR\')[0]:-999", float, doc="track pseudoangular distance from the jet axis", precision=10),
        DeepCSV_trackDeltaR_1 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDeltaR\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDeltaR\')[1]:-999", float, doc="track pseudoangular distance from the jet axis", precision=10),
        DeepCSV_trackDeltaR_2 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDeltaR\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDeltaR\')[2]:-999", float, doc="track pseudoangular distance from the jet axis", precision=10),
        DeepCSV_trackDeltaR_3 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDeltaR\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDeltaR\')[3]:-999", float, doc="track pseudoangular distance from the jet axis", precision=10),
        DeepCSV_trackDeltaR_4 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDeltaR\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDeltaR\')[4]:-999", float, doc="track pseudoangular distance from the jet axis", precision=10),
        DeepCSV_trackDeltaR_5 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDeltaR\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDeltaR\')[5]:-999", float, doc="track pseudoangular distance from the jet axis", precision=10),
        DeepCSV_trackPtRatio_0 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRatio\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRatio\')[0]:-999", float, doc="track transverse momentum, relative to the jet axis, normalized to its energy", precision=10),
        DeepCSV_trackPtRatio_1 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRatio\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRatio\')[1]:-999", float, doc="track transverse momentum, relative to the jet axis, normalized to its energy", precision=10),
        DeepCSV_trackPtRatio_2 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRatio\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRatio\')[2]:-999", float, doc="track transverse momentum, relative to the jet axis, normalized to its energy", precision=10),
        DeepCSV_trackPtRatio_3 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRatio\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRatio\')[3]:-999", float, doc="track transverse momentum, relative to the jet axis, normalized to its energy", precision=10),
        DeepCSV_trackPtRatio_4 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRatio\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRatio\')[4]:-999", float, doc="track transverse momentum, relative to the jet axis, normalized to its energy", precision=10),
        DeepCSV_trackPtRatio_5 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRatio\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRatio\')[5]:-999", float, doc="track transverse momentum, relative to the jet axis, normalized to its energy", precision=10),
        DeepCSV_trackSip3dSig_0 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip3dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip3dSig\')[0]:-999", float, doc="track 3D signed impact parameter significance", precision=10),
        DeepCSV_trackSip3dSig_1 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip3dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip3dSig\')[1]:-999", float, doc="track 3D signed impact parameter significance", precision=10),
        DeepCSV_trackSip3dSig_2 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip3dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip3dSig\')[2]:-999", float, doc="track 3D signed impact parameter significance", precision=10),
        DeepCSV_trackSip3dSig_3 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip3dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip3dSig\')[3]:-999", float, doc="track 3D signed impact parameter significance", precision=10),
        DeepCSV_trackSip3dSig_4 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip3dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip3dSig\')[4]:-999", float, doc="track 3D signed impact parameter significance", precision=10),
        DeepCSV_trackSip3dSig_5 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip3dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip3dSig\')[5]:-999", float, doc="track 3D signed impact parameter significance", precision=10),
        DeepCSV_trackSip2dSig_0 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip2dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip2dSig\')[0]:-999", float, doc="track 2D signed impact parameter significance", precision=10),
        DeepCSV_trackSip2dSig_1 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip2dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip2dSig\')[1]:-999", float, doc="track 2D signed impact parameter significance", precision=10),
        DeepCSV_trackSip2dSig_2 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip2dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip2dSig\')[2]:-999", float, doc="track 2D signed impact parameter significance", precision=10),
        DeepCSV_trackSip2dSig_3 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip2dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip2dSig\')[3]:-999", float, doc="track 2D signed impact parameter significance", precision=10),
        DeepCSV_trackSip2dSig_4 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip2dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip2dSig\')[4]:-999", float, doc="track 2D signed impact parameter significance", precision=10),
        DeepCSV_trackSip2dSig_5 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackSip2dSig\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackSip2dSig\')[5]:-999", float, doc="track 2D signed impact parameter significance", precision=10),
        DeepCSV_trackDecayLenVal_0 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDecayLenVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDecayLenVal\')[0]:-999", float, doc="track decay length", precision=10),
        DeepCSV_trackDecayLenVal_1 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDecayLenVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDecayLenVal\')[1]:-999", float, doc="track decay length", precision=10),
        DeepCSV_trackDecayLenVal_2 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDecayLenVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDecayLenVal\')[2]:-999", float, doc="track decay length", precision=10),
        DeepCSV_trackDecayLenVal_3 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDecayLenVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDecayLenVal\')[3]:-999", float, doc="track decay length", precision=10),
        DeepCSV_trackDecayLenVal_4 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDecayLenVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDecayLenVal\')[4]:-999", float, doc="track decay length", precision=10),
        DeepCSV_trackDecayLenVal_5 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackDecayLenVal\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackDecayLenVal\')[5]:-999", float, doc="track decay length", precision=10),
        DeepCSV_trackEtaRel_0 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackEtaRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackEtaRel\')[0]:-999", float, doc="track pseudorapidity, relative to the jet axis", precision=10),
        DeepCSV_trackEtaRel_1 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackEtaRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackEtaRel\')[1]:-999", float, doc="track pseudorapidity, relative to the jet axis", precision=10),
        DeepCSV_trackEtaRel_2 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackEtaRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackEtaRel\')[2]:-999", float, doc="track pseudorapidity, relative to the jet axis", precision=10),
        DeepCSV_trackEtaRel_3 = Var("?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackEtaRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackEtaRel\')[3]:-999", float, doc="track pseudorapidity, relative to the jet axis", precision=10),
        # Jet based
        DeepCSV_trackJetPt = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackJetPt\', -999)", float, doc="track-based jet transverse momentum", precision=10),
        DeepCSV_vertexCategory = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'vertexCategory\', -999)", float, doc="category of secondary vertex (Reco, Pseudo, No)", precision=10),
        DeepCSV_jetNSecondaryVertices = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'jetNSecondaryVertices\', -999)", int, doc="number of reconstructed possible secondary vertices in jet"),
        DeepCSV_jetNSelectedTracks = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'jetNSelectedTracks\', -999)", int, doc="selected tracks in the jet"), 
        DeepCSV_jetNTracksEtaRel = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'jetNTracksEtaRel\', -999)", int, doc="number of tracks for which etaRel is computed"), 
        DeepCSV_trackSumJetEtRatio = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSumJetEtRatio\', -999)", float, doc="ratio of track sum transverse energy over jet energy", precision=10),
        DeepCSV_trackSumJetDeltaR = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSumJetDeltaR\', -999)", float, doc="pseudoangular distance between jet axis and track fourvector sum", precision=10),
        DeepCSV_trackSip2dValAboveCharm = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSip2dValAboveCharm\', -999)", float, doc="track 2D signed impact parameter of first track lifting mass above charm", precision=10),
        DeepCSV_trackSip2dSigAboveCharm = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSip2dSigAboveCharm\', -999)", float, doc="track 2D signed impact parameter significance of first track lifting mass above charm", precision=10),
        DeepCSV_trackSip3dValAboveCharm = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSip3dValAboveCharm\', -999)", float, doc="track 3D signed impact parameter of first track lifting mass above charm", precision=10),
        DeepCSV_trackSip3dSigAboveCharm = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSip3dSigAboveCharm\', -999)", float, doc="track 3D signed impact parameter significance of first track lifting mass above charm", precision=10),
        DeepCSV_vertexMass = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'vertexMass\', -999)", float, doc="mass of track sum at secondary vertex", precision=10),
        DeepCSV_vertexNTracks = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'vertexNTracks\', -999)", int, doc="number of tracks at secondary vertex"),
        DeepCSV_vertexEnergyRatio = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'vertexEnergyRatio\', -999)", float, doc="ratio of energy at secondary vertex over total energy", precision=10),
        DeepCSV_vertexJetDeltaR = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'vertexJetDeltaR\', -999)", float, doc="pseudoangular distance between jet axis and secondary vertex direction", precision=10),
        DeepCSV_flightDistance2dVal = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'flightDistance2dVal\', -999)", float, doc="transverse distance between primary and secondary vertex", precision=10),
        DeepCSV_flightDistance2dSig = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'flightDistance2dSig\', -999)", float, doc="transverse distance significance between primary and secondary vertex", precision=10),
        DeepCSV_flightDistance3dVal = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'flightDistance3dVal\', -999)", float, doc="distance between primary and secondary vertex", precision=10),
        DeepCSV_flightDistance3dSig = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'flightDistance3dSig\', -999)", float, doc="distance significance between primary and secondary vertex", precision=10),
    )
    return DeepCSVVars




def get_DeepJet_vars():
    # inputs for DeepJet as described in 2008.10519 and https://github.com/DL4Jets/DeepJet/blob/master/modules/datastructures/TrainData_deepFlavour.py
    # producer: https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/plugins/DeepFlavourTagInfoProducer.cc
    # looking at https://github.com/cms-sw/cmssw/blob/master/DataFormats/BTauReco/interface/DeepFlavourFeatures.h the cpf, npf anf sv features are vectors, so it should be possible to access them via (examples for every case)
    #     tagInfo('\pfDeepFlavour\').features().jet_features.pt  # global (https://github.com/cms-sw/cmssw/blob/master/DataFormats/BTauReco/interface/JetFeatures.h)
    #     tagInfo('\pfDeepFlavour\').features().tag_info_features.trackSumJetEtRatio  # global, shallow (https://github.com/cms-sw/cmssw/blob/master/DataFormats/BTauReco/interface/ShallowTagInfoFeatures.h)
    #     tagInfo('\pfDeepFlavour\').features().c_pf_features.at(index).bTagPf_trackPtRel  # index must be taken care of, because this might need an offset according to n_jets and n_candidates in the event,
    #         not sure currently (otherwise one could simply try 0 to 24, if they exist, otherwise --> default)
    #         (https://github.com/cms-sw/cmssw/blob/master/DataFormats/BTauReco/interface/ChargedCandidateFeatures.h)
    #    tagInfo('\pfDeepFlavour\').features().n_pf_features.at(index).deltaR (https://github.com/cms-sw/cmssw/blob/master/DataFormats/BTauReco/interface/NeutralCandidateFeatures.h)
    #    tagInfo('\pfDeepFlavour\').features().sv_features.at(index).deltaR (https://github.com/cms-sw/cmssw/blob/master/DataFormats/BTauReco/interface/SecondaryVertexFeatures.h)
    #    maybe something extra for npv (std::size_t ?)
    #    there is also SeedingTrackFeatures used for DeepFlavourFeatures.h, but I never saw them beeing used for DeepJet? Are these features relevant as well?
    #    another thing I am not sure about: why does jet_pt occur isolated from the global jet features in https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/plugins/DeepFlavourONNXJetTagsProducer.cc ? 
    #
    #    might have to check first if the number of requested candidates is contained in the size of the container (e.g. if a jet only has 5 cpf candidates, it makes no sense to try to read any more or fill these with something other than a default value)
    #        for that: adapt the syntax found for DeepCSV (e.g. "?tagInfo(\'pfDeepCSV\').taggingVariables.checkTag(\'trackPtRel\')?tagInfo(\'pfDeepCSV\').taggingVariables.getList(\'trackPtRel\')[5]:-999") to
    #        e.g. ?....features.c_pf_features.size()>requested_n?then_read_it:some_negative_default_value
    #
    # (just for information, for DeepCSV, there is https://github.com/cms-sw/cmssw/blob/master/DataFormats/BTauReco/interface/BaseTagInfo.h and
    # and https://github.com/cms-sw/cmssw/blob/master/DataFormats/BTauReco/interface/TaggingVariable.h)
    #
    # another idea: wouldn't it be nice to have the features also here in a jagged format? Otherwise I'd write probably 600 lines of code for every input, which is of course possible, but...
    #     or one could prepare the many lines and exec() the pre-build command there, something like that (append _0, _1, ... and so on) or write a small script from which one can copy the lines that should go here :-)
    # will probably try the idea with one candidate first to see if that works, then add more candidates later to the PSet in some way.
    DeepJetVars = cms.PSet(
        # global (probably *not all* are necessary, because that's already part of DeepCSV Jet based, see above get_DeepCSV_vars, or of the normal jet collection) --> commented out
        #DeepJet_jet_pt = Var("tagInfo(\'pfDeepFlavour\').features().jet_features.pt", float, doc="pt of the jet", precision=10),
        #DeepJet_jet_eta = Var("tagInfo(\'pfDeepFlavour\').features().jet_features.eta", float, doc="eta of the jet", precision=10),
        DeepJet_nCpfcand = Var("tagInfo('\pfDeepFlavour\').features().c_pf_features.size()", int, doc="number of charged PFCands in jet", precision=10),
        DeepJet_nNpfcand = Var("tagInfo('\pfDeepFlavour\').features().n_pf_features.size()", int, doc="number of neutral PFCands in jet", precision=10),
        DeepJet_nsv = Var("tagInfo('\pfDeepFlavour\').features().sv_features.size()", int, doc="number of reconstructed possible secondary vertices in jet", precision=10),
        DeepJet_npv = Var("tagInfo(\'pfDeepFlavour\').features().npv", int, doc="number of primary vertices", precision=10),
        
        # global, shallow (probably not necessary, because that's already part of DeepCSV Jet based, see above get_DeepCSV_vars)
        ### this is one example to read shallow info via DeepFlavour
        # already in DeepCSV DeepJet_trackSumJetEtRatio = Var("tagInfo(\'pfDeepFlavour\').features().tag_info_features.trackSumJetEtRatio", float, doc="ratio of track sum transverse energy over jet energy", precision=10),
        ### and here are all others, but read via DeepCSV tag info (copied from above - there's the -999, where I don't know if I should keep it like that)
        # not part of DeepJet currently DeepJet_trackJetPt = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackJetPt\', -999)", float, doc="track-based jet transverse momentum", precision=10),
        # already in DeepCSV DeepJet_vertexCategory = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'vertexCategory\', -999)", float, doc="category of secondary vertex (Reco, Pseudo, No)", precision=10),
        # see global DeepJet_jetNSecondaryVertices = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'jetNSecondaryVertices\', -999)", int, doc="number of reconstructed possible secondary vertices in jet"),
        # already in DeepCSV DeepJet_jetNSelectedTracks = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'jetNSelectedTracks\', -999)", int, doc="selected tracks in the jet"), 
        # already in DeepCSV DeepJet_jetNTracksEtaRel = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'jetNTracksEtaRel\', -999)", int, doc="number of tracks for which etaRel is computed"), 
        # alraedy done above as DeepFlavour example DeepJet_trackSumJetEtRatio = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSumJetEtRatio\', -999)", float, doc="ratio of track sum transverse energy over jet energy", precision=10),
        # already in DeepCSV DeepJet_trackSumJetDeltaR = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSumJetDeltaR\', -999)", float, doc="pseudoangular distance between jet axis and track fourvector sum", precision=10),
        # already in DeepCSV DeepJet_trackSip2dValAboveCharm = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSip2dValAboveCharm\', -999)", float, doc="track 2D signed impact parameter of first track lifting mass above charm", precision=10),
        # already in DeepCSV DeepJet_trackSip2dSigAboveCharm = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSip2dSigAboveCharm\', -999)", float, doc="track 2D signed impact parameter significance of first track lifting mass above charm", precision=10),
        # already in DeepCSV DeepJet_trackSip3dValAboveCharm = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSip3dValAboveCharm\', -999)", float, doc="track 3D signed impact parameter of first track lifting mass above charm", precision=10),
        # already in DeepCSV DeepJet_trackSip3dSigAboveCharm = Var("tagInfo(\'pfDeepCSV\').taggingVariables.get(\'trackSip3dSigAboveCharm\', -999)", float, doc="track 3D signed impact parameter significance of first track lifting mass above charm", precision=10),
        
        # CPF
        #    what's inside https://algomez.web.cern.ch/algomez/testWeb/PFnano_content_v02.html#JetPFCands ? Wouldn't this result in duplicates?
        DeepJet_Cpfcan_BtagPf_trackEtaRel_0     = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackEtaRel:-1",     float, doc="track pseudorapidity, relative to the jet axis", precision=10),
        DeepJet_Cpfcan_BtagPf_trackPtRel_0      = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).bTagPf_trackPtRel:-1",      float, doc="track transverse momentum, relative to the jet axis", precision=10),
        DeepJet_Cpfcan_BtagPf_trackPPar_0       = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackPPar:-1",       float, doc="dot product of the jet and track momentum", precision=10),
        DeepJet_Cpfcan_BtagPf_trackPParRatio_0  = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackPParRatio:-1",  float, doc="dot product of the jet and track momentum divided by the magnitude of the jet momentum", precision=10),
        DeepJet_Cpfcan_BtagPf_trackDeltaR_0     = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackDeltaR:-1",     float, doc="track pseudoangular distance from the jet axis", precision=10),
        DeepJet_Cpfcan_BtagPf_trackSip2dVal_0   = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackSip2dVal:-1",   float, doc="track 2D signed impact parameter", precision=10),
        DeepJet_Cpfcan_BtagPf_trackSip2dSig_0   = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackSip2dSig:-1",   float, doc="track 2D signed impact parameter significance", precision=10),
        DeepJet_Cpfcan_BtagPf_trackSip3dVal_0   = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackSip3dVal:-1",   float, doc="track 3D signed impact parameter", precision=10),
        DeepJet_Cpfcan_BtagPf_trackSip3dSig_0   = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackSip3dSig:-1",   float, doc="track 3D signed impact parameter significance", precision=10),
        DeepJet_Cpfcan_BtagPf_trackJetDistVal_0 = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).btagPf_trackJetDistVal:-1", float, doc="minimum track approach distance to jet axis", precision=10),
        DeepJet_Cpfcan_ptrel_0                  = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).ptrel:-1",                  float, doc="fraction of the jet momentum carried by the track", precision=10),
        DeepJet_Cpfcan_drminsv_0                = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).drminsv:-1",                float, doc="track pseudoangular distance from the closest secondary vertex", precision=10),
        DeepJet_Cpfcan_vtx_ass_0                = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).vtx_ass:-1",                int,   doc="integer flag that indicates whether the track was used in the primary vertex fit", precision=10),
        DeepJet_Cpfcan_puppiw_0                 = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).puppiw:-1",                 float, doc="charged candidate PUPPI weight", precision=10),
        DeepJet_Cpfcan_chi2_0                   = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).chi2:-1",                   float, doc="chi2 of the charged track fit", precision=10),
        DeepJet_Cpfcan_quality_0                = Var("?tagInfo('\pfDeepFlavour\').features().c_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().c_pf_features.at(0).quality:-1",                int,   doc="integer flag which indicates the quality of the fitted track, based on number of detector hits used for the reconstruction as well as the overall chi2 of the charged track fit", precision=10),
        
        # NPF
        DeepJet_Npfcan_ptrel_0   = Var("?tagInfo('\pfDeepFlavour\').features().n_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().n_pf_features.at(0).ptrel:-1",   float, doc="fraction of the jet momentum carried by the neutral candidate", precision=10),
        DeepJet_Npfcan_deltaR_0  = Var("?tagInfo('\pfDeepFlavour\').features().n_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().n_pf_features.at(0).deltaR:-1",  float, doc="pseudoangular distance between the neutral candidate and the jet axis", precision=10),
        DeepJet_Npfcan_isGamma_0 = Var("?tagInfo('\pfDeepFlavour\').features().n_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().n_pf_features.at(0).isGamma:-1", int,   doc="integer flag indicating whether the neutral candidate is a photon", precision=10),
        DeepJet_Npfcan_HadFrac_0 = Var("?tagInfo('\pfDeepFlavour\').features().n_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().n_pf_features.at(0).hadFrac:-1", float, doc="fraction of the neutral candidate energy deposited in the hadronic calorimeter", precision=10),
        DeepJet_Npfcan_drminsv_0 = Var("?tagInfo('\pfDeepFlavour\').features().n_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().n_pf_features.at(0).drminsv:-1", float, doc="pseudoangular distance between the neutral candidate and the closest secondary vertex", precision=10),
        DeepJet_Npfcan_puppiw_0  = Var("?tagInfo('\pfDeepFlavour\').features().n_pf_features.size()>0?tagInfo('\pfDeepFlavour\').features().n_pf_features.at(0).puppiw:-1",  float, doc="neutral candidate PUPPI weight", precision=10),
        
        # SV
        #    same question, there is already JetSVs https://algomez.web.cern.ch/algomez/testWeb/PFnano_content_v02.html#JetSVs
        DeepJet_sv_pt_0           = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).pt:-1",           float, doc="pt of the secondary vertex", precision=10),
        DeepJet_sv_deltaR_0       = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).deltaR:-1",       float, doc="pseudoangular distance between jet axis and secondary vertex direction", precision=10),
        DeepJet_sv_mass_0         = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).mass:-1",         float, doc="secondary vertex mass", precision=10),
        DeepJet_sv_ntracks_0      = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).ntracks:-1",      int,   doc="number of tracks in the secondary vertex", precision=10),
        DeepJet_sv_chi2_0         = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).chi2:-1",         float, doc="chi2 of the secondary vertex fit", precision=10),
        DeepJet_sv_normchi2_0     = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).normchi2:-1",     float, doc="reduced chi2 of the secondary vertex fit", precision=10),
        DeepJet_sv_dxy_0          = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).dxyy:-1",         float, doc="secondary vertex 2D impact parameter value", precision=10),
        DeepJet_sv_dxysig_0       = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).dxysig:-1",       float, doc="secondary vertex 2D impact parameter significance", precision=10),
        DeepJet_sv_d3d_0          = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).d3d:-1",          float, doc="secondary vertex 3D impact parameter value", precision=10),
        DeepJet_sv_d3dsig_0       = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).d3dsig:-1",       float, doc="secondary vertex 3D impact parameter significance", precision=10),
        DeepJet_sv_costhetasvpv_0 = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).costhetasvpv:-1", float, doc="cosine of the angle between the secondary vertex flight direction and the direction of the secondary vertex momentum", precision=10),
        DeepJet_sv_enratio_0      = Var("?tagInfo('\pfDeepFlavour\').features().sv_features.size()>0?tagInfo('\pfDeepFlavour\').features().sv_features.at(0).enratio:-1",      float, doc="ratio of the secondary vertex energy ratio to the jet energy", precision=10),
        
    )
    return DeepJetVars



def add_BTV(process, runOnMC=False, onlyAK4=False, onlyAK8=False, keepInputs=True):
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
        btagDeepB_b=Var("bDiscriminator('pfDeepCSVJetTags:probb')",
                        float,
                        doc="DeepCSV b tag discriminator",
                        precision=10),
        btagDeepB_bb=Var("bDiscriminator('pfDeepCSVJetTags:probbb')",
                         float,
                         doc="DeepCSV bb tag discriminator",
                         precision=10),
        btagDeepL=Var("bDiscriminator('pfDeepCSVJetTags:probudsg')",  # actually, where is this inside https://algomez.web.cern.ch/algomez/testWeb/PFnano_content_v02.html#Jet ? It appears that instead, bTagDeepC is there
                      float,
                      doc="DeepCSV light btag discriminator",
                      precision=10),
        # in case DeepJet should also be part of CommonVars (which I don't know), then DeepJet starts here
        # ToDo: ask if that's the right place and what needs to go here and what not
        #btagDeepFlavB_b=Var("bDiscriminator('pfDeepFlavourJetTags:probb')",
        #                float,
        #                doc="DeepJet b tag discriminator",
        #                precision=10),
        #btagDeepFlavB_bb=Var("bDiscriminator('pfDeepFlavourJetTags:probbb')",
        #                 float,
        #                 doc="DeepJet bb tag discriminator",
        #                 precision=10),
        #btagDeepFlavB_lepb=Var("bDiscriminator('pfDeepFlavourJetTags:problepb')",
        #                 float,
        #                 doc="DeepJet lepb tag discriminator",
        #                 precision=10),
        #btagDeepFlavUDS=Var("bDiscriminator('pfDeepFlavourJetTags:probuds')",
        #              float,
        #              doc="DeepJet uds discriminator",
        #              precision=10),
        #btagDeepFlavG=Var("bDiscriminator('pfDeepFlavourJetTags:probg')",
        #              float,
        #              doc="DeepJet gluon discriminator",
        #              precision=10),
        #btagDeepFlavL=Var("bDiscriminator('pfDeepFlavourJetTags:probuds')+bDiscriminator('pfDeepFlavourJetTags:probg')",
        #              float,
        #              doc="DeepJet light discriminator",
        #              precision=10),
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
            get_DeepCSV_vars() if keepInputs else cms.PSet(),
            # similar to https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/custom_jme_cff.py#L185 (?) - they have an additional pt cut, which might not be necessary here (sets defaults to -1 for pt<15),
            #     in the same source: ParticleNetAK4 equivalent placed below DeepJet
            # ToDo: ask if that's the right place and what needs to go here and what not, some are already there in PFNano it seems
            cms.PSet(
                btagDeepFlavB_b    = Var("?pt>15?bDiscriminator('pfDeepFlavourJetTags:probb'):-1",float,doc="DeepJet b tag raw score",precision=10),
                btagDeepFlavB_bb   = Var("?pt>15?bDiscriminator('pfDeepFlavourJetTags:probbb'):-1",float,doc="DeepJet bb tag raw score",precision=10),
                btagDeepFlavB_lepb = Var("?pt>15?bDiscriminator('pfDeepFlavourJetTags:problepb'):-1",float,doc="DeepJet lepb tag raw score",precision=10),
                #btagDeepFlavB      = Var("?pt>15?bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb'):-1",float,doc="DeepJet b+bb+lepb tag discriminator",precision=10),
                #btagDeepFlavC      = Var("?pt>15?bDiscriminator('pfDeepFlavourJetTags:probc'):-1",float,doc="DeepJet c tag raw score",precision=10),
                btagDeepFlavG      = Var("?pt>15?bDiscriminator('pfDeepFlavourJetTags:probg'):-1",float,doc="DeepJet gluon tag raw score",precision=10),
                btagDeepFlavUDS    = Var("?pt>15?bDiscriminator('pfDeepFlavourJetTags:probuds'):-1",float,doc="DeepJet uds tag raw score",precision=10),
                btagDeepFlavL      = Var("?pt>15?(bDiscriminator('pfDeepFlavourJetTags:probuds')+bDiscriminator('pfDeepFlavourJetTags:probg')):-1",float,doc="DeepJet uds+g tag discriminator",precision=10),
                # for DeepCSV it looks like probc is not included in CommonVars, probably because it can easily be calculated from 1-(sum of others)
                # should this be done for DeepJet as well (like, not include probc here?)
                # looked at  https://algomez.web.cern.ch/algomez/testWeb/PFnano_content_v02.html#Jet to find out what's currently already there, commented out accordingly
                # also don't know if derived discriminators (fractions of raw scores) are necessary or only raw scores (+ basic sums to get B and L) - again found in custom_jme_cff, copied from there, lacks B versus X
                #btagDeepFlavCvL = Var("?(pt>15)&&(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probuds')+bDiscriminator('pfDeepFlavourJetTags:probg'))>0?bDiscriminator('pfDeepFlavourJetTags:probc')/(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probuds')+bDiscriminator('pfDeepFlavourJetTags:probg')):-1",float,doc="DeepJet c vs uds+g discriminator",precision=10),
                #btagDeepFlavCvB = Var("?(pt>15)&&(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb'))>0?bDiscriminator('pfDeepFlavourJetTags:probc')/(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb')):-1",float,doc="DeepJet c vs b+bb+lepb discriminator",precision=10),
                #btagDeepFlavQG  = Var("?(pt>15)&&(bDiscriminator('pfDeepFlavourJetTags:probg')+bDiscriminator('pfDeepFlavourJetTags:probuds'))>0?bDiscriminator('pfDeepFlavourJetTags:probg')/(bDiscriminator('pfDeepFlavourJetTags:probg')+bDiscriminator('pfDeepFlavourJetTags:probuds')):-1",float,doc="DeepJet g vs uds discriminator",precision=10),
            ),
            get_DeepJet_vars() if keepInputs else cms.PSet(),
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
            cms.PSet(
                btagDDBvLV2 = Var("bDiscriminator('pfMassIndependentDeepDoubleBvLV2JetTags:probHbb')",float,doc="DeepDoubleX V2 discriminator for H(Z)->bb vs QCD",precision=10),
                btagDDCvLV2 = Var("bDiscriminator('pfMassIndependentDeepDoubleCvLV2JetTags:probHcc')",float,doc="DeepDoubleX V2 discriminator for H(Z)->cc vs QCD",precision=10),
                btagDDCvBV2 = Var("bDiscriminator('pfMassIndependentDeepDoubleCvBV2JetTags:probHcc')",float,doc="DeepDoubleX V2 discriminator for H(Z)->cc vs H(Z)->bb",precision=10),
            ),
            get_DDX_vars() if keepInputs else cms.PSet(),
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
             btagDeepC = Var("bDiscriminator('pfDeepCSVJetTags:probc')",
                        float,
                        doc="DeepCSV charm btag discriminator",
                        precision=10),

    ))

    process.customSubJetMCExtTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = subJetTable.src,
    cut = subJetTable.cut,
        name = subJetTable.name,
        doc=subJetTable.doc,
    singleton = cms.bool(False),
       extension = cms.bool(True),
        variables = cms.PSet(
            subGenJetAK8Idx = Var("?genJetFwdRef().backRef().isNonnull()?genJetFwdRef().backRef().key():-1",
        int,
        doc="index of matched gen Sub jet"),
       )
    )

    if addAK4:
        process.customizeJetTask.add(process.customJetExtTable)
    
        '''
        process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituents")
        process.customConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                            src = candInput,
                                                            cut = cms.string(""), #we should not filter after pruning
                                                            name = cms.string("PFCands"),
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
        '''
        process.customAK4ConstituentsTableDJ = cms.EDProducer("PatJetConstituentTableProducerDeepJet",
                                                            #candidates = cms.InputTag("packedPFCandidates"),
                                                            candidates = candInput,
                                                            jets = cms.InputTag("finalJets"),
                                                            jet_radius = cms.double(0.4),
                                                            name = cms.string("JetPFCands"),
                                                            idx_name = cms.string("pFCandsIdx"),
                                                            nameSV = cms.string("JetSVs"),
                                                            idx_nameSV = cms.string("sVIdx"),
                                                            nameDeepJet = cms.string("JetDeepJet"),
                                                            idx_nameDeepJet = cms.string("djIdx"),
                                                            )

    
    
    if addAK8:
        process.customizeJetTask.add(process.customFatJetExtTable)
        process.customizeJetTask.add(process.customSubJetExtTable)
        #process.customizedPFCandsTask.add(process.customConstituentsExtTable)
        process.customizedPFCandsTask.add(process.customAK4ConstituentsTableDJ)
    if runOnMC and addAK8:
        process.customizeJetTask.add(process.customSubJetMCExtTable)

    return process
