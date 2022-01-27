import FWCore.ParameterSet.Config as cms
# from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.NanoAOD.jets_cff import jetTable, fatJetTable, subJetTable
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from PhysicsTools.PatAlgos.tools.helpers import addToProcessAndTask, getPatAlgosToolsTask


def update_jets_AK4(process):
    # Based on ``nanoAOD_addDeepInfo``
    # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
    # DeepJet flav_names as found in
    # https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/plugins/DeepFlavourONNXJetTagsProducer.cc#L86
    # and https://twiki.cern.ch/twiki/bin/view/CMS/DeepJet
    _btagDiscriminators = [
        'pfJetProbabilityBJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probbb',
        'pfDeepCSVJetTags:probudsg',       
        'pfDeepFlavourJetTags:probb',
        'pfDeepFlavourJetTags:probbb',
        'pfDeepFlavourJetTags:problepb',
        'pfDeepFlavourJetTags:probc',
        'pfDeepFlavourJetTags:probuds',
        'pfDeepFlavourJetTags:probg'
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
    process.updatedPatJetsTransientCorrectedWithDeepInfo.tagInfoSources.append(cms.InputTag("pfDeepFlavourTagInfosWithDeepInfo"))
    
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


def get_DeepJet_outputs():
    DeepJetOutputVars = cms.PSet(
        btagDeepFlavB_b=Var("bDiscriminator('pfDeepFlavourJetTags:probb')",
                            float,
                            doc="DeepJet b tag probability",
                            precision=10),
        btagDeepFlavB_bb=Var("bDiscriminator('pfDeepFlavourJetTags:probbb')",
                             float,
                             doc="DeepJet bb tag probability",
                             precision=10),
        btagDeepFlavB_lepb=Var("bDiscriminator('pfDeepFlavourJetTags:problepb')",
                               float,
                               doc="DeepJet lepb tag probability",
                               precision=10),
        btagDeepFlavC=Var("bDiscriminator('pfDeepFlavourJetTags:probc')",
                            float,
                            doc="DeepJet c tag probability",
                            precision=10),
        btagDeepFlavUDS=Var("bDiscriminator('pfDeepFlavourJetTags:probuds')",
                            float,
                            doc="DeepJet uds tag probability",
                            precision=10),
        btagDeepFlavG=Var("bDiscriminator('pfDeepFlavourJetTags:probg')",
                          float,
                          doc="DeepJet gluon tag probability",
                          precision=10)
        # discriminators are already part of jets_cff.py from NanoAOD and therefore not added here     
    )
    return DeepJetOutputVars


def add_BTV(process, runOnMC=False, onlyAK4=False, onlyAK8=False, keepInputs=['DeepCSV','DDX']):
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
        btagDeepL=Var("bDiscriminator('pfDeepCSVJetTags:probudsg')",
                      float,
                      doc="DeepCSV light btag discriminator",
                      precision=10),
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
            get_DeepCSV_vars() if ('DeepCSV' in keepInputs) else cms.PSet(),
            get_DeepJet_outputs()  # outputs are added in any case, inputs only if requested
        ))
    
    if ('DeepJet' in keepInputs):
        process.customAK4ConstituentsForDeepJetTable = cms.EDProducer("PatJetDeepJetTableProducer",
                                                                      jets = cms.InputTag("finalJets")
                                                                      )
    
    
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
            get_DDX_vars() if ('DDX' in keepInputs) else cms.PSet(),
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
        if ('DeepJet' in keepInputs):
            process.customizeJetTask.add(process.customAK4ConstituentsForDeepJetTable)

    if addAK8:
        process.customizeJetTask.add(process.customFatJetExtTable)
        process.customizeJetTask.add(process.customSubJetExtTable)
    if runOnMC and addAK8:
        process.customizeJetTask.add(process.customSubJetMCExtTable)

    return process
