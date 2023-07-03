import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.jetTools import setupPackedPuppi
from PhysicsTools.PatAlgos.tools.jetTools import supportedJetAlgos, addJetCollection, updateJetCollection
from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask, addToProcessAndTask

###################
#
# Setup puppi
#
###################
def setupPuppi(process, useExistingWeights):
  puppiLabel = setupPackedPuppi(process)
  task = getPatAlgosToolsTask(process)
  getattr(process,puppiLabel).useExistingWeights = useExistingWeights
  return puppiLabel

def puppiJetMETReclusterFromMiniAOD(process, useExistingWeights, runOnMC):

  task = getPatAlgosToolsTask(process)

  pfLabel = "packedPFCandidates"
  pvLabel = "offlineSlimmedPrimaryVertices"
  svLabel = "slimmedSecondaryVertices"
  muLabel = "slimmedMuons"
  elLabel = "slimmedElectrons"
  gpLabel = "prunedGenParticles"

  genJetsCollection = "slimmedGenJets"
  genJetsCollectionAK8 = "slimmedGenJetsAK8"
  genSubJetsForAK8Collection = "slimmedGenJetsAK8SoftDropSubJets"

  JETCorrLevels = ["L2Relative", "L3Absolute"]

  #########################
  #
  # Setup puppi weights
  # Two instances of PuppiProducer:
  # 1) packedpuppi (for jet reclustering)
  # 2) puppiNoLep (for MET reclustering)
  #
  ########################
  puppiLabel = setupPuppi(process, useExistingWeights)
  setattr(process, "puppiNoLep", getattr(process, puppiLabel).clone(
      puppiNoLep = True
    )
  )
  task.add(getattr(process, "puppiNoLep"))

  #########################
  #
  # AK4 Puppi jets
  #
  ########################
  #
  # Recluster jets
  #
  process.load("RecoJets.JetProducers.ak4PFJets_cfi")
  task.add(process.ak4PFJetsPuppi)
  process.ak4PFJetsPuppi.src = pfLabel
  process.ak4PFJetsPuppi.srcWeights = puppiLabel

  from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import j2tParametersVX
  process.ak4PFJetsPuppiTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    j2tParametersVX,
    jets = cms.InputTag("ak4PFJetsPuppi")
  )
  task.add(process.ak4PFJetsPuppiTracksAssociatorAtVertex)
  process.patJetPuppiCharge = cms.EDProducer("JetChargeProducer",
     src = cms.InputTag("ak4PFJetsPuppiTracksAssociatorAtVertex"),
     var = cms.string("Pt"),
     exp = cms.double(1.0)
  )

  #
  # PATify jets
  #
  addJetCollection(
    process,
    postfix            = "",
    labelName          = "Puppi",
    jetSource          = cms.InputTag("ak4PFJetsPuppi"),
    algo               = "ak",
    rParam             = 0.4,
    pfCandidates       = cms.InputTag(pfLabel),
    pvSource           = cms.InputTag(pvLabel),
    svSource           = cms.InputTag(svLabel),
    muSource           = cms.InputTag(muLabel),
    elSource           = cms.InputTag(elLabel),
    genJetCollection   = cms.InputTag(genJetsCollection),
    genParticles       = cms.InputTag(gpLabel),
    jetCorrections     = ('AK4PFPuppi', cms.vstring(["L2Relative", "L3Absolute"]), ''),
    getJetMCFlavour    = runOnMC
  )
  process.patJetsPuppi.jetChargeSource = cms.InputTag("patJetPuppiCharge")
  process.selectedPatJetsPuppi.cut = cms.string("pt > 10")
  if hasattr(process,"patJetFlavourAssociationPuppi"):
    process.patJetFlavourAssociationPuppi.weights = cms.InputTag(puppiLabel)

  #=============================================
  #
  # Update the selectedPatJet collection.
  # This is where we setup
  # -  JEC
  # -  b-tagging discriminators
  #
  #=============================================
  # update slimmedJetsPuppi to include deep taggers
  from PhysicsTools.PatAlgos.slimming.slimmedJets_cfi import slimmedJets
  addToProcessAndTask('slimmedJetsPuppiNoDeepTags', slimmedJets.clone(
      src = "selectedPatJetsPuppi",
      packedPFCandidates = "packedPFCandidates",
      dropDaughters = "0",
      rekeyDaughters = "0",
    ),
    process, task
  )

  from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll as pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll
  from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll as pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll
  from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4CHSCentralJetTagsAll as pfParticleNetFromMiniAODAK4CHSCentralJetTagsAll
  from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4CHSForwardJetTagsAll as pfParticleNetFromMiniAODAK4CHSForwardJetTagsAll
  bTagDeepCSV  = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfDeepCSVJetTags:probc','pfDeepCSVJetTags:probudsg']

  _btagDiscriminatorsAK4Puppi = cms.PSet(
    names=cms.vstring(
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
      'pfDeepFlavourJetTags:probc',
      'pfDeepFlavourJetTags:probuds',
      'pfDeepFlavourJetTags:probg'
    ) + pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll
    + pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll + bTagDeepCSV
  )

  updateJetCollection(
    process,
    jetSource = cms.InputTag("slimmedJetsPuppiNoDeepTags"),
    # updateJetCollection defaults to MiniAOD inputs but
    # here it is made explicit (as in training or MINIAOD redoing)
    pfCandidates = cms.InputTag(pfLabel),
    pvSource = cms.InputTag(pvLabel),
    svSource = cms.InputTag(svLabel),
    muSource = cms.InputTag(muLabel),
    elSource = cms.InputTag(elLabel),
    jetCorrections = ("AK4PFPuppi", cms.vstring(["L2Relative", "L3Absolute"]), "None"),
    btagDiscriminators = _btagDiscriminatorsAK4Puppi.names.value(),
    postfix = 'SlimmedDeepFlavour',
    printWarning = False
  )

  #
  #
  #
  addToProcessAndTask("slimmedJetsPuppi", process.selectedUpdatedPatJetsSlimmedDeepFlavour.clone(), process, task)
  del process.selectedUpdatedPatJetsSlimmedDeepFlavour

  ########################
  #
  # AK8 Puppi jets
  #
  ########################

  #
  # Recluster jets and do soft-drop grooming
  #
  process.load("RecoJets.JetProducers.ak8PFJets_cfi")
  task.add(process.ak8PFJetsPuppi)
  task.add(process.ak8PFJetsPuppiSoftDrop)

  # AK8 jet constituents for softdrop
  process.ak8PFJetsPuppi.src = pfLabel
  process.ak8PFJetsPuppi.srcWeights = puppiLabel

  # AK8 jet constituents for softdrop
  from CommonTools.RecoAlgos.miniAODJetConstituentSelector_cfi import miniAODJetConstituentSelector
  addToProcessAndTask("ak8PFJetsPuppiConstituents", miniAODJetConstituentSelector.clone(
      src = cms.InputTag("ak8PFJetsPuppi"),
      cut = cms.string("pt > 100.0 && abs(rapidity()) < 2.4")
    ),
    process, task
  )

  # Soft-drop grooming
  process.ak8PFJetsPuppiSoftDrop.src = "ak8PFJetsPuppiConstituents:constituents"
  process.ak8PFJetsPuppiSoftDrop.srcWeights = puppiLabel

  # Soft-drop mass
  process.load("RecoJets.JetProducers.ak8PFJetsPuppi_groomingValueMaps_cfi")
  task.add(process.ak8PFJetsPuppiSoftDropMass)
  process.ak8PFJetsPuppiSoftDropMass.src = "ak8PFJetsPuppi"
  process.ak8PFJetsPuppiSoftDropMass.matched = "ak8PFJetsPuppiSoftDrop"
  #=============================================
  #
  # PATify
  #
  #=============================================
  #
  # AK8 jets
  #
  addJetCollection(
    process,
    labelName          = "AK8Puppi",
    jetSource          = cms.InputTag("ak8PFJetsPuppi"),
    algo               = "ak",
    rParam             = 0.8,
    pfCandidates       = cms.InputTag(pfLabel),
    pvSource           = cms.InputTag(pvLabel),
    svSource           = cms.InputTag(svLabel),
    muSource           = cms.InputTag(muLabel),
    elSource           = cms.InputTag(elLabel),
    genJetCollection   = cms.InputTag(genJetsCollectionAK8),
    genParticles       = cms.InputTag(gpLabel),
    jetCorrections     = ("AK8PFPuppi", cms.vstring(["L2Relative", "L3Absolute"]), "None"),
    getJetMCFlavour    = runOnMC,
    btagDiscriminators = ([
        # "pfCombinedSecondaryVertexV2BJetTags",
        # "pfCombinedInclusiveSecondaryVertexV2BJetTags",
        # "pfCombinedMVAV2BJetTags",
        "pfDeepCSVJetTags:probb",
        "pfDeepCSVJetTags:probc",
        "pfDeepCSVJetTags:probudsg",
        "pfDeepCSVJetTags:probbb",
        "pfBoostedDoubleSecondaryVertexAK8BJetTags"
      ]
    ),
  )
  if hasattr(process,"patJetFlavourAssociationAK8Puppi"):
    process.patJetFlavourAssociationAK8Puppi.weights = cms.InputTag(puppiLabel)

  process.patJetsAK8Puppi.userData.userFloats.src = [] # start with empty list of user floats
  process.patJetsAK8Puppi.userData.userFloats.src += ["ak8PFJetsPuppiSoftDropMass"]
  process.patJetsAK8Puppi.addTagInfos = cms.bool(False)

  process.selectedPatJetsAK8Puppi.cut = cms.string("pt > 100")
  process.selectedPatJetsAK8Puppi.cutLoose = cms.string("pt > 30")
  process.selectedPatJetsAK8Puppi.nLoose = cms.uint32(3)

  #
  # Add AK8 Njetiness
  #
  from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
  addToProcessAndTask("NjettinessAK8Puppi", Njettiness.clone(
      src = "ak8PFJetsPuppi",
      srcWeights = puppiLabel
    ),
    process, task
  )
  process.patJetsAK8Puppi.userData.userFloats.src += [
    "NjettinessAK8Puppi:tau1",
    "NjettinessAK8Puppi:tau2",
    "NjettinessAK8Puppi:tau3",
    "NjettinessAK8Puppi:tau4"
  ]

  #
  # AK8 soft-drop jets
  #
  addJetCollection(
    process,
    labelName = "AK8PFPuppiSoftDrop",
    jetSource = cms.InputTag("ak8PFJetsPuppiSoftDrop"),
    btagDiscriminators = ["None"],
    pfCandidates = cms.InputTag(pfLabel),
    pvSource = cms.InputTag(pvLabel),
    svSource = cms.InputTag(svLabel),
    muSource = cms.InputTag(muLabel),
    elSource = cms.InputTag(elLabel),
    genJetCollection = cms.InputTag(genJetsCollectionAK8),
    genParticles = cms.InputTag(gpLabel),
    jetCorrections = ("AK8PFPuppi", cms.vstring(["L2Relative", "L3Absolute"]), "None"),
    getJetMCFlavour = False # jet flavor disabled regardless if running on MC or data
  )

  #
  # Soft-drop subjets
  #
  addJetCollection(
    process,
    labelName = "AK8PFPuppiSoftDropSubjets",
    jetSource = cms.InputTag("ak8PFJetsPuppiSoftDrop", "SubJets"),
    algo = "ak",  # needed for subjet flavor clustering
    rParam = 0.8, # needed for subjet flavor clustering
    explicitJTA = True,  # needed for subjet b tagging
    svClustering = True, # needed for subjet b tagging
    pfCandidates = cms.InputTag(pfLabel),
    pvSource = cms.InputTag(pvLabel),
    svSource = cms.InputTag(svLabel),
    muSource = cms.InputTag(muLabel),
    elSource = cms.InputTag(elLabel),
    genJetCollection = cms.InputTag(genSubJetsForAK8Collection),
    genParticles = cms.InputTag(gpLabel),
    fatJets = cms.InputTag("ak8PFJetsPuppi"),               # needed for subjet flavor clustering
    groomedFatJets = cms.InputTag("ak8PFJetsPuppiSoftDrop"), # needed for subjet flavor clustering
    jetCorrections = ("AK4PFPuppi", cms.vstring(["L2Relative", "L3Absolute"]), "None"),
    getJetMCFlavour = runOnMC,
    btagDiscriminators = [
      "pfDeepCSVJetTags:probb",
      "pfDeepCSVJetTags:probbb",
      # "pfCombinedInclusiveSecondaryVertexV2BJetTags",
      # "pfCombinedMVAV2BJetTags"
    ],
  )
  if hasattr(process,"patJetFlavourAssociationAK8PFPuppiSoftDropSubjets"):
    process.patJetFlavourAssociationAK8PFPuppiSoftDropSubjets.weights = cms.InputTag(puppiLabel)

  #=============================================
  #
  #
  #
  #=============================================
  #
  # add groomed ECFs and N-subjettiness to soft dropped pat::Jets for fat jets and subjets
  #
  process.load('RecoJets.JetProducers.ECF_cff')

  addToProcessAndTask('nb1AK8PuppiSoftDrop', process.ecfNbeta1.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop"),
      srcWeights = puppiLabel,
      cuts = cms.vstring('', '', 'pt > 250')
    ),
    process, task
  )
  process.patJetsAK8PFPuppiSoftDrop.userData.userFloats.src += [
    'nb1AK8PuppiSoftDrop:ecfN2',
    'nb1AK8PuppiSoftDrop:ecfN3',
  ]

  addToProcessAndTask('nb2AK8PuppiSoftDrop', process.ecfNbeta2.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop"),
      srcWeights = puppiLabel,
      cuts = cms.vstring('', '', 'pt > 250')
    ),
    process, task
  )
  process.patJetsAK8PFPuppiSoftDrop.userData.userFloats.src += [
    'nb2AK8PuppiSoftDrop:ecfN2',
    'nb2AK8PuppiSoftDrop:ecfN3',
  ]

  #
  # add groomed ECFs and N-subjettiness to soft drop subjets
  #
  addToProcessAndTask("nb1AK8PuppiSoftDropSubjets", process.ecfNbeta1.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop", "SubJets"),
      srcWeights = puppiLabel,
    ),
    process, task
  )

  process.patJetsAK8PFPuppiSoftDropSubjets.userData.userFloats.src += [
    'nb1AK8PuppiSoftDropSubjets:ecfN2',
    'nb1AK8PuppiSoftDropSubjets:ecfN3'
  ]

  addToProcessAndTask("nb2AK8PuppiSoftDropSubjets", process.ecfNbeta2.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop", "SubJets"),
      srcWeights = puppiLabel,
    ),
    process, task
  )

  process.patJetsAK8PFPuppiSoftDropSubjets.userData.userFloats.src += [
    'nb2AK8PuppiSoftDropSubjets:ecfN2',
    'nb2AK8PuppiSoftDropSubjets:ecfN3'
  ]

  addToProcessAndTask("NjettinessAK8Subjets", Njettiness.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop", "SubJets"),
      srcWeights = puppiLabel
    ),
    process, task
  )
  process.patJetsAK8PFPuppiSoftDropSubjets.userData.userFloats.src += [
    "NjettinessAK8Subjets:tau1",
    "NjettinessAK8Subjets:tau2",
    "NjettinessAK8Subjets:tau3",
    "NjettinessAK8Subjets:tau4",
  ]

  addToProcessAndTask("slimmedJetsAK8PFPuppiSoftDropSubjets", cms.EDProducer("PATJetSlimmer",
      src = cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropSubjets"),
      packedPFCandidates = cms.InputTag("packedPFCandidates"),
      dropJetVars = cms.string("1"),
      dropDaughters = cms.string("0"),
      rekeyDaughters = cms.string("0"),
      dropTrackRefs = cms.string("1"),
      dropSpecific = cms.string("1"),
      dropTagInfos = cms.string("1"),
      modifyJets = cms.bool(True),
      mixedDaughters = cms.bool(False),
      modifierConfig = cms.PSet( modifications = cms.VPSet() )
    ),
    process, task
  )

  ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
  ## FIKRI: Take subjets and put it in the grommed jet
  addToProcessAndTask("slimmedJetsAK8PFPuppiSoftDropPacked", cms.EDProducer("BoostedJetMerger",
      jetSrc    = cms.InputTag("selectedPatJetsAK8PFPuppiSoftDrop"),
      subjetSrc = cms.InputTag("slimmedJetsAK8PFPuppiSoftDropSubjets")
    ),
    process, task
  )

  addToProcessAndTask("packedPatJetsAK8", cms.EDProducer("JetSubstructurePacker",
      jetSrc = cms.InputTag("selectedPatJetsAK8Puppi"),
      distMax = cms.double(0.8),
      algoTags = cms.VInputTag(
        cms.InputTag("slimmedJetsAK8PFPuppiSoftDropPacked")
      ),
      algoLabels = cms.vstring(
        'SoftDropPuppi'
      ),
      fixDaughters = cms.bool(False),
      packedPFCandidates = cms.InputTag("packedPFCandidates"),
    ),
    process, task
  )

  #=============================================
  #
  # Update the selectedPatJet collection.
  # This is where we setup
  # -  JEC
  # -  b-tagging discriminators
  #
  #=============================================
  from PhysicsTools.PatAlgos.slimming.slimmedJets_cfi import slimmedJetsAK8
  addToProcessAndTask("slimmedJetsAK8NoDeepTags", slimmedJetsAK8.clone(rekeyDaughters = "0"), process, task)
  # Reconfigure the slimmedAK8 jet information to keep
  process.slimmedJetsAK8NoDeepTags.dropDaughters = cms.string("pt < 170")
  process.slimmedJetsAK8NoDeepTags.dropSpecific = cms.string("pt < 170")
  process.slimmedJetsAK8NoDeepTags.dropTagInfos = cms.string("pt < 170")

  from RecoBTag.ONNXRuntime.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsAll as pfDeepBoostedJetTagsAll
  from RecoBTag.ONNXRuntime.pfHiggsInteractionNet_cff import _pfHiggsInteractionNetTagsProbs as pfHiggsInteractionNetTagsProbs
  from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
  from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetMassRegressionOutputs as pfParticleNetMassRegressionOutputs
  from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK8_cff import _pfParticleNetFromMiniAODAK8JetTagsAll as pfParticleNetFromMiniAODAK8JetTagsAll

  _btagDiscriminatorsAK8 = cms.PSet(names = cms.vstring(
      'pfMassIndependentDeepDoubleBvLV2JetTags:probHbb',
      'pfMassIndependentDeepDoubleCvLV2JetTags:probHcc',
      'pfMassIndependentDeepDoubleCvBV2JetTags:probHcc',
    )  + pfDeepBoostedJetTagsAll +  pfParticleNetJetTagsAll + pfParticleNetMassRegressionOutputs
       + pfHiggsInteractionNetTagsProbs + pfParticleNetFromMiniAODAK8JetTagsAll
  )

  updateJetCollection(
    process,
    jetSource = cms.InputTag("slimmedJetsAK8NoDeepTags"),
    # updateJetCollection defaults to MiniAOD inputs but
    # here it is made explicit (as in training or MINIAOD redoing)
    pfCandidates = cms.InputTag(pfLabel),
    pvSource = cms.InputTag(pvLabel),
    svSource = cms.InputTag(svLabel),
    muSource = cms.InputTag(muLabel),
    elSource = cms.InputTag(elLabel),
    rParam = 0.8,
    jetCorrections = ('AK8PFPuppi', cms.vstring(["L2Relative", "L3Absolute"]), 'None'),
    btagDiscriminators = _btagDiscriminatorsAK8.names.value(),
    postfix = "SlimmedAK8DeepTags",
    printWarning = False
  )

  addToProcessAndTask("slimmedJetsAK8", process.selectedUpdatedPatJetsSlimmedAK8DeepTags.clone(), process, task)
  del process.selectedUpdatedPatJetsSlimmedAK8DeepTags

  #
  #
  #
  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
  runMetCorAndUncFromMiniAOD(process,
    isData=not(runOnMC),
    jetCollUnskimmed="slimmedJetsPuppi",
    metType="Puppi",
    postfix="Puppi",
    jetFlavor="AK4PFPuppi",
    recoMetFromPFCs=True
  )

  #
  # Modify JECs when processing real Data
  # Disable any MC-only features.
  #
  if not(runOnMC):
    from PhysicsTools.PatAlgos.tools.coreTools import runOnData
    runOnData( process, names=["Jets"], outputModules = [])

  return process
