import FWCore.ParameterSet.Config as cms
from PhysicsTools.PFNano.addPFCands_cff import addPFCands
from PhysicsTools.PFNano.addBTV import add_BTV
#from PhysicsTools.PFNano.addBTV_WIP import add_BTV # modified script that should also be able to include DeepJet is here (currently not working)
from PhysicsTools.NanoAOD.common_cff import Var


def PFnano_customizeMC(process):
    addPFCands(process, True)
    add_BTV(process, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_add_DeepJet(process):
    addPFCands(process, True)
    add_BTV(process, True, add_DeepJet=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_add_DeepJet_noclip(process):
    addPFCands(process, True)
    add_BTV(process, True, add_DeepJet=True, add_DeepJet_noclip=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_allPF(process):
    addPFCands(process, True, True)
    add_BTV(process, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_allPF_add_DeepJet(process):
    addPFCands(process, True, True)
    add_BTV(process, True, add_DeepJet=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_allPF_add_DeepJet_noclip(process):
    addPFCands(process, True, True)
    add_BTV(process, True, add_DeepJet=True, add_DeepJet_noclip=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_AK4JetsOnly(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_AK4JetsOnly_add_DeepJet(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, True, add_DeepJet=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_AK4JetsOnly_add_DeepJet_noclip(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, True, add_DeepJet=True, add_DeepJet_noclip=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_AK8JetsOnly(process):
    addPFCands(process, True, False, False, True)
    add_BTV(process, True, False, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_noInputs(process):
    add_BTV(process, True, keepInputs=False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

#### DATA customization
def PFnano_customizeData(process):
    addPFCands(process, False)
    add_BTV(process, False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_add_DeepJet(process):
    addPFCands(process, False)
    add_BTV(process, False, add_DeepJet=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_add_DeepJet_noclip(process):
    addPFCands(process, False)
    add_BTV(process, False, add_DeepJet=True, add_DeepJet_noclip=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_allPF(process):
    addPFCands(process, False, True)
    add_BTV(process, False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_allPF_add_DeepJet(process):
    addPFCands(process, False, True)
    add_BTV(process, False, add_DeepJet=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_allPF_add_DeepJet_noclip(process):
    addPFCands(process, False, True)
    add_BTV(process, False, add_DeepJet=True, add_DeepJet_noclip=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK4JetsOnly(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK4JetsOnly_add_DeepJet(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, True, add_DeepJet=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK4JetsOnly_add_DeepJet_noclip(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, True, add_DeepJet=True, add_DeepJet_noclip=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK8JetsOnly(process):
    addPFCands(process, False, False, False, True)
    add_BTV(process, False, False, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_noInputs(process):
    add_BTV(process, False, keepInputs=False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process
