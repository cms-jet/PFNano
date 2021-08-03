import FWCore.ParameterSet.Config as cms
from PhysicsTools.PFNano.addPFCands_cff import addPFCands
from PhysicsTools.PFNano.addBTV import add_BTV
from PhysicsTools.NanoAOD.common_cff import Var


def PFnano_customizeMC(process):
    addPFCands(process, runOnMC=True, addAK4=True, addAK8=True)
    add_BTV(process, runOnMC=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_allPF(process):
    addPFCands(process, runOnMC=True, saveAll=True, addAK4=True, addAK8=True) # All PFCands, and tables for AK4 and AK8
    add_BTV(process, runOnMC=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process


def PFnano_customizeMC_AK4JetsOnly(process):
    addPFCands(process, runOnMC=True, addAK4=True)
    add_BTV(process, runOnMC=True, addAK4=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_AK8JetsOnly(process):
    addPFCands(process, runOnMC=True, addAK8=True)
    add_BTV(process, runOnMC=True, addAK4=False, addAK8=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeMC_noInputs(process):
    add_BTV(process, runOnMC=True, keepInputs=False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

#### DATA customization
def PFnano_customizeData(process):
    addPFCands(process, runOnMC=False, addAK4=True, addAK8=True)
    add_BTV(process, runOnMC=False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_allPF(process):
    addPFCands(process, runOnMC=False, saveAll=True, addAK4=True, addAK8=True)
    add_BTV(process, runOnMC=False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK4JetsOnly(process):
    addPFCands(process, runOnMC=False, addAK4=True)
    add_BTV(process, runOnMC=False, addAK4=True, addAK8=False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_AK8JetsOnly(process):
    addPFCands(process, runOnMC=False, addAK8=True)
    add_BTV(process, runOnMC=False, addAK4=False, addAK8=True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def PFnano_customizeData_noInputs(process):
    add_BTV(process, runOnMC=False, keepInputs=False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process
