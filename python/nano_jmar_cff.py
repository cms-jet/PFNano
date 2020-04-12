import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAODJMAR.addPFCands_cff import addPFCands
from PhysicsTools.NanoAOD.common_cff import Var


def JMARnano_customizeMC(process):
    addPFCands(process, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def JMARnano_customizeMC_AK4JetsOnly(process):
    addPFCands(process, True, True, False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def JMARnano_customizeMC_AK8JetsOnly(process):
    addPFCands(process, True, False, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

#### DATA customization
def JMARnano_customizeData(process):
    addPFCands(process, False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def JMARnano_customizeData_AK4JetsOnly(process):
    addPFCands(process, False, True, False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def JMARnano_customizeData_AK8JetsOnly(process):
    addPFCands(process, False, False, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process
