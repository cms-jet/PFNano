import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAODJMAR.addPFCands_cff import addPFCands
from PhysicsTools.NanoAODJMAR.addBTV import add_BTV
from PhysicsTools.NanoAOD.common_cff import Var


def JMARnano_customizeMC(process):
    addPFCands(process, True)
    add_BTV(process, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def JMARnano_customizeMC_AK4JetsOnly(process):
    addPFCands(process, True, True, False)
    add_BTV(process, True, True, False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def JMARnano_customizeMC_AK8JetsOnly(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, False, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

#### DATA customization
def JMARnano_customizeData(process):
    addPFCands(process, False)
    add_BTV(process, False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def JMARnano_customizeData_AK4JetsOnly(process):
    addPFCands(process, False, True, False)
    add_BTV(process, False, True, False)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process

def JMARnano_customizeData_AK8JetsOnly(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, False, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process
