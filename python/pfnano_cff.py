import FWCore.ParameterSet.Config as cms
from PhysicsTools.PFNano.addPFCands_cff import addPFCands
from PhysicsTools.PFNano.addBTV import add_BTV
from PhysicsTools.NanoAOD.common_cff import Var

# keepInputs can take DeepCSV, DeepJet and DDX (any combination, or use empty placeholder list if no inputs are required)
# from 12_6_X onwards, the default changes to:
#  - for AK4: keep DeepCSV (only jet-based vars), and DeepJet (leading 3 nPF/cPF candidates)
#  - for AK8: keep DDX
def PFnano_customizeMC(process):
    addPFCands(process, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeMC_allPF(process):
    addPFCands(process, True, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeMC_AK4JetsOnly(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, True, keepInputs=['DeepCSV','DeepJet'])
    return process

def PFnano_customizeMC_AK8JetsOnly(process):
    addPFCands(process, True, False, False, True)
    add_BTV(process, True, False, True, keepInputs=['DDX'])
    return process

def PFnano_customizeMC_noPF(process):
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeMC_noInputs(process):
    add_BTV(process, True, keepInputs=[])
    return process


#### DATA customization
def PFnano_customizeData(process):
    addPFCands(process, False)
    add_BTV(process, False, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeData_allPF(process):
    addPFCands(process, False, True)
    add_BTV(process, False, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeData_AK4JetsOnly(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, True, keepInputs=['DeepCSV','DeepJet'])
    return process

def PFnano_customizeData_AK8JetsOnly(process):
    addPFCands(process, False, False, False, True)
    add_BTV(process, False, False, True, keepInputs=['DDX'])
    return process

def PFnano_customizeData_noPF(process):
    add_BTV(process, False, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeData_noInputs(process):
    add_BTV(process, False, keepInputs=[])
    return process
