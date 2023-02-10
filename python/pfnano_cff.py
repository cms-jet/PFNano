import FWCore.ParameterSet.Config as cms
from PhysicsTools.PFNano.addPFCands_cff import addPFCands
from PhysicsTools.PFNano.addBTV import add_BTV
from PhysicsTools.NanoAOD.common_cff import Var

# keepInputs can take DeepCSV, DeepJet and DDX (any combination, or use empty placeholder list if no inputs are required)
def PFnano_customizeMC(process):
    addPFCands(process, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DDX'])
    return process

def PFnano_customizeMC_add_DeepJet(process):
    addPFCands(process, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeMC_add_DeepJet_and_Truth(process):
    addPFCands(process, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'], storeAK4Truth="yes")
    return process

def PFnano_customizeMC_allPF(process):
    addPFCands(process, True, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DDX'])
    return process

def PFnano_customizeMC_allPF_add_DeepJet(process):
    addPFCands(process, True, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeMC_allPF_add_DeepJet_and_Truth(process):
    addPFCands(process, True, True)
    add_BTV(process, True, keepInputs=['DeepCSV','DeepJet','DDX'], storeAK4Truth="yes")
    return process

def PFnano_customizeMC_AK4JetsOnly(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, True, keepInputs=['DeepCSV'])
    return process

def PFnano_customizeMC_AK4JetsOnly_add_DeepJet(process):
    addPFCands(process, True, False, True)
    add_BTV(process, True, True, keepInputs=['DeepCSV','DeepJet'])
    return process

def PFnano_customizeMC_AK8JetsOnly(process):
    addPFCands(process, True, False, False, True)
    add_BTV(process, True, False, True, keepInputs=['DDX'])
    return process

def PFnano_customizeMC_noInputs(process):
    add_BTV(process, True, keepInputs=[])
    return process


#### DATA customization
def PFnano_customizeData(process):
    addPFCands(process, False)
    add_BTV(process, False, keepInputs=['DeepCSV','DDX'])
    return process

def PFnano_customizeData_add_DeepJet(process):
    addPFCands(process, False)
    add_BTV(process, False, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeData_allPF(process):
    addPFCands(process, False, True)
    add_BTV(process, False, keepInputs=['DeepCSV','DDX'])
    return process

def PFnano_customizeData_allPF_add_DeepJet(process):
    addPFCands(process, False, True)
    add_BTV(process, False, keepInputs=['DeepCSV','DeepJet','DDX'])
    return process

def PFnano_customizeData_AK4JetsOnly(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, True, keepInputs=['DeepCSV'])
    return process

def PFnano_customizeData_AK4JetsOnly_add_DeepJet(process):
    addPFCands(process, False, False, True)
    add_BTV(process, False, True, keepInputs=['DeepCSV','DeepJet'])
    return process

def PFnano_customizeData_AK8JetsOnly(process):
    addPFCands(process, False, False, False, True)
    add_BTV(process, False, False, True, keepInputs=['DDX'])
    return process

def PFnano_customizeData_noInputs(process):
    add_BTV(process, False, keepInputs=[])
    return process
