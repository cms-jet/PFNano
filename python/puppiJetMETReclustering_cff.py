from PhysicsTools.PFNano.puppiJetMETReclusteringTools import puppiJetMETReclusterFromMiniAOD

def nanoAOD_puppiRecluster(process, useExistingWeights=False, runOnMC=True):
    process = puppiJetMETReclusterFromMiniAOD(process, useExistingWeights, runOnMC)
    return process

def nanoPuppiReclusterCustomize_MC(process):
    process = nanoAOD_puppiRecluster(process, runOnMC=True)
    return process

def nanoPuppiReclusterCustomize_Data(process):
    process = nanoAOD_puppiRecluster(process, runOnMC=False)
    return process
