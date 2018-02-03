import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAODJMAR.recluster_cff import *


nanoSequence += finalJetsAK8Constituents + finalJetsAK8ConstituentsTable

nanoSequenceMC += finalJetsAK8Constituents + genJetsAK8Constituents + finalJetsAK8ConstituentsTable + genJetsAK8ParticleTable



def nanoAOD_customizeData_JMAR(process):
    process = nanoAOD_customizeData(process)
    return process
    
def nanoAOD_customizeMC_JMAR(process):
    process = nanoAOD_customizeMC(process)
    return process
