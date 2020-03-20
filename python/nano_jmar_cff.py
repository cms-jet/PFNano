import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAODJMAR.recluster_cff import *
#from PhysicsTools.NanoAOD.custom_jme_cff import *


nanoSequenceAK4 = finalJetsAK4Constituents + finalJetsAK4ConstituentsTable
nanoSequenceAK8 = finalJetsAK8Constituents + finalJetsAK8ConstituentsTable

nanoSequenceGenAK4 = genJetsAK4Constituents + genJetsAK4ParticleTable
nanoSequenceGenAK8 = genJetsAK8Constituents + genJetsAK8ParticleTable

nanoSequence += nanoSequenceAK4 + nanoSequenceAK8
nanoSequenceMC += nanoSequenceGenAK4 + nanoSequenceGenAK8 + nanoSequenceAK4 + nanoSequenceAK8

def nanoAOD_customizeData_JMAR(process):
    #process = PrepJMECustomNanoAOD(process, runOnMC=True)
    process = nanoAOD_customizeData(process)
    return process

def nanoAOD_customizeMC_JMAR(process):
    #process = PrepJMECustomNanoAOD(process, runOnMC=False)
    process = nanoAOD_customizeMC(process)
    return process
