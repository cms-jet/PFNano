#!/usr/bin/env python2

import os
import sys
#from input_crab_data import dataset_files
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser
from http.client import HTTPException
from multiprocessing import Process
import copy

from CRABClient.UserUtilities import config, ClientException, getUsername
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException

sys.path.append(".")
from production_tag import production_tag # Get from a text file
# Get from git tag (tbd)
#production_tag = "vTEST9" # Specify by hand
requestname_base = "pfnano"
output_site = "T2_US_MIT"
output_lfn_base = "/store/user/{username}/pfnano/{production_tag}".format(
                                                    username=getUsername(), 
                                                    production_tag=production_tag)

if __name__ == '__main__':

    def submit(config):
        print("DEBUG : In submit()")
        try:
            crabCommand('submit', config=config)
        except HTTPException as hte:
            print("Failed submitting task: {}".format(hte.headers))
            print(hte)
        except ClientException as cle:
            print("Failed submitting task: {}".format(cle))

    parser = ArgumentParser()
    parser.add_argument('-y', '--yaml', default = 'samples_datatest.yaml', help = 'File with dataset descriptions')
    args = parser.parse_args()

    with open(args.yaml) as f:
        doc = yaml.load(f) # Parse YAML file
        defaults = doc['defaults'] if 'defaults' in doc else {}

        for sample in sorted(doc["samples"].keys()):
            info = copy.deepcopy(defaults)
            info.update(doc["samples"][sample])
            # defaults.update(doc["samples"][sample])
            print("\n\n*** Sample {} ***".format(sample))

            for dataset_shortname, dataset in info['datasets'].iteritems():            
                print("\n*** Submitting {}: {}".format(dataset_shortname, dataset))

                #isMC = info.get("isMC", defaults.get("isMC", None))
                isMC = info.get("isMC", None)
                if isMC == None:
                    raise ValueError("Please specify parameter isMC")

                this_config = config()

                this_config.section_('General')
                this_config.General.transferOutputs = True
                this_config.General.transferLogs = True
                this_config.General.workArea = "crab/{}_{}/".format(requestname_base, production_tag)
                this_config.General.requestName = "{}_{}_{}_{}".format(requestname_base, production_tag, info["year"], dataset_shortname)

                this_config.section_('JobType')
                this_config.JobType.pluginName = 'Analysis'
                this_config.JobType.psetName = os.path.expandvars(info.get("pset", None))
                #this_config.JobType.maxJobRuntimeMin = 3000
                this_config.JobType.allowUndistributedCMSSW = True
                this_config.JobType.numCores = 4
                this_config.JobType.maxMemoryMB = 8000
                #this_config.JobType.outputFiles = ["_".join("nano", "mc" if isMC else "data")]
                #this_config.JobType.outputFiles = ["nanoskim.root", "hists.root"]
                #this_config.JobType.outputFiles = ['_'.join(['DijetSkim', 'mc' if isMC else 'data', production_tag])+'.root']
                #this_config.JobType.sendPythonFolder  = True
                globaltag = info.get("globaltag", None)
                this_config.JobType.pyCfgParams = [
                        'isMC={}'.format(isMC), 
                        'reportEvery=1000',
                        'tag={}'.format(production_tag),
                        'globalTag={}'.format(globaltag),
                ]

                this_config.section_('User')
                this_config.section_('Site')
                this_config.Site.storageSite = output_site

                this_config.section_('Data')
                this_config.Data.publication = True
                this_config.Data.outLFNDirBase = "{}/{}/{}".format(output_lfn_base, info["year"], sample)
                this_config.Data.outputDatasetTag = dataset_shortname
                # Outputs land at outLFNDirBase/outputDatasetTag
                this_config.Data.inputDBS = 'global'
                if "private" in info.keys():
                    if info["private"]:
                        this_config.Data.inputDBS = 'phys03'
                this_config.Data.inputDataset = dataset
                splitting_mode = info.get("splitting", "Automatic")
                if not splitting_mode in ["Automatic", "FileBased", "LumiBased"]:
                    raise ValueError("Unrecognized splitting mode: {}".format(splitting_mode))
                this_config.Data.splitting = splitting_mode

                if not isMC:
                        this_config.Data.lumiMask = info.get('lumimask', None)
                else:
                        this_config.Data.lumiMask = ''

                unitsPerJob = info.get("unitsPerJob", None)
                if unitsPerJob is not None:
                    this_config.Data.unitsPerJob = unitsPerJob

                totalUnits = info.get("totalUnits", None)
                if totalUnits is not None:
                    this_config.Data.totalUnits = totalUnits

                allowInvalid = info.get("allowInvalid", False)
                if allowInvalid:
                    this_config.Data.allowNonValidInputDataset = True
                
                print(this_config)
                p = Process(target=submit, args=(this_config,))
                p.start()
                p.join()
                #submit(this_config)
            print("*** Done with Sample {} ***\n\n".format(sample))

