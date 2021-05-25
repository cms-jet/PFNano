#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
from __future__ import print_function
import os
import sys
import re
import argparse

def make_list(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('To see usage run: python submit_all.py -h')

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("-c", "--config", dest="cfg", required=True,
        help=("CMSSW job configuration file"))
    parser.add_argument("-d", "--dir", dest="dir", required=True,
        help=("Crab working dir, e.g. 'crab_dir'"))
    parser.add_argument("-f", "--datasets", dest="datasets", required=True,
        help=("File listing datasets to run over"))
    parser.add_argument("-o", "--output", dest="out", 
        default="/store/user/${USER}/PFNano",
        help=("T2/T2 storage path"))    
    parser.add_argument("-s", "--storageSite", dest="storageSite", required=True,
        help=('Storage site, example `T3_US_FNALLPC`'))
    parser.add_argument("-l", "--lumiMask", dest="lumiMask",
        help=("Lumi Mask JSON file"))
    parser.add_argument("--ext", '--extension', dest="extension",
        help=("Extension string to include in publication name"))

    parser.add_argument("--test-only", dest="test_only", action='store_true',
        help=("Submit only a few files to test process. Disable publication"))

    # For debugging
    parser.add_argument("--remote", dest="remote", action='store_true',
        help=("Run with ignore locality (add whitelist manually)")) 
    parser.add_argument("--whitelist", dest="whitelist", default="",
        help=("Whitelist to use with --remote, passed as 'T2_DE_RWTH,T2_DE_DESY'")),

    # For Germans
    parser.add_argument("--dcms", dest="dcms", action='store_true',
        help=("Running with dcms proxy"))

    options = parser.parse_args()

    # Check choices and emit warnings
    if options.cfg == None or options.dir == None or options.datasets == None or options.storageSite == None:
        parser.error(usage)

    if options.out == parser.get_default('out'):
        print("Output directory not specified, files will be stored in {}".format(options.out))
        if raw_input("Are you sure? (y/n)") != "y":
            exit()
    
    if options.lumiMask is None:
        if raw_input("No lumi mask specified. (OK for MC) Continue? (y/n)") != "y":
            exit()

    if options.test_only is False:
        if raw_input("`--test` is set to False. About to run a full production. Continue? (y/n)") != "y":
            exit()

    if options.extension is None:
        if raw_input('`--extension` is not specified. "PFNano" will be appended by default. Continue? (y/n)') != "y":
            exit()
        options.extension = 'PFNanoAOD'
    return options


def main():
    options = getOptions()
    from CRABAPI.RawCommand import crabCommand
    from WMCore.Configuration import Configuration
    config = Configuration()
    from httplib import HTTPException


    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.section_("General")
    config.General.workArea = options.dir
    config.General.transferLogs = True

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.cfg
    config.JobType.maxMemoryMB = 5000 # Default is 2500 : Max I have used is 13000
    # minutes tied to not automatic splitting
    # config.JobType.maxJobRuntimeMin = 2750 #Default is 1315; 2750 minutes guaranteed to be available; Max I have used is 9000 
    config.JobType.numCores = 4
    config.JobType.allowUndistributedCMSSW = True

    config.section_("Debug")
    config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

    config.section_("Data")
    config.Data.inputDataset = None
    config.Data.splitting = ''
    config.Data.ignoreLocality = options.remote
    config.Data.publication = True if not options.test_only else False
    if options.test_only:
        config.Data.totalUnits = 1
    config.Data.publishDBS = 'phys03'

    config.section_("Site")
    if options.remote:
        config.Site.whitelist = options.whitelist.split(',')
    config.Site.storageSite = options.storageSite

    if options.dcms:
        config.section_("User")
        config.User.voGroup = "dcms"

    print('Using config ' + options.cfg)
    print('Writing to directory ' + options.dir)

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print('Cannot execute commend')
            print(hte.headers)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []
    for ijob in jobsLines :
        s = ijob.rstrip()
        if (len(s)==0 or s[0][0]=='#'): continue
        s = ijob.rstrip()
        jobs.append( s )
        print('  --> added ' + s)

    for ijob, job in enumerate(jobs) :

        ptbin = job.split('/')[1]
        cond = job.split('/')[2]
        datatier = job.split('/')[3]
        requestname = ptbin + '_' + cond
        if len(requestname) > 93:
            requestname = ''.join((requestname[:93-len(requestname)]).split('_')[:-1])
            if 'ext' in cond and not 'ext' in requestname:
                requestname = requestname + '_' + cond.split('_')[-1]
        print('requestname = ', requestname)
        config.General.requestName = requestname
        config.Data.inputDataset = job
        config.Data.outputDatasetTag = re.sub(r'MiniAOD[v]?[0-9]?', options.extension, cond) if cond.startswith('RunII') else cond+'_'+options.extension
        print(config.Data.outputDatasetTag)
        config.Data.outLFNDirBase = options.out 
        config.Data.allowNonValidInputDataset = True
        if datatier == 'MINIAODSIM':
        #   config.Data.splitting = 'FileBased'
        #   config.Data.unitsPerJob = 10
            config.Data.splitting = 'Automatic'
        elif datatier == 'MINIAOD':
          config.Data.splitting = 'LumiBased'
          config.Data.lumiMask = options.lumiMask
          config.Data.unitsPerJob = 50 #10  - at 100, some jobs will run out of time
          config.JobType.maxJobRuntimeMin = 2750 #Default is 1315; 2750 minutes guaranteed to be available; Max I have used is 9000
        print('Submitting ' + config.General.requestName + ', dataset = ' + job)
        print('Configuration :')
        print(config)
        try :
            from multiprocessing import Process

            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
            #submit(config)
        except :
            print('Not submitted.')



if __name__ == '__main__':
    main()
