#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser

def make_list(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: python submit_all.py -c CONFIG -d DIR -f DATASETS_FILE')

    parser = OptionParser(usage=usage)    
    parser.add_option("-c", "--config", dest="cfg", default="test94X_NANO.py",
        help=("The crab script you want to submit "),
        metavar="CONFIG")
    parser.add_option("-d", "--dir", dest="dir", default="NANO",
        help=("The crab directory you want to use "),
        metavar="DIR")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")
    parser.add_option("-s", "--storageSite", dest="storageSite", default="T3_US_FNALLPC",
        help=("Site"),
        metavar="SITE")
    parser.add_option("-l", "--lumiMask", dest="lumiMask",
        help=("Lumi Mask JSON file"),
        metavar="LUMIMASK")

    (options, args) = parser.parse_args()

    if options.cfg == None or options.dir == None or options.datasets == None or options.storageSite == None:
        parser.error(usage)
    
    return options
    

def main():

    options = getOptions()

    from WMCore.Configuration import Configuration
    config = Configuration()

    from CRABAPI.RawCommand import crabCommand
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
    config.JobType.maxJobRuntimeMin = 2750 #Default is 1315; 2750 minutes guaranteed to be available; Max I have used is 9000

    config.section_("Debug")
    config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

    config.section_("Data")
    config.Data.inputDataset = None
    config.Data.splitting = ''
    #config.Data.unitsPerJob = 1
    config.Data.ignoreLocality = False
    config.Data.publication = True    
    config.Data.publishDBS = 'phys03'

    config.section_("Site")
    config.Site.blacklist = ['T2_IN_TIFR','T2_US_Caltech']
    #config.Site.whitelist = ['T2_US_UCSD','T2_DE_DESY', 'T1_US_FNAL','T2_UK_SGrid_RALPP','T2_PL_Swierk','T2_TW_NCHC','T2_BR_SPRACE']
    config.Site.storageSite = options.storageSite

    print 'Using config ' + options.cfg
    print 'Writing to directory ' + options.dir
    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute commend'
            print hte.headers

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
        print '  --> added ' + s
        
    for ijob, job in enumerate(jobs) :

        ptbin = job.split('/')[1]
        cond = job.split('/')[2]
        datatier = job.split('/')[3]
        requestname = ptbin + '_' + cond
        if len(requestname) > 93: 
            requestname = ''.join((requestname[:93-len(requestname)]).split('_')[:-1])
            if 'ext' in cond and not 'ext' in requestname:
                requestname = requestname + '_' + cond.split('_')[-1]
        #requestname = requestname + '_try3'
        print 'requestname = ', requestname
        config.General.requestName = requestname
        config.Data.inputDataset = job
        config.Data.outputDatasetTag = requestname 
        config.Data.outLFNDirBase    = '/store/group/lpctlbsm/NanoAODJMAR_2019_V1/Production/CRAB/'
        if datatier == 'MINIAODSIM': 
          config.Data.splitting = 'FileBased'
          config.Data.unitsPerJob = 1
        elif datatier == 'MINIAOD': 
          config.Data.splitting = 'LumiBased'
          config.Data.lumiMask = options.lumiMask
          config.Data.unitsPerJob = 50 #10 # 200
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print config
        try :
            from multiprocessing import Process
            
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
            #submit(config)
        except :
            print 'Not submitted.'



if __name__ == '__main__':
    main()            
