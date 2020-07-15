# NanoAODJMAR

This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for advance developments of jet algorithms. 
The current content of this branch can be seen [here](http://algomez.web.cern.ch/algomez/testWeb/PFnanoplusBtag_content_106x_v01.html) and the size [here](http://algomez.web.cern.ch/algomez/testWeb/PFnanoplusBtag_size_106x_v01.html).

## Recipe

For **UL** 2016, 2017 and 2018 data and MC **NanoAODv6** according to the [XPOG](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv6) and [PPD](https://twiki.cern.ch/twiki/bin/view/CMS/PdmVLegacy2017Analysis) recommendations:

```
cmsrel  CMSSW_10_6_14
cd  CMSSW_10_6_14/src
cmsenv
git cms-addpkg PhysicsTools/NanoAOD
git cms-addpkg RecoBTag/Combined
git cms-merge-topic andrzejnovak:DDXV2_106
git clone https://github.com/cms-data/RecoBTag-Combined.git RecoBTag/Combined/data
git clone https://github.com/cms-jet/NanoAODJMAR.git PhysicsTools/NanoAODJMAR -b 106x_v01
scram b -j 10
cd PhysicsTools/NanoAODJMAR/test
```
Note: This configuration has been tested for this combination of CMSSW release, global tag, era and dataset. When running over a new dataset you should check with [the nanoAOD workbook twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Running_on_various_datasets_from) to see if the era modifiers in the CRAB configuration files are correct. The jet correction versions are taken from the global tag.

## Local Usage:

UL2017 MC:
```
cmsRun nano106X_on_mini106X_2017_mc_NANO.py
```

UL2017 DATA:
```
cmsRun nano106X_on_mini106X_2017_data_NANO.py
```

## Submission to CRAB

The `submit_all.py` script submit several crab jobs depending on the python config (specify with `-c`). A list of datasets can be included in a txt file after option `-f`, your T2/T3 site should be specify with `-s` and `-d` will create a folder for the crab jobs. For data, `-l` is for the JSON file. For example:

 * For MC:
```
python submit_all.py -c nano106X_on_mini106X_2017_mc_NANO.py -f datasets/MC_UL2017.txt -d NANOUL2017MC -s T3_US_FNALLPC
```

 * For Data:
```
python submit_all.py -c nano106X_on_mini106X_2017_data_NANO.py -f datasets/SingleMuon_UL2017.txt -d NANOUL2017 -s T3_US_FNALLPC -l /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt  
```

## Documenting the Extended NanoAOD Samples

Please document the input and output datasets on the following [google spreadsheet](https://docs.google.com/spreadsheets/d/1puQgU7gn7qeU2SfL_EUkKTG1X_EZ3SbVjiGC-QlkJCI/edit#gid=0). For the MC, the number of events can be found by looking up the output dataset in DAS, or with the `report` option of crab. For the data, you will need to run brilcalc to get the total luminosity of the dataset. See the instructions below. 

## How to create python files using cmsDriver

All the previous python config files were produced with `cmsDriver.py`. Two imporant parameters that one needs to verify in the central nanoAOD documentation are `--conditions` and `--era`. Then, an example of how to create those file, if needed, is shown below:

```
cmsDriver.py nano106X_on_mini106X_2017_mc -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mc2017_realistic_v6 --era Run2_2017 --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 8 --customise PhysicsTools/NanoAODJMAR/nano_jmar_cff.JMARnano_customizeMC -n 100 --filein /store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6_ext2-v2/40000/DC8FEA2E-191B-834F-96C3-53FCF0291ADB.root
```

## How to create website with nanoAOD content

To create nice websites like [this one](http://algomez.web.cern.ch/algomez/testWeb/JMECustomNano102x_mc_v01.html#Jet) with the content of nanoAOD, use the `inspectNanoFile.py` file from the `PhysicsTools/nanoAOD` package as:
```
python PhysicsTools/NanoAOD/test/inspectNanoFile.py NANOAOD.root -s website_with_collectionsize.html -d website_with_collectiondescription.html
```

## Running brilcalc
These are condensed instructions from the lumi POG TWiki (https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM). Also see the brilcalc quickstart guide: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BrilcalcQuickStart.

Note: brilcalc should be run on lxplus. It does not work on the lpc.

Instructions:

1.) Add the following lines to your .bashrc file (or equivalent for your shell). Don't forget to source this file afterwards!

    export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
    export PATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$HOME/.local/bin:$PATH
    
2.) Install brilws:

    pip install --install-option="--prefix=$HOME/.local" brilws
    
3.) Get the json file for your output dataset. In the area in which you submitted your jobs:

    crab report -d [your crab directory]
    
The processedLumis.json file will tell you which lumi sections you successfully ran over. The lumi sections for incomplete, failed, or unpublished jobs are listed in notFinishedLumis.json, failedLumis.json, and notPublishedLumis.json. More info can be found at https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3Commands#crab_report.
    
4.) Run brilcalc on lxplus:

    brilcalc lumi -i processedLumis.json -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -b "STABLE BEAMS"
    
The luminosity of interest will be listed under "totrecorded(/fb)." You can also run this over the other previously mentioned json files.
    
Note: '-b "STABLE BEAMS"' is optional if you've already run over the golden json. 
        Using the normtag is NOT OPTIONAL, as it defines the final calibrations and detectors that are used for a given run.
