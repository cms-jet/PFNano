# PFNano

This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for advance developments of jet algorithms. 
The current full content of this development branch can be seen [here](http://algomez.web.cern.ch/algomez/testWeb/PFnano_content_v02.html) and the size [here](http://algomez.web.cern.ch/algomez/testWeb/PFnano_size_v02.html).
In this version, PFcandidates can be saved for AK4 only, AK8 only, or all the PF candidates. More below.
This format can be used with [fastjet](http://fastjet.fr) directly.

## Recipe

**THIS IS A DEVELOPMENT BRANCH**

For **UL** 2016, 2017 and 2018 data and MC **NanoAODv8** according to the [XPOG](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv8) and [PPD](https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysisSummaryTable) recommendations:

```
cmsrel  CMSSW_10_6_20 # in principle not a constraint
cd  CMSSW_10_6_20/src
cmsenv
git cms-rebase-topic andrzejnovak:614nosort
git clone https://github.com/cms-jet/PFNano.git PhysicsTools/PFNano
scram b -j 10
cd PhysicsTools/PFNano/test
```
Note: When running over a new dataset you should check with [the nanoAOD workbook twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Running_on_various_datasets_from) to see if the era modifiers in the CRAB configuration files are correct. The jet correction versions are taken from the global tag.

## Local Usage:

There are python config files ready to run in `PhysicsTools/PFNano/test/` for the UL campaign of nanoAODv8, named `nano106Xv8_on_mini106X_201*_data_NANO.py`. Notice that the current version can create 4 types of files depending on the PF candidates content. 
In this files, for simplicity, the 4 options are included but only one is commented out for use. For instance:
```
process = PFnano_customizeMC(process)
#process = PFnano_customizeMC_allPF(process)            ##### PFcands will content ALL the PF Cands
#process = PFnano_customizeMC_AK4JetsOnly(process)      ##### PFcands will content only the AK4 jets PF cands
#process = PFnano_customizeMC_AK8JetsOnly(process)      ##### PFcands will content only the AK8 jets PF cands
#process = PFnano_customizeMC_noInputs(process)         ##### No PFcands but all the other content is available.
```



### How to create python files using cmsDriver

All  python config files were produced with `cmsDriver.py`.

Two imporant parameters that one needs to verify in the central nanoAOD documentation are `--conditions` and `--era`. 
- `--era` options from [WorkBookNanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) or [XPOG](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv8)
- `--conditions` can be found here [PdMV](https://twiki.cern.ch/twiki/bin/view/CMS/PdmV)

A an example set of `cmsDriver.py` commands to create configs can be found in `make_configs_preUL.sh`

```
bash make_configs_preUL.sh  # run to produce configs
bash make_configs_preUL.sh  -e # run to actually execute configs on 1000 events
```

## How to create website with nanoAOD content

To create nice websites like [this one](http://algomez.web.cern.ch/algomez/testWeb/JMECustomNano102x_mc_v01.html#Jet) with the content of nanoAOD, use the `inspectNanoFile.py` file from the `PhysicsTools/nanoAOD` package as:
```
python PhysicsTools/NanoAOD/test/inspectNanoFile.py NANOAOD.root -s website_with_collectionsize.html -d website_with_collectiondescription.html
```

## Submission to CRAB

Samples can be submitted to crab using the `submit_all.py` script. Run with `-h` option to see usage. Example can look like this:

```
python submit_all.py -c nano_config.py -s T2_DE_RWTH -f datasets/text_list.txt  -o /store/user/$USER/PFNano/  --ext test --test -d crab_noinpts

```

For the UL datasets:
```
##python submit_all.py -c nano102x_on_mini94x_2016_mc_NANO.py  -f 2016mc_miniAODv3_DY.txt  -d NANO2016MC

python submit_all.py -c nano106Xv8_on_mini106X_2017_mc_NANO.py -f 2017mc_miniAODv2_DY.txt  -d NANO2017MC

python submit_all.py -c nano106Xv8_on_mini106X_2018_mc_NANO.py -f 2018mc_DY.txt  -d NANO2018MC


##python submit_all.py -c nano102x_on_mini94x_2016_data_NANO.py -f 2016data_17Jul2018.txt -d NANO2016 -l Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt

python submit_all.py -c nano106Xv8_on_mini106X_2017_data_NANO.py -f 2017data_31Mar2018.txt  -d NANO2017 -l /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt 

python submit_all.py -c nano106Xv8_on_mini106X_2018_data_NANO.py -f datasets_2018D.txt -d NANO2018 -l /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt 

```



## Processing data

When processing data, a lumi mask should be applied. The so called golden JSON should be applicable in most cases. Should also be checked here https://twiki.cern.ch/twiki/bin/view/CMS/PdmV

 * Golden JSON, UL
```
# 2017: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt
# 2018: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
#
```

 * Golden JSON, pre-UL
```
# 2016
jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
# 2017 
jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
# 2018
jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
```

Include in crab submission as `--lumiMask jsons/...txt`

## Documenting the Extended NanoAOD Samples

Please document the input and output datasets on the following twiki: https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1. For the MC, the number of events can be found by looking up the output dataset in DAS. For the data, you will need to run brilcalc to get the total luminosity of the dataset. See the instructions below. 


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
