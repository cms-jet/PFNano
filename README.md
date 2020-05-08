# NanoAODJMAR

This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for advance developments of jet algorithms. 
The current preliminary content of this development branch can be seen [here](http://algomez.web.cern.ch/algomez/testWeb/JMARNanoAODv6_106X_v00_size.html).
This format can be used with [fastjet](http://fastjet.fr) directly.

## Recipe

**THIS IS A DEVELOPMENT BRANCH dev_106x**

For **UL** 2016, 2017 and 2018 data and MC **NanoAODv6** according to the [XPOG recommendations](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv6):

```
cmsrel  CMSSW_10_6_5
cd  CMSSW_10_6_5/src
cmsenv
git cms-addpkg PhysicsTools/NanoAOD
git cms-merge-topic laurenhay:setmerger_10_6_5
git clone git@github.com:cms-jet/NanoAODJMAR.git -b dev_106x PhysicsTools/NanoAODJMAR
scram b -j 10
cd PhysicsTools/NanoAODJMAR/test
```
Note: This configuration has been tested for this combination of CMSSW release, global tag, era and dataset. When running over a new dataset you should check with [the nanoAOD workbook twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Running_on_various_datasets_from) to see if the era modifiers in the CRAB configuration files are correct. The jet correction versions are taken from the global tag.

## Local Usage:

2017 MC:
```
cmsRun nano106X_on_mini106X_2017_mc_NANO.py
```

2017 DATA:
```
cmsRun nano106X_on_mini106X_2017_data_NANO.py
```

<!--
## Submission to CRAB

```
python submit_all.py -c nano102x_on_mini94x_2016_mc_NANO.py  -f 2016mc_miniAODv3_DY.txt  -d NANO2016MC

python submit_all.py -c nano102x_on_mini94x_2017_mc_NANO.py -f 2017mc_miniAODv2_DY.txt  -d NANO2017MC

python submit_all.py -c nano102x_on_mini102x_2018_mc_NANO.py -f 2018mc_DY.txt  -d NANO2018MC


python submit_all.py -c nano102x_on_mini94x_2016_data_NANO.py -f 2016data_17Jul2018.txt -d NANO2016 -l Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt

python submit_all.py -c nano102x_on_mini94x_2017_data_NANO.py  -f 2017data_31Mar2018.txt  -d NANO2017 -l Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt



python submit_all.py -c nano102x_on_mini102x_2018_data_abc_NANO.py  -f  2018data_17Sep2018.txt  -d NANO2018 -l Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt

python submit_all.py -c nano102x_on_mini102x_2018_data_d_NANO.py  -f datasets_2018D.txt -d NANO2018 -l Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt

```

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
-->
