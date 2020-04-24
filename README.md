# NanoAODJMAR

This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for advance developments of jet algorithms.
The current content of this branch can be seen [here](http://algomez.web.cern.ch/algomez/testWeb/JMARNanoAODv6_102X_v02_size.html).
This format can be used with [fastjet](http://fastjet.fr) directly.

## Recipes

Note: These configurations have been tested for these combinations of CMSSW releases, global tag, era and dataset. When running over a new dataset you should check with [the nanoAOD workbook twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Running_on_various_datasets_from) to see if the era modifiers in the CRAB configuration files are correct. The jet correction versions are taken from the global tag.

### For NANOAODv6

For 2016, 2017 and 2018 data and MC **NANOAODv6** according to the [XPOG recommendations](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv6):

```
cmsrel  CMSSW_10_2_18
cd  CMSSW_10_2_18/src
cmsenv
git cms-addpkg PhysicsTools/NanoAOD
git clone git@github.com:cms-jet/NanoAODJMAR.git -b 102x PhysicsTools/NanoAODJMAR
scram b -j 10
cd PhysicsTools/NanoAODJMAR/test
```

### For NANOAODv7

**This recipe includes the JMECustomNano collections.**
For 2016, 2017 and 2018 data and MC **NANOAODv7** according to the [XPOG recommendations](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv7):

```
cmsrel  CMSSW_10_2_22
cd  CMSSW_10_2_22/src
cmsenv
git cms-init
git cms-addpkg PhysicsTools/NanoAOD
git cms-merge-topic nurfikri89:nanojme_forNanov7_withDeepJetQG
git clone git@github.com:cms-jet/NanoAODJMAR.git -b 102x PhysicsTools/NanoAODJMAR
scram b -j 10
cd PhysicsTools/NanoAODJMAR/test
```

## How to run

All the python files to run the framework are locate in [PhysicsTools/NanoAODJMAR/test/](PhysicsTools/NanoAODJMAR/test/). Then, use cmsRun to the files described below as:
```
cmsRun nano102x_on_mini94x_2016_mc_NANO.py
```

### Local MC Usage:

2016
  * For JMARnano only (nanoAODv6): `nanov6102x_on_mini94x_2016_mc_NANO.py`
  * For JMARnano only (nanoAODv7): `nano102x_on_mini94x_2016_mc_NANO.py`
  * For JMECustomNano (only for nanoAODv7): `JMECustomNano102x_on_mini94x_2016_mc_NANO.py`
  * For JMECustomNano plus JMARnano (only for nanoAODv7): `JMECustomplusJMARnano102x_on_mini94x_2016_mc_NANO.py`

2017
  * For JMARnano only (nanoAODv6): `nanov6102x_on_mini94x_2017_mc_NANO.py`
  * For JMARnano only (nanoAODv7): `nano102x_on_mini94x_2017_mc_NANO.py`

2018
  * For JMARnano only (nanoAODv6): `nanov6102x_on_mini102x_2018_mc_NANO.py`
  * For JMARnano only (nanoAODv7): `nano102x_on_mini102x_2018_mc_NANO.py`

### Local Data usage:

2016
  * For JMARnano only (nanoAODv6): `nanov6102x_on_mini94x_2016_data_NANO.py`
  * For JMARnano only (nanoAODv7): `nano102x_on_mini94x_2016_data_NANO.py`
  * For JMECustomNano (only for nanoAODv7): `JMECustomNano102x_on_mini94x_2016_data_NANO.py`
  * For JMECustomNano plus JMARnano (only for nanoAODv7): `JMECustomplusJMARnano102x_on_mini94x_2016_data_NANO.py`

2017
  * For JMARnano only (nanoAODv6): `nanov6102x_on_mini94x_2017_data_NANO.py`
  * For JMARnano only (nanoAODv7): `nano102x_on_mini94x_2017_data_NANO.py`

2018
  * For JMARnano only (nanoAODv6): `nanov6102x_on_mini102x_2018_data_abc_NANO.py`, `nanov6102x_on_mini102x_2018_data_d_NANO.py`
  * For JMARnano only (nanoAODv7): `nano102x_on_mini102x_2018_data_abc_NANO.py`, `nano102x_on_mini102x_2018_data_d_NANO.py`

### How to create python files using cmsDriver

All the previous python config files were produced with `cmsDriver.py`. Two imporant parameters that one needs to verify in the central nanoAOD documentation are `--conditions` and `--era`. Then, an example of how to create those file, if needed, is shown below:

```
cmsDriver.py JMECustomNano_plusJMARnano102x_on_mini94x_2016_mc --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_mcRun2_asymptotic_v8 --step NANO --era Run2_2016,run2_nanoAOD_94X2016 --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False))) \n from PhysicsTools.NanoAOD.custom_jme_cff import PrepJMECustomNanoAOD_MC; PrepJMECustomNanoAOD_MC(process)\n" -n 100 --filein /store/mc/RunIISummer16MiniAODv3/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/120000/086BBC76-7AEA-E811-AA5A-6CC2173D9FB0.root --nThreads 2  --customise PhysicsTools/NanoAODJMAR/nano_jmar_cff.JMARnano_customizeMC
```
This complete example includes the JMECustomNano and the JMARnano collections.


### Submission to CRAB

```
python submit_all.py -c nano102x_on_mini94x_2016_mc_NANO.py  -f datasets/2016mc_miniAODv3_DY.txt  -d NANO2016MC

python submit_all.py -c nano102x_on_mini94x_2017_mc_NANO.py -f datasets/2017mc_miniAODv2_DY.txt  -d NANO2017MC

python submit_all.py -c nano102x_on_mini102x_2018_mc_NANO.py -f datasets/2018mc_DY.txt  -d NANO2018MC


python submit_all.py -c nano102x_on_mini94x_2016_data_NANO.py -f datasets/2016data_17Jul2018.txt -d NANO2016 

python submit_all.py -c nano102x_on_mini94x_2017_data_NANO.py  -f datasets/2017data_31Mar2018.txt  -d NANO2017 

python submit_all.py -c nano102x_on_mini102x_2018_data_abc_NANO.py  -f  datasets/2018data_17Sep2018.txt  -d NANO2018 

python submit_all.py -c nano102x_on_mini102x_2018_data_d_NANO.py  -f datasets/datasets_2018D.txt -d NANO2018 

```
The lumimask for data is hardcoded in the `submit_all.py` script.

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

<!--
## Documenting the Extended NanoAOD Samples

Please document the input and output datasets on the following twiki: https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1. For the MC, the number of events can be found by looking up the output dataset in DAS. For the data, you will need to run brilcalc to get the total luminosity of the dataset. See the instructions below. 

-->
