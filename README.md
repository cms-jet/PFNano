# PFNano

**You are currently viewing a development branch**  
Uses PUPPI Jets as default for Run3.

Tested with 2022 data (ReRecoCDE and PromptRecoFG), MC for Run3 (Run3Summer22 and Run3Summer22EE, nanoAODv11).

If you are searching for a recipe to run with Run2 samples, please have a look at the master branch (106X).

This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for advance developments of jet algorithms. 
For data, the current full content of this development branch can be seen [here](https://annika-stein.web.cern.ch/PFNano/desc_data2022.html) and the size [here](https://annika-stein.web.cern.ch/PFNano/size_data2022.html). For MC (124X) the description can be accessed [here](https://annika-stein.web.cern.ch/PFNano/desc_mc2022.html) and the size [here](https://annika-stein.web.cern.ch/PFNano/size_mc2022.html), the same for 122X MC: [description](https://annika-stein.web.cern.ch/PFNano/desc_mc2022_122X.html) and [size](https://annika-stein.web.cern.ch/PFNano/size_mc2022_122X.html).
In this version, PFcandidates can be saved for AK4 only, AK8 only, or all the PF candidates. More below.
This format can be used with [fastjet](http://fastjet.fr) directly.

## Recipe

For 2022 data and MC **NanoAOD v11** according to the [XPOG](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv11) and [PPD](https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis) recommendations:

```
cmsrel CMSSW_12_6_0_patch1
cd CMSSW_12_6_0_patch1/src
cmsenv
git clone https://github.com/cms-jet/PFNano.git PhysicsTools/PFNano
cd PhysicsTools/PFNano
git fetch
git switch 12_6_0
cd ../..
scram b -j 10
cd PhysicsTools/PFNano/test
```
Note: When running over a new dataset you should check with [the nanoAOD workbook twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Running_on_various_datasets_from) to see if the era modifiers in the CRAB configuration files are correct. The jet correction versions are taken from the global tag.

## Local Usage:

There are python config files ready to run in `PhysicsTools/PFNano/test/`.

@BTV-Commissioning-Team: the recommended PFNano customization for commissioning (as of October 2022) for data is `PFnano_customizeData_add_DeepJet` and for MC `PFnano_customizeMC_add_DeepJet_and_Truth`, so without the `allPF` tag.

The list of options that are currently implemented inside `pfnano_cff.py` (e.g. for MC) looks like that:
```
process = PFnano_customizeMC(process)
#process = PFnano_customizeMC_add_DeepJet(process)                  ##### DeepJet inputs are added to the Jet collection
#process = PFnano_customizeMC_add_DeepJet_and_Truth(process)        ##### DeepJet inputs as well as a truth branch with fine-grained labels
#process = PFnano_customizeMC_allPF(process)                        ##### PFcands will contain ALL the PF Cands
#process = PFnano_customizeMC_allPF_add_DeepJet(process)            ##### PFcands will contain ALL the PF Cands; + DeepJet inputs for Jets
#process = PFnano_customizeMC_allPF_add_DeepJet_and_Truth(process)  ##### PFcands will contain ALL the PF Cands; + DeepJet inputs + truth labels for Jets
#process = PFnano_customizeMC_AK4JetsOnly(process)                  ##### PFcands will contain only the AK4 jets PF cands
#process = PFnano_customizeMC_AK4JetsOnly_add_DeepJet(process)      ##### PFcands will contain only the AK4 jets PF cands; + DeepJet inputs for Jets
#process = PFnano_customizeMC_AK8JetsOnly(process)                  ##### PFcands will contain only the AK8 jets PF cands
#process = PFnano_customizeMC_noInputs(process)                     ##### No PFcands but all the other content is available.
```
In general, whenever `_add_DeepJet` is specified (does not apply to `AK8JetsOnly` and `noInputs`), the DeepJet inputs are added to the Jet collection. For all other cases that involve adding tagger inputs, only DeepCSV and / or DDX are taken into account as default (= the old behaviour when `keepInputs=True`). Internally, this is handled by selecting a list of taggers, namely choosing from `DeepCSV`, `DeepJet`, and `DDX` (or an empty list for the `noInputs`-case, formerly done by setting `keepInputs=False`, now set `keepInputs=[]`). This refers to a change of the logic inside `pfnano_cff.py` and `addBTV.py`. If one wants to use this new flexibility, one can also define new customization functions with other combinations of taggers. Currently, there are all configurations to reproduce the ones that were available previously, and all configuations that extend the old ones by adding DeepJet inputs. DeepJet outputs, on top of the discriminators already present in NanoAOD, are added in any case where AK4Jets are added, i.e. there is no need to require the full set of inputs to get the individual output nodes / probabilities. The updated description using `PFnano_customizeMC_allPF_add_DeepJet_and_Truth` can be viewed [here](https://annika-stein.web.cern.ch/PFNano/desc_mc2022.html) and the size [here](https://annika-stein.web.cern.ch/PFNano/size_mc2022.html).

### How to create python files using cmsDriver

All python config files were produced with `cmsDriver.py`.

Two imporant parameters that one needs to verify in the central nanoAOD documentation are `--conditions` and `--era`. 
- `--era` options from [WorkBookNanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) or [XPOG](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv11)
- `--conditions` can be found here [PdMV](https://twiki.cern.ch/twiki/bin/view/CMS/PdmV)

@BTV-Commissioning-Team: the recommended PFNano customization for data is `PFnano_customizeData_add_DeepJet` and for MC `PFnano_customizeMC_add_DeepJet_and_Truth`.

<details>
    <summary>Here are five example commands, with which the runnable configs in `test` have been created:</summary>
    
    
```
cmsDriver.py nano_data_2022ABCD --data --eventcontent NANOAOD --datatier NANOAOD --step NANO \
--conditions 124X_dataRun3_v11   --era Run3,run3_nanoAOD_124 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)));process.MessageLogger.cerr.FwkReport.reportEvery=1000;process.NANOAODoutput.fakeNameForCrab = cms.untracked.bool(True)" --nThreads 4 \
-n -1 --filein "/store/data/Run2022C/DoubleMuon/MINIAOD/10Dec2022-v1/2820000/dea1757f-d2ef-467a-9062-714775d00e45.root" --fileout file:nano_data2022ABCD.root \
--customise="PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData_add_DeepJet"  --no_exec
```
<br>
```
cmsDriver.py nano_data_2022E --data --eventcontent NANOAOD --datatier NANOAOD --step NANO \
--conditions 124X_dataRun3_v14   --era Run3,run3_nanoAOD_124 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)));process.MessageLogger.cerr.FwkReport.reportEvery=1000;process.NANOAODoutput.fakeNameForCrab = cms.untracked.bool(True)" --nThreads 4 \
-n -1 --filein "/store/data/Run2022E/Muon/MINIAOD/10Dec2022-v2/2560000/35d857c4-6c18-4bac-b686-b05540331c10.root" --fileout file:nano_data2022E.root \
--customise="PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData_add_DeepJet"  --no_exec
```
<br>
```
cmsDriver.py nano_data_2022FG --data --eventcontent NANOAOD --datatier NANOAOD --step NANO \
--conditions 24X_dataRun3_Prompt_v10   --era Run3,run3_nanoAOD_124 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)));process.MessageLogger.cerr.FwkReport.reportEvery=1000;process.NANOAODoutput.fakeNameForCrab = cms.untracked.bool(True)" --nThreads 4 \
-n -1 --filein "/store/data/Run2022F/Muon/MINIAOD/PromptReco-v1/000/360/381/00000/0736ad9a-2b1d-4375-9493-9e7e01538978.root" --fileout file:nano_data2022FG.root \
--customise="PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData_add_DeepJet"  --no_exec
```
<br>    
```    
cmsDriver.py nano_mc_Run3 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 126X_mcRun3_2022_realistic_v2   --era Run3,run3_nanoAOD_124 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)));process.MessageLogger.cerr.FwkReport.reportEvery=1000;process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)" --nThreads 4 \
-n -1 --filein "/store/mc/Run3Summer22MiniAODv3/QCD_PT-15to20_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_v12-v1/30000/8590bc1e-abd3-4be4-a068-16f4cb6b4994.root" --fileout file:nano_mcRun3.root \
--customise="PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC_add_DeepJet_and_Truth"  --no_exec
```
<br>
    
```    
cmsDriver.py nano_mc_Run3_EE --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 126X_mcRun3_2022_realistic_postEE_v1   --era Run3,run3_nanoAOD_124 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)));process.MessageLogger.cerr.FwkReport.reportEvery=1000;process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)" --nThreads 4 \
-n -1 --filein "/store/mc/Run3Summer22EEMiniAODv3/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v1/2550000/eddaff63-eb30-4155-afdc-3db5b07105b8.root" --fileout file:nano_mcRun3_EE.root \
--customise="PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC_add_DeepJet_and_Truth"  --no_exec
```
    
</details>


## Submission to CRAB

For crab submission a handler script `crabby.py`, a crab baseline template `template_crab.py` and an example 
submission yaml card `card_example_data.yml` are provided. Fill out the individual entries for each new submission, e.g. dataset from DAS. @BTV-Commissioning-Team: this is also the file to put "BTV_Run3_2022_Comm_v2" for the output folder.

- A single campaign (data/mc, year, config, output path) should be configured statically in a copy of `card_example_dataABCD.yml`.
- To submit:
  ```
  source /cvmfs/grid.cern.ch/centos7-umd4-ui-4_200423/etc/profile.d/setup-c7-ui-example.sh
  source /cvmfs/cms.cern.ch/common/crab-setup.sh prod # note: this is new w.r.t. 106X instructions
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  voms-proxy-init --voms cms --valid 192:00
  cd  CMSSW_12_6_0_patch1/src
  cmsenv
  cd PhysicsTools/PFNano/test
  python3 crabby.py -c card_example_dataABCD.yml --make --submit
  ```


  Or alternatively, split creation and submission of config which allows manual inspection before submission:
  ```
  python3 crabby.py -c card_example_dataABCD.yml --make
  ```
  then inspect manually if configuration is correct, and if all is fine:
  ```
  python3 crabby.py -c card_example_dataABCD.yml --submit
  ```
- Add `--test True` to disable publication on otherwise publishable config and produce a single file per dataset


## Processing data

When processing data, a lumi mask should be applied. The so called golden JSON should be applicable in most cases. Should also be checked here https://twiki.cern.ch/twiki/bin/view/CMS/PdmV

 * Golden JSON re-reco
```
# 2022: TBA
```


 * Golden JSON prompt
```
# 2022: /eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json
```


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

Include in `card.yml` for `crabby.py` submission.


## How to create website with nanoAOD content

To create nice websites like [this one](http://algomez.web.cern.ch/algomez/testWeb/JMECustomNano102x_mc_v01.html#Jet) with the content of nanoAOD, use the `inspectNanoFile.py` file from the `PhysicsTools/nanoAOD` package as:
```
python PhysicsTools/NanoAOD/test/inspectNanoFile.py NANOAOD.root -s website_with_collectionsize.html -d website_with_collectiondescription.html
```

## Documenting the Extended NanoAOD Samples

Please document the input and output datasets on the following twiki: https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1. For the MC, the number of events can be found by looking up the output dataset in DAS. For the data, you will need to run brilcalc to get the total luminosity of the dataset. See the instructions below. 


## Running brilcalc
These are condensed instructions from the lumi POG TWiki (https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM). Also see the brilcalc quickstart guide: https://twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart.

Note: brilcalc should be run on lxplus. It does not work on the lpc.

Instructions:

1.) Add the following lines to your .bashrc file (or equivalent for your shell). Don't forget to source this file afterwards!

    export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
    export PATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$HOME/.local/bin:$PATH
    
2.) Install brilws:

    pip install --install-option="--prefix=$HOME/.local" brilws
    
(Optional: upgrade brilws:)

    pip install --user --upgrade brilws
    
3.) Get the json file for your output dataset. In the area in which you submitted your jobs:

    crab report -d [your crab directory]
    
The processedLumis.json file will tell you which lumi sections you successfully ran over. The lumi sections for incomplete, failed, or unpublished jobs are listed in notFinishedLumis.json, failedLumis.json, and notPublishedLumis.json. More info can be found at https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3Commands#crab_report.
    
4.) Run brilcalc on lxplus:
Note: for Run3, there is no PHYSICS normtag available as of Oct 20, 2022 -> use normtag_BRIL

    brilcalc lumi -i processedLumis.json -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json -b "STABLE BEAMS"
    
The luminosity of interest will be listed under "totrecorded(/fb)." You can also run this over the other previously mentioned json files.
    
Note: '-b "STABLE BEAMS"' is optional if you've already run over the golden json. 
        Using the normtag is NOT OPTIONAL, as it defines the final calibrations and detectors that are used for a given run.
