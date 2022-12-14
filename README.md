# PFNano

This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for advanced developments of jet algorithms. 

The repository consists of multiple branches which are each dedicated to specific releases of [CMSSW](https://github.com/cms-sw/cmssw). If you came here to run over Run3 samples, please checkout the most up-to-date 12_4_8 branch (e.g. from the dropdown menu above). The branch you are viewing right now is optimized to run over Run2 UL samples, using the 106X release cycle.


The current full content of this development branch can be seen [here](https://annika-stein.web.cern.ch/PFNano/AddDeepJetTagInfo_desc.html) and the size [here](https://annika-stein.web.cern.ch/PFNano/AddDeepJetTagInfo_size.html).
In this version, PFcandidates can be saved for AK4 only, AK8 only, or all the PF candidates. More below.
This format can be used with [fastjet](http://fastjet.fr) directly.

## Recipe

**THIS IS A DEVELOPMENT BRANCH**

For **UL** 2016, 2017 and 2018 data and MC **NanoAODv8/v9** according to the [XPOG](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/home) and [PPD](https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysisSummaryTable) recommendations:

### Prerequisites
If not already done in your `.bashrc` or similar:
```
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh or /cvmfs/grid.cern.ch/centos7-umd4-ui-4_200423/etc/profile.d/setup-c7-ui-example.sh
source /cvmfs/cms.cern.ch/common/crab-setup.sh prod
source /cvmfs/cms.cern.ch/cmsset_default.sh
```
### Setup PFNano
```
cmsrel  CMSSW_10_6_30
cd  CMSSW_10_6_30/src
cmsenv
git clone https://github.com/AnnikaStein/PFNano.git PhysicsTools/PFNano #change once it's officially included: git clone https://github.com/cms-jet/PFNano.git PhysicsTools/PFNano
cd PhysicsTools/PFNano
git fetch
git switch ParT_106X
cd ../..
scram b -j 10
cd PhysicsTools/PFNano/test
voms-proxy-init --voms cms:/cms/dcms --valid 192:00 # cms:/cms/dcms if you have a German grid certificate, prio to run at German sites; use cms only otherwise
```
Note: When running over a new dataset you should check with [the nanoAOD workbook twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Running_on_various_datasets_from) to see if the era modifiers in the CRAB configuration files are correct. The jet correction versions are taken from the global tag.

## Local Usage:

There are python config files ready to run in `PhysicsTools/PFNano/test/` for the UL campaign of nanoAODv8(v9), named `nano106Xv8_on_mini106X_201*_data_NANO.py` (`nano_data_2017_ULv2_allPF_ParT_NANO.py`, similar for MC). Notice that the current version can create different types of files depending on the PF candidates content. Run a configuration via `cmsRun nano_data_2017_ULv2_allPF_ParT_NANO.py` to test on one file locally. If that works, configure your crab submission (more below).

### Different use cases:

New since December 2022: Particle Transformer inputs. Particle Transformer variables can be added by specifying one of the suiting customizations `*add_DeepJet_ParT*`.

Now the list of options that are currently implemented inside `pfnano_cff.py` (e.g. for MC) looks like that:
```
process = PFnano_customizeMC(process)
#process = PFnano_customizeMC_add_DeepJet(process)                  ##### DeepJet inputs are added to the Jet collection
#process = PFnano_customizeMC_add_DeepJet_and_Truth(process)        ##### DeepJet inputs as well as a truth branch with fine-grained labels
#process = PFnano_customizeMC_add_DeepJet_ParT_and_Truth(process)        ##### DeepJet & ParT inputs as well as a truth branch with fine-grained labels
#process = PFnano_customizeMC_allPF(process)                        ##### PFcands will contain ALL the PF Cands
#process = PFnano_customizeMC_allPF_add_DeepJet(process)            ##### PFcands will contain ALL the PF Cands; + DeepJet inputs for Jets
#process = PFnano_customizeMC_allPF_add_DeepJet_and_Truth(process)  ##### PFcands will contain ALL the PF Cands; + DeepJet inputs + truth labels for Jets
#process = PFnano_customizeMC_allPF_add_DeepJet_ParT_and_Truth(process)  ##### PFcands will contain ALL the PF Cands; + DeepJet & ParT inputs + truth labels for Jets
#process = PFnano_customizeMC_AK4JetsOnly(process)                  ##### PFcands will contain only the AK4 jets PF cands
#process = PFnano_customizeMC_AK4JetsOnly_add_DeepJet(process)      ##### PFcands will contain only the AK4 jets PF cands; + DeepJet inputs for Jets
#process = PFnano_customizeMC_AK8JetsOnly(process)                  ##### PFcands will contain only the AK8 jets PF cands
#process = PFnano_customizeMC_noInputs(process)                     ##### No PFcands but all the other content is available.
```
In general, whenever `_add_DeepJet` is specified (does not apply to `AK8JetsOnly` and `noInputs`), the DeepJet inputs are added to the Jet collection. For all other cases that involve adding tagger inputs, only DeepCSV and / or DDX are taken into account as default (= the old behaviour when `keepInputs=True`). Internally, this is handled by selecting a list of taggers, namely choosing from `DeepCSV`, `DeepJet`, and `DDX` (or an empty list for the `noInputs`-case, formerly done by setting `keepInputs=False`, now set `keepInputs=[]`). This refers to a change of the logic inside `pfnano_cff.py` and `addBTV.py`. If one wants to use this new flexibility, one can also define new customization functions with other combinations of taggers. Currently, there are all configurations to reproduce the ones that were available previously, and all configuations that extend the old ones by adding DeepJet inputs. DeepJet outputs, on top of the discriminators already present in NanoAOD, are added in any case where AK4Jets are added, i.e. there is no need to require the full set of inputs to get the individual output nodes / probabilities. The updated description using `PFnano_customizeMC_add_DeepJet` can be viewed [here](https://annika-stein.web.cern.ch/PFNano/AddDeepJetTagInfo_desc.html) and the size [here](https://annika-stein.web.cern.ch/PFNano/AddDeepJetTagInfo_size.html).

The latest addition before moving to the Run3 recipe was the inclusion of a fine-grained truth flavour branch, which encodes various bits known from DeepNTuples into one integer. It can be activated for simulated samples by adding the `_and_Truth` flag to the customizer, i.e. using `PFnano_customizeMC_add_DeepJet_and_Truth`.

### How to create python files using cmsDriver

(You can skip this step, if the existing configurations are sufficient for your use case.)

All python config files were produced with `cmsDriver.py`.

Two imporant parameters that one needs to verify in the central nanoAOD documentation are `--conditions` and `--era`. 
- `--era` options from [WorkBookNanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) or [XPOG](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/home)
- `--conditions` can be found here [PdMV](https://twiki.cern.ch/twiki/bin/view/CMS/PdmV)

Pre UL `cmsRun` python config files are generated by running `make_configs_preUL.sh`

```
bash make_configs_preUL.sh  # run to only produce configs
bash make_configs_preUL.sh  -e # run to actually execute configs on 1000 events
```

**UL `cmsRun` python config files are generated by running `make_configs_UL.sh`**

```
bash make_configs_UL.sh  # run to only produce configs
bash make_configs_UL.sh  -e # run to actually execute configs on 1000 events
```

## Submission to CRAB

For crab submission a handler script `crabby.py`, a crab baseline template `template_crab.py` and an example 
submission yaml card `card_example.yml` are provided.

- A single campaign (data/mc, year, config, output path) should be configured statically in a copy of `card_example.yml`.
- To submit:
  ```
  source crab.sh
  python crabby.py -c card.yml --make --submit
  ```
- `--make` and `--submit` calls are independent, allowing manual inspection of submit configs
- Add `--test True` to disable publication on otherwise publishable config and produce a single file per dataset

<details>
    <summary>If experiencing problems with crab submission using the above instructions, e.g. on NAF-DESY</summary>
    
    
    ```
    source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh or /cvmfs/grid.cern.ch/centos7-umd4-ui-4_200423/etc/profile.d/setup-c7-ui-example.sh
    source /cvmfs/cms.cern.ch/common/crab-setup.sh prod
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    < navigate to CMSSW_X_Y_Z/src >
    cmsenv
    cd PhysicsTools/PFNano/test
    ```
    
    
    and proceed with `crabby.py` as explained above, activate proxy for submission to be successful.
</details>

<details>
    <summary>Useful commands to get paths to individual processed files</summary>
    
    This is to get a list of files stored at the respective site.
    ```
    xrdfs [insert redirector to site] ls /store/path/to/your/crab/output/serialnumber > filelist.txt
    # (example for T2_DE_RWTH: redirector would be grid-cms-xrootd.physik.rwth-aachen.de)
    # ( if there is more than one serial number (more than 1k files processed) repeat command but append to textfile using >> instead of > )
    # ( clean textfile for log entries )
    # ( then append the redirector (needs modification by you for specific site) using this helper: )
    python dataset_paths.py name_of_txt_file T2_DE_RWTH
    ```
    
</details>

<details>
  <summary>Deprecated submission.</summary>
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
</details>


## Processing data

When processing data, a lumi mask should be applied. The so called golden JSON should be applicable in most cases. Should also be checked here https://twiki.cern.ch/twiki/bin/view/CMS/PdmV

 * Golden JSON, UL
```
# 2ß16: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
# 2017: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt
# 2018: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
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

Include in `card.yml` for `crabby.py` submission. (In deprecated interactive submissiong add `--lumiMask jsons/...txt`)


## How to create website with nanoAOD content

To create nice websites like [this one](http://algomez.web.cern.ch/algomez/testWeb/JMECustomNano102x_mc_v01.html#Jet) with the content of nanoAOD, use the `inspectNanoFile.py` file from the `PhysicsTools/nanoAOD` package as:
```
python PhysicsTools/NanoAOD/test/inspectNanoFile.py NANOAOD.root -s website_with_collectionsize.html -d website_with_collectiondescription.html
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
