# NanoAODJMAR
This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for testing jet algorithms. It takes the AK8 `pat::Jets` and `reco::GenJets` from `MINIAOD`, unpacks their constituents correctly, and writes a NANOAOD flat table of their four vector, pdgid, and charge. This is intended for advanced development using [fastjet](http://fastjet.fr) directly.

## Recipe


For 2016 and 2017 data (NOT TESTED YET):


```
cmsrel  CMSSW_9_4_4
cd  CMSSW_9_4_4/src
cmsenv
git cms-merge-topic cms-nanoAOD:master-94X
git clone https://github.com/UBParker/NanoAODJMAR.git PhysicsTools/NanoAODJMAR
scram b -j 10
cd PhysicsTools/NanoAODJMAR/test
```

For 2018 data:


```
cmsrel  CMSSW_10_2_9
cd  CMSSW_10_2_9/src
cmsenv
git cms-merge-topic cms-nanoAOD:master-102X
git clone https://github.com/UBParker/NanoAODJMAR.git PhysicsTools/NanoAODJMAR
scram b -j 10
cd PhysicsTools/NanoAODJMAR/test
```


## MC Usage:

```
cmsRun test102X_NANO_recluster.py
cmsRun test94X_NANO_recluster.py
cmsRun test80X_NANO_recluster.py
```

## Data usage:
```
cmsRun test_data_102X_NANO_recluster.py
cmsRun test_data_94X_NANO_recluster.py
cmsRun test_data_80X_NANO_recluster.py

```


## Submission to CRAB
```
python submit_all.py -c test80X_NANO_recluster.py -f datasets_2016mc.txt  -d NANO2016MC
python submit_all.py -c test80X_NANO_recluster.py -f datasets_2016mc_ext.txt  -d NANO2016MCEXT
python submit_all.py -c test94X_NANO_recluster.py -f datasets_2017mc.txt  -d NANO2017MC
python submit_all.py -c test94X_NANO_recluster.py -f datasets_2017mc_ext.txt  -d NANO2017MCEXT
python submit_all.py -c test_data_80X_NANO_recluster.py -f datasets_2016rereco.txt -d NANO2016 -l Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
python submit_all.py -c test_data_94X_NANO_recluster.py -f datasets_2017rereco.txt -d NANO2017 -l Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt
python submit_all.py -c test_data_102X_NANO_recluster.py -f datasets_2018.txt -d NANO2018 -l Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt

```

