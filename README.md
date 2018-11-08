# NanoAODJMAR
This is a [NanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) framework for testing jet algorithms. It takes the AK8 `pat::Jets` and `reco::GenJets` from `MINIAOD`, unpacks their constituents correctly, and writes a NANOAOD flat table of their four vector, pdgid, and charge. This is intended for advanced development using [fastjet](http://fastjet.fr) directly.

## Recipe



```
cmsrel CMSSW_9_4_11_cand1
cd CMSSW_9_4_11_cand1/src
git cms-merge-topic cms-nanoAOD:master-94X
git clone https://github.com/cms-jet/NanoAODJMAR.git PhysicsTools/NanoAODJMAR
scram b -j 10
cd PhysicsTools/NanoAODJMAR/test
```


## MC Usage:

```
cmsRun test94X_NANO_recluster.py
cmsRun test80X_NANO.py
```

## Data usage:
```
cmsRun test_data_80X_NANO.py
cmsRun test_data_92X_NANO.py
```


## Submission to CRAB
```
python submit_all.py -c test80X_NANO_recluster.py -f datasets_2016mc.txt  -d NANO2016MC
python submit_all.py -c test80X_NANO_recluster.py -f datasets_2016mc_ext.txt  -d NANO2016MCEXT
python submit_all.py -c test94X_NANO_recluster.py -f datasets_2017mc.txt  -d NANO2017MC
python submit_all.py -c test94X_NANO_recluster.py -f datasets_2017mc_ext.txt  -d NANO2017MCEXT
python submit_all.py -c test_data_80X_NANO_recluster.py -f datasets_2016rereco.txt -d NANO2016 -l Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
python submit_all.py -c test_data_94X_NANO_recluster.py -f datasets_2017rereco.txt -d NANO2017 -l Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt
```

