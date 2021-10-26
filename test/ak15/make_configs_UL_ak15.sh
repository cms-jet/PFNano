#!/bin/bash

# Super overblown parser
PARAMS=""
# Set defaults
NO_EXEC="--no_exec "
# Parser
while (( "$#" )); do
  # Parse
  case "$1" in
    -e|--exec)
      NO_EXEC=""
      shift
      ;;
    -p|--parallel)
      PARALLEL=" &"
      shift
      ;;
    --noInputs)
      NOINPUTS="_noInputs"
      shift
      ;;
    -b|--my-flag-with-argument)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        MY_FLAG_ARG=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
# set positional arguments in their proper place
eval set -- "$PARAMS"

# Example files
# file dataset=/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM site=T1_US_FNAL_Disk
export FILE_MC_2016="/store/mc/RunIISummer20UL16MiniAODv2/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/280000/560DDA6F-A784-FE43-B3ED-DC961B1DEA0D.root"
# file dataset=/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM site=T1_US_FNAL_Disk
export FILE_MC_2016PREVFP="/store/mc/RunIISummer20UL16MiniAODAPVv2/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/280000/CF2F010D-2073-DF40-BABB-9B767F347919.root"
# file dataset=/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM site=T1_US_FNAL_Disk
export FILE_MC_2017="/store/mc/RunIISummer20UL17MiniAODv2/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/100000/1421B213-2454-7D42-9103-7841A8763110.root"
# file dataset=/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM site=T1_US_FNAL_Disk
export FILE_MC_2018="/store/mc/RunIISummer20UL18MiniAODv2/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/230000/3CC1A90C-F5D0-084D-9C72-0AE3D658B09B.root"

# file dataset=/JetHT/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD site=T1_US_FNAL_Disk
export FILE_DATA_2016="/store/data/Run2016D/JetHT/MINIAOD/HIPM_UL2016_MiniAODv2-v2/70000/EF073836-B5FC-F342-9907-922315729D4A.root"
# file dataset=/JetHT/Run2017D-UL2017_MiniAODv2-v1/MINIAOD site=T1_US_FNAL_Disk
export FILE_DATA_2017="/store/data/Run2017D/JetHT/MINIAOD/UL2017_MiniAODv2-v1/40000/23C67454-5942-FB42-9D3E-63A8E286C6A4.root"
# file dataset=/JetHT/Run2018C-UL2018_MiniAODv2-v1/MINIAOD site=T1_US_FNAL_Disk
export FILE_DATA_2018="/store/data/Run2018C/JetHT/MINIAOD/UL2018_MiniAODv2-v1/50000/F8CC7E78-1F96-2545-8718-61A3E3AA9AB6.root"


# MC (2016, preVFP):
# cmsDriver.py  --python_filename HIG-RunIISummer20UL16NanoAODAPVv9-00527_1_cfg.py --eventcontent NANOEDMAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --fileout file:HIG-RunIISummer20UL16NanoAODAPVv9-00527.root --conditions 106X_mcRun2_asymptotic_preVFP_v11 --step NANO --filein "dbs:/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM" --era Run2_2016_HIPM,run2_nanoAOD_106Xv2 --no_exec --mc -n $EVENTS || exit $? ;
cmsDriver.py nano_mc_2016_ULPreVFP --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_mcRun2_asymptotic_preVFP_v11  --era Run2_2016_HIPM,run2_nanoAOD_106Xv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein  $FILE_MC_2016PREVFP --fileout file:nano_mc2016pre.root \
--customise PhysicsTools/PFNano/ak15/addAK15_cff.setupPFNanoAK15_mc  $NO_EXEC $PARALLEL 
# PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS

# MC (2016, postVFP):
# cmsDriver.py  --python_filename HIG-RunIISummer20UL16NanoAODv9-00093_1_cfg.py --eventcontent NANOEDMAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --fileout file:HIG-RunIISummer20UL16NanoAODv9-00093.root --conditions 106X_mcRun2_asymptotic_v17 --step NANO --filein "dbs:/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM" --era Run2_2016,run2_nanoAOD_106Xv2 --no_exec --mc -n $EVENTS || exit $? ;
cmsDriver.py nano_mc_2016_ULPostVFP --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_mcRun2_asymptotic_v17  --era Run2_2016,run2_nanoAOD_106Xv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein $FILE_MC_2016 --fileout file:nano_mc2016post.root \
--customise PhysicsTools/PFNano/ak15/addAK15_cff.setupPFNanoAK15_mc  $NO_EXEC $PARALLEL 
# PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS

# Data (2016):
cmsDriver.py nano_data_2016_UL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_dataRun2_v35   --era Run2_2016,run2_nanoAOD_106Xv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein $FILE_DATA_2016 --fileout file:nano_data2016.root \
--customise PhysicsTools/PFNano/ak15/addAK15_cff.setupPFNanoAK15_data  $NO_EXEC $PARALLEL
# ,PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS

# MC (2017):
# cmsDriver.py  --python_filename HIG-RunIISummer20UL17NanoAODv9-00678_1_cfg.py --eventcontent NANOEDMAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --fileout file:HIG-RunIISummer20UL17NanoAODv9-00678.root --conditions 106X_mc2017_realistic_v9 --step NANO --filein "dbs:/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM" --era Run2_2017,run2_nanoAOD_106Xv2 --no_exec --mc -n $EVENTS || exit $? ;
cmsDriver.py nano_mc_2017_UL --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_mc2017_realistic_v9   --era Run2_2017,run2_nanoAOD_106Xv2  \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein $FILE_MC_2017 --fileout file:nano_mc2017.root \
--customise PhysicsTools/PFNano/ak15/addAK15_cff.setupPFNanoAK15_mc  $NO_EXEC $PARALLEL 
# PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS

# Data (2017):
cmsDriver.py nano_data_2017_UL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_dataRun2_v35    --era Run2_2017,run2_nanoAOD_106Xv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein  $FILE_DATA_2017 --fileout file:nano_data2017.root \
--customise PhysicsTools/PFNano/ak15/addAK15_cff.setupPFNanoAK15_data  $NO_EXEC $PARALLEL
# ,PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS

# MC (2018):
# cmsDriver.py  --python_filename HIG-RunIISummer20UL18NanoAODv9-00669_1_cfg.py --eventcontent NANOEDMAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --fileout file:HIG-RunIISummer20UL18NanoAODv9-00669.root --conditions 106X_upgrade2018_realistic_v16_L1v1 --step NANO --filein "dbs:/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM" --era Run2_2018,run2_nanoAOD_106Xv2 --no_exec --mc -n $EVENTS || exit $? ;
cmsDriver.py nano_mc_2018_UL --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_upgrade2018_realistic_v16_L1v1   --era Run2_2018,run2_nanoAOD_106Xv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein $FILE_MC_2018 --fileout file:nano_mc2018.root \
--customise PhysicsTools/PFNano/ak15/addAK15_cff.setupPFNanoAK15_mc  $NO_EXEC $PARALLEL 
# PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS

# Data (2018):
cmsDriver.py nano_data_2018_UL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_dataRun2_v35    --era Run2_2018,run2_nanoAOD_106Xv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein $FILE_DATA_2018 --fileout file:nano_data2018.root \
--customise PhysicsTools/PFNano/ak15/addAK15_cff.setupPFNanoAK15_data  $NO_EXEC $PARALLEL
# ,PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS
