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


# MC (2016, preVFP):
cmsDriver.py nano_mc_2016_ULPreVFP --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_mcRun2_asymptotic_preVFP_v9  --era Run2_2016,run2_nanoAOD_106Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein  /store/mc/RunIISummer19UL16MiniAOD/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v2/100000/EF482706-F402-E648-AB6A-D5B2700C8023.root --fileout file:nano_mc2016pre.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS  $NO_EXEC $PARALLEL

# MC (2016, postVFP):
cmsDriver.py nano_mc_2016_ULPostVFP --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_mcRun2_asymptotic_v15  --era Run2_2016,run2_nanoAOD_106Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein /store/mc/RunIISummer19UL16MiniAOD/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v2/100000/EF482706-F402-E648-AB6A-D5B2700C8023.root --fileout file:nano_mc2016post.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS  $NO_EXEC $PARALLEL

# Data (2016):
cmsDriver.py nano_data_2016_UL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_dataRun2_v32   --era Run2_2016,run2_nanoAOD_106Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein /store/data/Run2016E/JetHT/MINIAOD/21Feb2020_UL2016_HIPM-v1/40000/BA1C74E6-FBB0-AE40-982F-8C7D3E73E0DC.root --fileout file:nano_data2016.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS  $NO_EXEC $PARALLEL

# MC (2017):
cmsDriver.py nano_mc_2017_UL --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_mc2017_realistic_v8   --era Run2_2017,run2_nanoAOD_106Xv1  \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein /store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/BFAAC85A-F5C5-8843-8D2A-76A9E873E24B.root --fileout file:nano_mc2017.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS  $NO_EXEC $PARALLEL

# Data (2017):
cmsDriver.py nano_data_2017_UL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_dataRun2_v32    --era Run2_2017,run2_nanoAOD_106Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein  /store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/50000/421E15F6-C6DA-D848-9537-FEC70D67C61C.root --fileout file:nano_data2017.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS  $NO_EXEC $PARALLEL

# MC (2018):
cmsDriver.py nano_mc_2018_UL --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_upgrade2018_realistic_v15_L1v1   --era Run2_2018,run2_nanoAOD_106Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein /store/mc/RunIISummer19UL18MiniAOD/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2310000/C42FF9DE-BC2E-EA48-A79B-131083B72BFC.root --fileout file:nano_mc2018.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS  $NO_EXEC $PARALLEL

# Data (2018):
cmsDriver.py nano_data_2018_UL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 106X_dataRun2_v32    --era Run2_2018,run2_nanoAOD_106Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 100 --filein /store/data/Run2018A/JetHT/MINIAOD/12Nov2019_UL2018-v2/100000/BBE577F0-95A1-A542-B86C-58FD46079404.root --fileout file:nano_data2018.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS  $NO_EXEC $PARALLEL
