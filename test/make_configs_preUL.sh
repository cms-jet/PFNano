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


# MC (2016, 94X, MiniAODv3):
cmsDriver.py nano_mc_2016_preUL --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 102X_mcRun2_asymptotic_v8  --era Run2_2016,run2_nanoAOD_94X2016 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 1000 --filein /store/mc/RunIISummer16MiniAODv3/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/10000/60BC7248-5B7D-E911-AEE2-FA163E83A98B.root --fileout file:nano_mc2016.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS  $NO_EXEC $PARALLEL

# Data (2016, 94X, MiniAODv3):
cmsDriver.py nano_data_2016_preUL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 102X_dataRun2_v13   --era Run2_2016,run2_nanoAOD_94X2016 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 1000 --filein /store/data/Run2016F/JetHT/MINIAOD/17Jul2018-v1/20000/9801F194-178E-E811-988A-A0369FE2C210.root --fileout file:nano_data2016.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS  $NO_EXEC $PARALLEL

# MC (2017, 94X, MiniAODv2):
cmsDriver.py nano_mc_2017_preUL --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 102X_mc2017_realistic_v8   --era Run2_2017,run2_nanoAOD_94XMiniAODv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 1000 --filein /store/mc/RunIIFall17MiniAODv2/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/270000/B4131826-50C9-E811-852D-001E675811CC.root --fileout file:nano_mc2017.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS  $NO_EXEC $PARALLEL

# Data (2017, 94X, MiniAODv2):
cmsDriver.py nano_data_2017_preUL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 102X_dataRun2_v13    --era Run2_2017,run2_nanoAOD_94XMiniAODv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 1000 --filein /store/data/Run2017D/JetHT/MINIAOD/31Mar2018-v1/20000/C8217B8D-1C43-E811-B1AF-0CC47A4C8E3C.root --fileout file:nano_data2017.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS  $NO_EXEC $PARALLEL

# MC (2018, 102X):
cmsDriver.py nano_mc_2018_preUL --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 102X_upgrade2018_realistic_v21   --era Run2_2018,run2_nanoAOD_102Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 1000 --filein /store/mc/RunIIAutumn18MiniAOD/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/20000/DFDD5227-3EFE-004C-93E2-F3FA2F8D2389.root --fileout file:nano_mc2018.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeMC$NOINPUTS  $NO_EXEC $PARALLEL

# Data (2018ABC, 102X):
cmsDriver.py nano_data_2018abc_preUL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 102X_dataRun2_v13    --era Run2_2018,run2_nanoAOD_102Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 1000 --filein /store/data/Run2018B/JetHT/MINIAOD/17Sep2018-v1/00000/D88AC7EE-3193-8446-9056-08CB6977C73B.root --fileout file:nano_data2018abc.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS  $NO_EXEC $PARALLEL

# Data (2018D, 102X):
cmsDriver.py nano_data_2018d_preUL --data --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO \
--conditions 102X_dataRun2_Prompt_v16   --era Run2_2018,run2_nanoAOD_102Xv1 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --nThreads 4 \
-n 1000 --filein /store/data/Run2018D/JetHT/MINIAOD/PromptReco-v2/000/325/057/00000/B58EBF0D-5D25-7842-A7FB-8449EE1A2262.root --fileout file:nano_data2018d.root \
--customise PhysicsTools/PFNano/pfnano_cff.PFnano_customizeData$NOINPUTS  $NO_EXEC $PARALLEL
