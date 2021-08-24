#!/bin/bash
source make_configs_preUL_ak15.sh
edmConfigDump nano_data_2016_preUL_NANO.py >& dump_nano_data_2016_preUL_NANO.py
edmConfigDump nano_data_2017_preUL_NANO.py >& dump_nano_data_2017_preUL_NANO.py
edmConfigDump nano_data_2018abc_preUL_NANO.py >& dump_nano_data_2018abc_preUL_NANO.py
edmConfigDump nano_data_2018d_preUL_NANO.py >& dump_nano_data_2018d_preUL_NANO.py
edmConfigDump nano_mc_2016_preUL_NANO.py >& dump_nano_mc_2016_preUL_NANO.py
edmConfigDump nano_mc_2017_preUL_NANO.py >& dump_nano_mc_2017_preUL_NANO.py
edmConfigDump nano_mc_2018_preUL_NANO.py >& dump_nano_mc_2018_preUL_NANO.py

source make_configs_preUL_ak15.sh
edmConfigDump nano_data_2016_UL_NANO.py >& dump_nano_data_2016_UL_NANO.py
edmConfigDump nano_data_2017_UL_NANO.py >& dump_nano_data_2017_UL_NANO.py
edmConfigDump nano_data_2018_UL_NANO.py >& dump_nano_data_2018_UL_NANO.py
edmConfigDump nano_mc_2016_ULPreVFP_NANO.py >& dump_nano_mc_2016_ULPreVFP_NANO.py
edmConfigDump nano_mc_2016_ULPostVFP_NANO.py >& dump_nano_mc_2016_ULPostVFP_NANO.py
edmConfigDump nano_mc_2017_UL_NANO.py >& dump_nano_mc_2017_UL_NANO.py
edmConfigDump nano_mc_2018_UL_NANO.py >& dump_nano_mc_2018_UL_NANO.py

diff dump_nano_data_2016_preUL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_data_2016_preUL_NANO.py >& diff_nano_data_2016_preUL_NANO.txt
diff dump_nano_data_2017_preUL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_data_2017_preUL_NANO.py >& diff_nano_data_2017_preUL_NANO.txt
diff dump_nano_data_2018abc_preUL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_data_2018abc_preUL_NANO.py >& diff_nano_data_2018abc_preUL_NANO.txt
diff dump_nano_data_2018d_preUL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_data_2018d_preUL_NANO.py >& diff_nano_data_2018d_preUL_NANO.txt
diff dump_nano_mc_2016_preUL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_mc_2016_preUL_NANO.py >& diff_nano_mc_2016_preUL_NANO.txt
diff dump_nano_mc_2017_preUL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_mc_2017_preUL_NANO.py >& diff_nano_mc_2017_preUL_NANO.txt
diff dump_nano_mc_2018_preUL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_mc_2018_preUL_NANO.py >& diff_nano_mc_2018_preUL_NANO.txt
diff dump_nano_data_2016_UL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_data_2016_UL_NANO.py >& diff_nano_data_2016_UL_NANO.txt
diff dump_nano_data_2017_UL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_data_2017_UL_NANO.py >& diff_nano_data_2017_UL_NANO.txt
diff dump_nano_data_2018_UL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_data_2018_UL_NANO.py >& diff_nano_data_2018_UL_NANO.txt
diff dump_nano_mc_2016_ULPreVFP_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_mc_2016_ULPreVFP_NANO.py >& diff_nano_mc_2016_ULPreVFP_NANO.txt
diff dump_nano_mc_2016_ULPostVFP_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_mc_2016_ULPostVFP_NANO.py >& diff_nano_mc_2016_ULPostVFP_NANO.txt
diff dump_nano_mc_2017_UL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_mc_2017_UL_NANO.py >& diff_nano_mc_2017_UL_NANO.txt
diff dump_nano_mc_2018_UL_NANO.py ../../../../../CMSSW_10_6_20/src/PhysicsTools/PFNano/test/dump_nano_mc_2018_UL_NANO.py >& diff_nano_mc_2018_UL_NANO.txt
