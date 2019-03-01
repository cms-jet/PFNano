# Test MC production
python submit_all.py -c  test94X_NANO_recluster.py -d 94X_JMARNANO -f  DY_94X_test.txt

python submit_all.py -c  test80X_NANO_recluster.py -d 80X_JMARNANO -f  DY_80X_test.txt

python submit_all.py -c  test102X_NANO_recluster.py -d 102X_JMARNANO -f  DY_102X_test.txt

# Test Data production
python submit_all.py -c  test_data_94X_NANO_recluster.py -d 94X_JMARNANO -f  SingleMuon_94X_test.txt -l Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt

python submit_all.py -c  test_data_80X_NANO_recluster.py -d 80X_JMARNANO -f  SingleMuon_80X_test.txt -l Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt

python submit_all.py -c  test_data_102X_NANO_recluster.py -d 102X_JMARNANO -f  SingleMuon_102X_test.txt -l Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt