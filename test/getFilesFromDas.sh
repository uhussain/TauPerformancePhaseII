#usage->
#bash getFilesFromDas.sh /RelValZTT_13/CMSSW_8_1_0_pre15-PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13ESPU200-v1/GEN-SIM-RECO

input=$1
python dascli.py --query="file dataset=$1" --limit=0
