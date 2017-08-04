#!/bin/sh                                                                                                                                                                                                                                                                                                                     
#voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat farmout_reRunTaus_miniAOD.py>> MichalRecipe_July19.py 
cat submit.py >> MichalRecipe_July19.py


#for dir in QCD-Flat_911_pu0 QCD-Flat_911_pu200 Ztt_911_pu0 Ztt_911_pu200 
for dir in RelValTTbar_300 RelValTTbar_3000 
do 
    echo " "
    echo "====================" $dir "========================="

    rm -r /data/uhussain/$dir-MichalRecipe_July19/
    mkdir /data/uhussain/$dir-MichalRecipe_July19/
    
#make dag dir
    mkdir -p /data/uhussain/$dir-MichalRecipe_July19/dags
    mkdir -p /data/uhussain/$dir-MichalRecipe_July19/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/$1/server?SFN=/hdfs/store/user/uhussain/$dir-MichalRecipe_July19/
    
    farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=$dir.txt --input-files-per-job=1 \
	--submit-dir=/data/uhussain/$dir-MichalRecipe_July19/submit \
	--output-dag-file=/data/uhussain/$dir-MichalRecipe_July19/dags/dag \
	$dir  \
	$CMSSW_BASE  \
	$CMSSW_BASE/src/RecoTauTag/phase2Taus/test/MichalRecipe_July19.py 

done

rm MichalRecipe_July19.py
