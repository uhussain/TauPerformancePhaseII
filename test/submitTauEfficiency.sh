#!/bin/sh                                                                                                                                                                                                                                                                                                                     
#voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat runTauEfficiency.py > TwoTausEff_911_patch3_July19.py
cat submit.py >> TwoTausEff_911_patch3_July19.py


#for dir in VBFHToTauTau_pu0 HToTauTau_pu0
#for dir in QCD-Flat-pu0 QCD-Flat-pu140 QCD-Flat-pu200 ZEE-RelVal-pu0 ZEE-RelVal-pu140 ZEE-RelVal-pu200 Ztt_911_pu0 Ztt_911_pu140 Ztt_911_pu200
for dir in RelValTTbar_miniAOD_300 RelValTTbar_miniAOD_3000 
do 
    echo " "
    echo "====================" $dir "========================="

    rm -r /data/uhussain/$dir-TwoTausEff_911_patch3_July19/
    mkdir /data/uhussain/$dir-TwoTausEff_911_patch3_July19/
    
#make dag dir
    mkdir -p /data/uhussain/$dir-TwoTausEff_911_patch3_July19/dags
    mkdir -p /data/uhussain/$dir-TwoTausEff_911_patch3_July19/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/$1/server?SFN=/hdfs/store/user/uhussain/$dir-TauEff_July19/
    
    farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=$dir.txt \
	--submit-dir=/data/uhussain/$dir-TwoTausEff_911_patch3_July19/submit \
	--output-dag-file=/data/uhussain/$dir-TwoTausEff_911_patch3_July19/dags/dag \
	$dir  \
	$CMSSW_BASE  \
	$CMSSW_BASE/src/RecoTauTag/phase2Taus/test/TwoTausEff_911_patch3_July19.py 

done

rm TwoTausEff_911_patch3_July19.py
