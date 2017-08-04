DATE=July4

rm -rf /data/uhussain/MichalTaus_${DATE}_hadd/
mkdir -p /data/uhussain/MichalTaus_${DATE}_hadd/

#for dir in Ztt_pre4_pu0 Ztt_pre4_pu140 Ztt_pre4_pu200; do
#for dir in Ztt_pre4_pu200; do

for dir in Ztt_pre4_pu0 Ztt_pre4_pu140; do

    hadd -f /data/uhussain/MichalTaus_${DATE}_hadd/${dir}.root /hdfs/store/user/uhussain/${dir}-MichalRecipe_${DATE}/MichalRecipe_${DATE}*.root

done

# Combined Data Groups
#hadd -f /data/uhussain/MichalTaus_${DATE}_hadd/ICHEPRuns.root /data/uhussain/MichalTaus_${DATE}_hadd/Run2016[B,C,D].root
#hadd -f /data/uhussain/MichalTaus_${DATE}_hadd/AllRuns.root /data/uhussain/MichalTaus_${DATE}_hadd/Run*.root

# MC Samples
#hadd -f /data/uhussain/TauTiming_${DATE}_hadd/DYJets.root /hdfs/store/user/uhussain/tagAndProbe${DATE}DYJets-ConfFile_MC_reHLT_cfg/*.root 
