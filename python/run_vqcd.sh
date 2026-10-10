#!/bin/bash
# [MC-GP] line and ATM_MC_QC_DETRAIN (default 0): byte check, 2 iterations from scratch, 1 thread.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATHAD_COND/python
( cd run_v_pre && OMP_NUM_THREADS=1 ../../cli/cond_pre config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_v_new && OMP_NUM_THREADS=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_v_on  && OMP_NUM_THREADS=1 ATM_SATADJ_DIAG=1 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_v_qcd && OMP_NUM_THREADS=1 ATM_MC_QC_DETRAIN=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
wait
cmpdir(){ n=0; s=0; for f in $(cd $1 && find . -type f); do n=$((n+1)); cmp -s $1/$f $2/$f && s=$((s+1)); done; echo "$1 = $2 : $s/$n files identical"; }
cmpdir run_v_pre/out run_v_new/out
cmpdir run_v_new/out run_v_on/out
cmpdir run_v_new/out run_v_qcd/out
echo "exits: $(cat run_v_*/exit | tr '\n' ' ')"
