#!/bin/bash
# ATM_MC_EVAP_LIMIT (default 0) with ATM_MC_QC_DETRAIN: three 20-iteration arms at 8 threads, two 2-iteration byte checks at 1 thread.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATHAD_COND/python
D="ATM_SATADJ_DIAG=1 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1"
( cd run_el_q  && env OMP_NUM_THREADS=8 $D ATM_MC_QC_DETRAIN=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_el_qe && env OMP_NUM_THREADS=8 $D ATM_MC_QC_DETRAIN=1 ATM_MC_EVAP_LIMIT=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_el_e  && env OMP_NUM_THREADS=8 $D ATM_MC_EVAP_LIMIT=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_el_vnew && OMP_NUM_THREADS=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_el_ve   && OMP_NUM_THREADS=1 ATM_MC_EVAP_LIMIT=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
wait
cmpdir(){ n=0; s=0; for f in $(cd $1 && find . -type f); do n=$((n+1)); cmp -s $1/$f $2/$f && s=$((s+1)); done; echo "$1 = $2 : $s/$n files identical"; }
cmpdir run_v_pre/out run_el_vnew/out
cmpdir run_el_vnew/out run_el_ve/out
cmpdir run_el_q/out run_el_qe/out
echo "exits: $(cat run_el_*/exit | tr '\n' ' ')"
