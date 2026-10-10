#!/bin/bash
# ATM_SATADJ_DIAG / ATM_MC_DIAG / ATM_MC_CAP_DIAG: print-only check. 3 arms x 4 it at 1 thread, 1 arm x 20 it at 8 threads.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATHAD_COND/python
( cd run_dg_pre && OMP_NUM_THREADS=1 ../../cli/cond_pre config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_dg_new && OMP_NUM_THREADS=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_dg_on  && OMP_NUM_THREADS=1 ATM_SATADJ_DIAG=1 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_dg_20  && OMP_NUM_THREADS=8 ATM_SATADJ_DIAG=1 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
wait
cmpdir(){ n=0; s=0; for f in $(cd $1 && find . -type f); do n=$((n+1)); cmp -s $1/$f $2/$f && s=$((s+1)) || echo "   differs: $f"; done; echo "$1 = $2 : $s/$n files identical"; }
cmpdir run_dg_pre/output_dg_pre run_dg_new/output_dg_new
cmpdir run_dg_new/output_dg_new run_dg_on/output_dg_on
echo "exits: $(cat run_dg_*/exit | tr '\n' ' ')"
