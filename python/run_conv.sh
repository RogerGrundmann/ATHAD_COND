#!/bin/bash
# ATM_CONV_SCHEME: adjust (default since 2026-10-10) against tiedtke. 20 iterations from scratch, 8 threads, instruments on;
# byte check that ATM_CONV_SCHEME=tiedtke reproduces the binary before 2026-10-10's changes, 2 iterations, 1 thread.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATHAD_COND/python
( cd run_ca_adj  && env OMP_NUM_THREADS=8 ATM_SATADJ_DIAG=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_ca_tied && env OMP_NUM_THREADS=8 ATM_SATADJ_DIAG=1 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_CONV_SCHEME=tiedtke ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_ca_vt   && OMP_NUM_THREADS=1 ATM_CONV_SCHEME=tiedtke ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
wait
cmpdir(){ n=0; s=0; for f in $(cd $1 && find . -type f); do n=$((n+1)); cmp -s $1/$f $2/$f && s=$((s+1)); done; echo "$1 = $2 : $s/$n files identical"; }
cmpdir run_v_pre/out run_ca_vt/out
cmpdir run_dg_20/output_dg_20 run_ca_tied/out
cmpdir run_ca_tied/out run_ca_adj/out
echo "exits: $(cat run_ca_*/exit | tr '\n' ' ')"
