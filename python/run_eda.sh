#!/bin/bash
# ATM_MC_ED_AREA (default 0) on top of ATM_MC_QC_DETRAIN=1 ATM_MC_EVAP_LIMIT=1, with ATM_MC_QVD 0 / 1 / 2: 20 iterations, 8 threads; byte check of the unset default, 2 iterations, 1 thread.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATHAD_COND/python
D="ATM_SATADJ_DIAG=1 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_MC_QC_DETRAIN=1 ATM_MC_EVAP_LIMIT=1 ATM_MC_ED_AREA=1"
for q in 0 1 2; do ( cd run_ea$q && env OMP_NUM_THREADS=8 $D ATM_MC_QVD=$q ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) & done
( cd run_ea_vnew && OMP_NUM_THREADS=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
wait
cmpdir(){ n=0; s=0; for f in $(cd $1 && find . -type f); do n=$((n+1)); cmp -s $1/$f $2/$f && s=$((s+1)); done; echo "$1 = $2 : $s/$n files identical"; }
cmpdir run_v_pre/out run_ea_vnew/out
cmpdir run_qv1/out run_ea1/out
cmpdir run_el_qe/out run_ea0/out
echo "exits: $(cat run_ea*/exit | tr '\n' ' ')"
