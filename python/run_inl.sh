#!/bin/bash
# ATM_MC_INLINE (default 0): byte check of the unset default against the binary before 2026-10-10's changes, 2 iterations, 1 thread.
# The pass-count arms (run_in<mode>_<iter_prec>) were: OMP_NUM_THREADS=3 ATM_MC_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_MC_QC_DETRAIN=1 ATM_MC_EVAP_LIMIT=1
#   ATM_MC_ED_AREA=1 ATM_MC_QVD=1 ATM_MC_INLINE=<1|2>, config_port_new.xml with nm = 2, checkpoint = 2 and <iter_prec> 3 / 4 / 5 / 6.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATHAD_COND/python
( cd run_in_vnew && OMP_NUM_THREADS=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit )
n=0; s=0; for f in $(cd run_v_pre/out && find . -type f); do n=$((n+1)); cmp -s run_v_pre/out/$f run_in_vnew/out/$f && s=$((s+1)); done
echo "run_v_pre/out = run_in_vnew/out : $s/$n files identical; exit $(cat run_in_vnew/exit)"
