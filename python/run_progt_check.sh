#!/bin/bash
# 2026-11-26 byte check of the ATM_PROGNOSTIC_T port (user: "port ATM_PROGNOSTIC_T to ATHAD_PERID and COND"):
#   run_pt_head  binary before the port (HEAD 29a0311), 2 it, 1 thread
#   run_pt_new   binary with the port, knob unset, 2 it, 1 thread      -> out/ must equal run_pt_head
#   run_pt_on    with the port, ATM_PROGNOSTIC_T=1 ATM_CONV_ADJ_MOIST=1, 4 it, 2 threads (smoke test)
cd "$(dirname "$0")"; H=$1
mk(){ rm -rf $1; mkdir -p $1/out; sed -e "s#<nm>[0-9]*</nm>#<nm>$2</nm>#" -e 's#<checkpoint>[0-9]*</checkpoint>#<checkpoint>2</checkpoint>#' \
      -e 's#<panorama_print>[0-9]*</panorama_print>#<panorama_print>100</panorama_print>#' run_long_adj/config.xml > $1/config.xml; }
mk run_pt_head 2; mk run_pt_new 2; mk run_pt_on 4
( cd run_pt_head && OMP_NUM_THREADS=1 $H config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_pt_new  && OMP_NUM_THREADS=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_pt_on   && OMP_NUM_THREADS=2 ATM_PROGNOSTIC_T=1 ATM_CONV_ADJ_MOIST=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
wait
n=0; d=0; for f in run_pt_head/out/*; do b=$(basename "$f"); n=$((n+1)); cmp -s "$f" "run_pt_new/out/$b" || { d=$((d+1)); echo "DIFF $b"; }; done
echo "exit head $(cat run_pt_head/exit) new $(cat run_pt_new/exit) on $(cat run_pt_on/exit); files compared $n, differing $d; new has $(ls run_pt_new/out | wc -l) files"
