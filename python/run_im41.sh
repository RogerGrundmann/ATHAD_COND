#!/bin/bash
# 2026-11-26 (user: "add the im option and start the test"):
#   run_im_default  default build (ATM_IM 61) after the im option, 2 it, 1 thread -> out/ must equal run_im_ref
#   run_im_ref      binary before the im option (the ATM_PROGNOSTIC_T commit ebdcd57), 2 it, 1 thread
#   run_im41        cli/cond_im41 (make clean && make cond IM=41), 40 it from scratch, 4 threads,
#                   ATM_GRID_PRESSURE=1 ATM_GRID_PTOP=1.68e-5 ATM_GRID_BETA=4.09 (bottom layer ~300 m), t_skin 221.05;
#                   compare with run_gridp1 (the same PTOP at im 61)
# The run folders must not exist yet; the script stops if one does (it deletes nothing).
cd "$(dirname "$0")"; REF=$1
mk(){ mkdir "$1" || exit 1; mkdir "$1/out"
      sed -e "s#<nm>[0-9]*</nm>#<nm>$2</nm>#" -e "s#<checkpoint>[0-9]*</checkpoint>#<checkpoint>$3</checkpoint>#" \
          -e 's#<panorama_print>[0-9]*</panorama_print>#<panorama_print>100</panorama_print>#' run_long_adj/config.xml > "$1/config.xml"; }
mk run_im_ref 2 2; mk run_im_default 2 2; mk run_im41 40 20
sed -i 's#<t_skin>254.0</t_skin>#<t_skin>221.05</t_skin>#' run_im41/config.xml; grep -c '<t_skin>221.05</t_skin>' run_im41/config.xml
( cd run_im_ref     && OMP_NUM_THREADS=1 $REF config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_im_default && OMP_NUM_THREADS=1 ../../cli/cond config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_im41       && OMP_NUM_THREADS=4 ATM_GRID_PRESSURE=1 ATM_GRID_PTOP=1.68e-5 ATM_GRID_BETA=4.09 ../../cli/cond_im41 config.xml > run.log 2>&1; echo $? > exit ) &
wait
n=0; dd=0; for f in run_im_ref/out/*; do n=$((n+1)); cmp -s "$f" "run_im_default/out/$(basename "$f")" || { dd=$((dd+1)); echo "DIFF $(basename "$f")"; }; done
echo "byte check: $n files, $dd differ; exits ref $(cat run_im_ref/exit) default $(cat run_im_default/exit) im41 $(cat run_im41/exit)"
