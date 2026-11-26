#!/bin/bash
# 2026-11-26 (user: "yes, implement (1) and (a)"): byte check of the two ATM_PROGNOSTIC_T fixes with the knob UNSET.
#   run_fx_old  the binary before the fixes ($1), 2 it, 1 thread
#   run_fx_new  cli/$2 with the fixes, 2 it, 1 thread            -> out/ must be byte-identical
# Fix (1): under the knob the radiatively heated t is carried into tn (cAtmosphereModel.cpp, radiation_mode 2).
# Fix (a): under the knob bcRadius pins the lid to t_skin (BC_Atm.h). Creates the folders; deletes nothing.
cd "$(dirname "$0")"; OLD=$1; B=$2
mk(){ mkdir "$1" || exit 1; mkdir "$1/out"
      sed -e "s#<nm>[0-9]*</nm>#<nm>2</nm>#" -e "s#<checkpoint>[0-9]*</checkpoint>#<checkpoint>2</checkpoint>#" \
          -e 's#<panorama_print>[0-9]*</panorama_print>#<panorama_print>100</panorama_print>#' run_long_adj/config.xml > "$1/config.xml"; }
mk run_fx_old; mk run_fx_new
( cd run_fx_old && OMP_NUM_THREADS=1 $OLD config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_fx_new && OMP_NUM_THREADS=1 ../../cli/$B config.xml > run.log 2>&1; echo $? > exit ) &
wait
n=0; dd=0; for f in run_fx_old/out/*; do n=$((n+1)); cmp -s "$f" "run_fx_new/out/$(basename "$f")" || { dd=$((dd+1)); echo "DIFF $(basename "$f")"; }; done
echo "byte check: $n files, $dd differ; exits old $(cat run_fx_old/exit) new $(cat run_fx_new/exit)"
