#!/bin/bash
# 2026-11-26 (user: "make im=41 the default for COND and PERID"): byte checks of the new default, 2 it, 1 thread.
#   run_d41_new  cli/<b> (new default build, ATM_IM 41)        = run_d41_old  cli/<b>_im41 (make IM=41 of the same code)
#   run_d61_new  $2 (make IM=61 of the new tree)               = run_d61_old  $1 (the 61-level binary before the flip)
# Usage: run_imdefault.sh <old61-binary> <new IM=61 binary> <b = perid|cond>.  Creates the folders; deletes nothing.
cd "$(dirname "$0")"; OLD61=$1; NEW61=$2; B=$3
mk(){ mkdir "$1" || exit 1; mkdir "$1/out"
      sed -e "s#<nm>[0-9]*</nm>#<nm>2</nm>#" -e "s#<checkpoint>[0-9]*</checkpoint>#<checkpoint>2</checkpoint>#" \
          -e 's#<panorama_print>[0-9]*</panorama_print>#<panorama_print>100</panorama_print>#' run_long_adj/config.xml > "$1/config.xml"; }
for a in run_d41_new run_d41_old run_d61_new run_d61_old; do mk $a; done
( cd run_d41_new && OMP_NUM_THREADS=1 ../../cli/$B config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_d41_old && OMP_NUM_THREADS=1 ../../cli/${B}_im41 config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_d61_new && OMP_NUM_THREADS=1 $NEW61 config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_d61_old && OMP_NUM_THREADS=1 $OLD61 config.xml > run.log 2>&1; echo $? > exit ) &
wait
for p in d41 d61; do n=0; dd=0
  for f in run_${p}_old/out/*; do n=$((n+1)); cmp -s "$f" "run_${p}_new/out/$(basename "$f")" || { dd=$((dd+1)); echo "DIFF $p $(basename "$f")"; }; done
  echo "$p: $n files, $dd differ; exits new $(cat run_${p}_new/exit) old $(cat run_${p}_old/exit); levels $(grep -cE '^ +[0-9]+ +[0-9.]+ +[0-9.]+ +[0-9.]+ +[0-9.]+ +[0-9.]+' run_${p}_new/run.log | head -1)"
done
