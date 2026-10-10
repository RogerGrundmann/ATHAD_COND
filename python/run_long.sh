#!/bin/bash
# THE LONG RUN (2026-10-10): does the convective transport that went with Tiedtke matter?
# Per tree, three arms side by side, 400 iterations from scratch, 8 threads each, checkpoint 50, no restart dumps:
#   run_long_adj    ATM_CONV_SCHEME unset (= adjust, the default since 2026-10-10)
#   run_long_tied   ATM_CONV_SCHEME=tiedtke
#   run_long_tied2  ATM_CONV_SCHEME=tiedtke again -- the noise floor: how far do two identical runs drift apart in 400 iterations
# Trees in sequence: ATHAD_COND (config_port_new.xml), ATHAD (config_chk40.xml), ATHAD_PERID (config_port_new.xml).
# Binaries are copies made before the run (cli/<name>_long), so rebuilding a tree does not touch a run in flight.
R=/home/roger/SynologyDrive/Cloudstation_Notebook
LOG=$R/ATHAD_COND/python/run_long.status
: > $LOG
for spec in "ATHAD_COND:cond_long" "ATHAD:had_long" "ATHAD_PERID:perid_long"; do
  T=${spec%:*}; B=${spec#*:}
  free=$(df --output=avail -BG / | tail -1 | tr -dc 0-9)
  if [ "$free" -lt 12 ]; then echo "$(date +%T) $T NOT STARTED: only ${free} GB free" >> $LOG; break; fi
  echo "$(date +%T) $T started, ${free} GB free" >> $LOG
  cd $R/$T/python
  X=""; [ $T = ATHAD_COND ] && X="ATM_SATADJ_DIAG=1"
  ( cd run_long_adj   && env OMP_NUM_THREADS=8 $X ../../cli/$B config.xml > run.log 2>&1; echo $? > exit ) &
  ( cd run_long_tied  && env OMP_NUM_THREADS=8 $X ATM_CONV_SCHEME=tiedtke ../../cli/$B config.xml > run.log 2>&1; echo $? > exit ) &
  ( cd run_long_tied2 && env OMP_NUM_THREADS=8 $X ATM_CONV_SCHEME=tiedtke ../../cli/$B config.xml > run.log 2>&1; echo $? > exit ) &
  wait
  find run_long_adj run_long_tied run_long_tied2 -name '*.bin' -size +100M -delete
  echo "$(date +%T) $T finished: exits $(cat run_long_adj/exit run_long_tied/exit run_long_tied2/exit | tr '\n' ' '), sizes $(du -sh run_long_adj run_long_tied run_long_tied2 | cut -f1 | tr '\n' ' ')" >> $LOG
done
echo "$(date +%T) ALL DONE, $(df -h / | tail -1 | awk '{print $4}') free" >> $LOG
