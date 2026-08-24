#!/usr/bin/env python3
"""Generate one config per ATM_PRESS_SWEEPS arm.

The sweep count is an ENV VAR, not a config field, so the arms differ ONLY in
output_path: every physics setting is byte-identical across the scan and the
comparison isolates the pressure-solver sweep count.

Moist physics is off automatically (moist_phys_start_iter > nm), which also keeps
item 61's trigger-set chaos amplifier out of the comparison.

usage: mksweeps.py <nm> <sweeps> [<sweeps> ...]
"""
import sys, re

def make(nm, s):
    src = open("config_cond.xml").read()
    subs = [
        (r"<output_path>[^<]*</output_path>", f"<output_path>output_sweeps_{s}/</output_path>"),
        (r"<nm>[0-9]*</nm>",                  f"<nm>{nm}</nm>"),
        (r"<moist_phys_start_iter>[0-9]*</moist_phys_start_iter>",
                                              f"<moist_phys_start_iter>{nm+1}</moist_phys_start_iter>"),
        (r"<restart_stride>[0-9-]*</restart_stride>",   "<restart_stride>0</restart_stride>"),
        (r"<checkpoint_save_iter>[0-9-]*</checkpoint_save_iter>", "<checkpoint_save_iter>-1</checkpoint_save_iter>"),
        (r"<panorama_print>[0-9-]*</panorama_print>",   f"<panorama_print>{nm+1}</panorama_print>"),
        (r"<checkpoint>[0-9]*</checkpoint>",            "<checkpoint>20</checkpoint>"),
        (r"<diagnostic_stride>[0-9]*</diagnostic_stride>", "<diagnostic_stride>20</diagnostic_stride>"),
    ]
    out = src
    for pat, rep in subs:
        out, n = re.subn(pat, rep, out)
        if n != 1:
            raise SystemExit(f"pattern {pat!r} matched {n} times")
    fn = f"config_sweeps_{s}.xml"
    open(fn, "w").write(out)
    print(f"{fn}  ATM_PRESS_SWEEPS={s}  nm={nm}")

if __name__ == "__main__":
    nm = int(sys.argv[1])
    for s in sys.argv[2:]:
        make(nm, int(s))
