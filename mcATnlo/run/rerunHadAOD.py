#!/usr/bin/env python

from glob import glob
import subprocess
import os
import time

dataset = "tt01j_bsbar_2l_FxFx"

torun = glob("/xrootd/store/user/iawatson/%s/MINIAODSIM/*.root" % dataset)
max_processes = 8
run_on_batch = True
delete_had = True
delete_nano = False

command = "./rerun.sh"
processes = set()

print "---- RUNNING ----"

if delete_had:
    os.system("rm -f /xrootd/store/user/iawatson/%s/HADAOD/*.root" % dataset)

if delete_nano:
    os.system("rm -f /xrootd/store/user/iawatson/%s/NANOAOD/*.root" % dataset)

for name in torun:
    if run_on_batch:
        n = name.split('/')[-1].split('.')[0]
        script = 'condor_submit -append arguments="%s %s" rerun.jds' % (dataset, n)
        print script
        os.system(script)
    else:
        # Run locally
        processes.add(subprocess.Popen([command, name]))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update([
                p for p in processes if p.poll() is not None])

print "---- FINSIHED ----"
