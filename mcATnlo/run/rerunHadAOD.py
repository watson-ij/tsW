#!/usr/bin/env python

from glob import glob
import subprocess
import os
import time

dataset = "tt012j_bbars_2l_FxFx_AOD"

torun = glob("/home/iwatson/tsW/mcATnlo/data/%s/MINIAODSIM/*.root" % dataset)
max_processes = 15
run_on_batch = False
delete_had = True
delete_nano = True

command = "./rerun.sh"
processes = set()

print "---- RUNNING ----"

if delete_had:
    os.system("rm -f /home/iwatson/tsW/mcATnlo/data/%s/HADAOD/*.root" % dataset)

if delete_nano:
    os.system("rm -f /home/iwatson/tsW/mcATnlo/data/%s/NANOAOD/*.root" % dataset)

for name in torun:
    n = name.split('/')[-1].split('.')[0]
    if run_on_batch:
        script = 'condor_submit -append arguments="%s %s" rerun.jds' % (dataset, n)
        print script
        os.system(script)
    else:
        # Run locally
        processes.add(subprocess.Popen([command, dataset, n]))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update([
                p for p in processes if p.poll() is not None])

print "---- FINSIHED ----"
