#!/usr/bin/env python

from glob import glob
import subprocess
import os
import time

dataset = "tt012j_bbars_2l_FxFx"

torun = glob("/xrootd/store/user/iawatson/%s/GEN/*.root" % dataset)
max_processes = 15
run_on_batch = False
delete_hadtruth = False
delete_had = False
delete_nano = True
delete_log = True

command = "./rerun.sh"
processes = set()

print "---- RUNNING ----"

if delete_had:
    os.system("rm -f /xrootd/store/user/iawatson/%s/HADAOD/*.root" % dataset)

if delete_hadtruth:
    os.system("rm -f /xrootd/store/user/iawatson/%s/HADTRUTHAOD/*.root" % dataset)

if delete_nano:
    os.system("rm -f /xrootd/store/user/iawatson/%s/NANOAOD/*.root" % dataset)

if delete_log:
    os.system("rm -f /xrootd/store/user/iawatson/%s/LOG/*" % dataset)

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
