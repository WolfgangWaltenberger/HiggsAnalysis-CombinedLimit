#!/usr/bin/env python
        
import subprocess
import numpy

factors=numpy.arange(0.25,3.76,.25).tolist()
factors.append(0.20)
factors.remove(1.00)

dry_run=True

for signal in [ "-1", "-2", "-12" ]:
    for factor in factors:
        f=open("crab.cfg.template")
        lines=f.readlines()
        f.close()
        g=open("crab.cfg","w")
        for line in lines:
            line=line.replace("@@SIGNAL@@",signal)
            line=line.replace("@@FACTOR@@",str(factor) )
            g.write ( line )
        g.close()
        if dry_run:
            continue
        o=subprocess.check_output( "crab -create -submit", shell=True )
        lines=o.split("\n")
        for line in lines:
            print("[crab_submit.py] %s" % line)

