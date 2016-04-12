#!/usr/bin/env python3

import sys
import os
import time
import subprocess

# case dict
# e.g {
#        "kernel/foo/foo.cl" : "",
#        "kernel/xxx.cl" : "kernel/xxx.opt"
#        ...
# }
cases = {}

# blacklist case
blacklist = [
    "kernel/luxmark/hotel-1",   #beignet takes too long time > 600s
    ]

def process_directories():
    for root, dirs, filenames in os.walk("kernel/"):
        filenames.sort()
        for f in filenames:
            item = os.path.splitext(f)
            if root + "/" + item[0] in blacklist:
                continue
            if item[1] == ".cl":
                cases[root+"/"+f] = ""
            elif item[1] == ".opt":
                cases[root+"/"+item[0]+".cl"] = root + "/" + f

def run_cases(drv):
    for k in cases:
        cmd = ['./run', '-p', drv, k]
        if cases[k] != '':
            cmd.append(cases[k])
        try:
            p = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            (stdout, stderr) = p.communicate()
            results = stdout.decode("UTF-8")
        except KeyboardInterrupt:
            sys.exit(1)
        print(results)

def main():

    try:
        os.stat("run")
    except OSError:
        print("First to build run by make")
        sys.exit(1)

    process_directories()

    run_cases("0")

if __name__ == "__main__":
    main()
