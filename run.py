#!/usr/bin/env python3

import sys
import os
import time
import subprocess
import argparse

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

def process_directories(target):
    for d in target:
        if os.path.isdir(d):
            for root, dirs, filenames in os.walk(d):
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
        print(results, end="", flush=True)

def main():

    parser = argparse.ArgumentParser(description='OCL kernel compiler evaluation program')
    parser.add_argument("--platform", nargs='?', default='0',
                        metavar='0 | 1',
                        help='Specify target OCL platform: 0 means beignet, 1 means VPG. Beignet is default.')
    parser.add_argument("kernel", nargs='*', default=['kernel/'],
                        metavar='kernel directory',
                        help='Specify OCL kernel run target directory.')
    args = parser.parse_args()
    
    try:
        os.stat("run")
    except OSError:
        print("First to build run by make")
        sys.exit(1)

    process_directories(args.kernel)

    run_cases(args.platform)

if __name__ == "__main__":
    main()
