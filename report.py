#!/usr/bin/env python3

import sys
import os
import csv

if len(sys.argv) != 3:
    print(sys.argv[0] + " <old> <new>")
    sys.exit(1)

old_f = open(sys.argv[1], newline='')
new_f = open(sys.argv[2], newline='')

o_list=[]
old_csv = csv.reader(old_f)
for o in old_csv:
    o_list.append(o)
o_list.sort()

n_list=[]
new_csv = csv.reader(new_f)
for n in new_csv:
    n_list.append(n)
n_list.sort()

print(o_list[0][1] + ": binary size | build time\n")

for o in o_list:
    for n in n_list:
        if o[0] != n[0]:
            continue
        if o[1] != n[1]:
            print("Non-equal platform")
            sys.exit(1)
        result=o[0] + ": "
        size_diff = float(n[2]) - float(o[2])
        size_c = abs(size_diff)/float(o[2])
        if size_diff < 0:
            result+="-"+"{:.2f}%".format(size_c)
        elif size_diff == 0:
            result+="="
        else:
            result+="+"+"{:.2f}%".format(size_c)
        result+=" | "
        time_diff = float(n[3]) - float(o[3])
        time_c = abs(time_diff)/float(o[3])
        if time_diff < 0:
            result+="-"+"{:.2f}%".format(time_c)
        elif time_diff == 0:
            result+="="
        else:
            result+="+"+"{:.2f}%".format(time_c)
        print(result)
            
