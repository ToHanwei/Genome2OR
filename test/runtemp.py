#!coding:utf-8

import os
import sys


indir = sys.argv[1]
hmm = sys.argv[2]
outdir = sys.argv[3]
iters = sys.argv[4]
index = int(sys.argv[5])

files = os.listdir(indir)
for i in range(index, len(files)):
    _file = files[i]
    genome = os.path.join(indir, _file)
    base = os.path.splitext(_file)[0]
    command = ("python Iteration.py"
        + " " + hmm
        + " " + outdir
        + " " + genome
        + " -p " + base
        + " -i " + iters
        + " -c 36")
    os.system(command)
    print(base + "is done")
    print("index: " + str(i))
