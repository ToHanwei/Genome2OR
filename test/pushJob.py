#!coding:utf-8

import re
import os
import sys
import time

sys.setrecursionlimit(3000)

def getfiles(indir):
    """
    get input files
    genaration input list
    """
    outlist = []
    outdir = os.path.basename(indir)
    logdir = outdir + '_outlog/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    files = os.listdir(indir)
    for _file in files:
        infile = os.path.join(indir, _file)
        log = logdir + _file
        outlist.append((log, outdir, infile))
    return outlist


def checkjob():
    """
    check bhosts no jobs node
    """
    while True:
        hosts = os.popen('bhosts').read()
        hosts = hosts.strip().split('\n')
        hosts = [line.split() for line in hosts[1:-1]]
        hosts = [(ele[0], ele[4]) for ele in hosts]
        empty = [ele for ele in hosts if ele[1] =='0']
        if empty:
            break
        else:
            time.sleep(10)
            empty = checkjob()
    return empty

def generate():
    indir = sys.argv[1]
    hmm = sys.argv[2]
    iters = sys.argv[3]
    joblist = getfiles(indir)
    empty = None
    have_done = []
    for i, job in enumerate(joblist):
        if not empty:
            empty = checkjob()
            time.sleep(10)
        node = empty.pop()[0]
        log = job[0]
        base = os.path.split(log)[1]
        base = os.path.splitext(base)[0]
        outdir = job[1]
        infile = job[2]
        jobfile = []
        fout = open('jobfile.lsf', 'w')
        jobfile.append("#!/bin/bash\n")
        jobfile.append("#BSUB -J "+log+"\n")
        jobfile.append("#BSUB -e "+log+".err\n")
        jobfile.append("#BSUB -o "+log+".out\n")
        jobfile.append("#BSUB -q zhaolab\n")
        jobfile.append("#BSUB -m "+node+"\n")
        jobfile.append("#BSUB -n 36\n\n\n")
        jobfile.append(
            "python Iteration.py "
            + hmm + " "
            + outdir + " "
            + infile
            + " -p " + base
            + " -i " + iters
            + " -c 36"
            + "\n")
        fout.writelines(jobfile)
        fout.close()
        print(log)
        os.system('bsub < jobfile.lsf')
        print('COUNT: ', i)


if __name__ == "__main__":
    generate()

