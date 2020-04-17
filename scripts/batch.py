#!coding:utf-8

import os
import sys
import logging

from src.Functions import logger
from src.ParseArgs import ParseCommand
from src.Functions import platform_info

if len(sys.argv) == 1:
    os.system("python batch.py -h")
    sys.exit(1)

# Parse command line
Parse = ParseCommand()
args = Parse.batchparse()
# profile nhmmer need(hmm, [un]alignment file)
pofile = args.profile
# input file path
indir = args.inputdir
# run nhmmer program result
nhmmerout = args.nhmmerout
# run FindOR program result
findoutdir = args.findoutdir
# run IdentityFunc program result
identityout = args.identityout
# number of parallel CPU workers to use
cpus = str(args.cpus)
# Print verbose information
verbose = args.verbose
# Sequence similarity threshold
EvalueLimit = str(args.EvalueLimit)
# An artificially set OR's sequence length threshold
SeqLengthLimit = str(args.SeqLengthLimit)

# open logging
logger('IdentityFunc')

# check platform
platform_info(verbose)

files = os.listdir(indir)
# number of genomic in inputdir Folder
nums = len(files)
inpaths = [os.path.join(indir, fi) for fi in files]
npaths = [os.path.join(nhmmerout, fi + ".tblout") for fi in files]

# batch run nhmmer and FindOR program
for infile in files:
    genome = os.path.join(indir, infile)
    nhmmout = os.path.join(nhmmerout, infile + ".tblout")
    prefix = infile.split('.')[0]
    # run nhmmer program
    hmmcom = "python nhmmer.py " \
             + pofile + " " \
             + genome + " " \
             + nhmmout \
             + " -e " + EvalueLimit \
             + " -c " + cpus
    if verbose:
        hmmcom += " -v"
    os.system(hmmcom)

    # run FindOR program
    findcom = "time python FindOR.py " \
              + nhmmout + " "\
              + genome \
              + " -o " + findoutdir \
              + " -p " + prefix \
              + " -e " + EvalueLimit \
              + " -l " + SeqLengthLimit
    if verbose:
        findcom += " -v"
    os.system(findcom)

    # run IdentityFunc program
    findout = prefix + '_ORs_pro.fa'
    findout = os.path.join(findoutdir, findout)
    print(findout)
    identitycom = "time python IdentifyFunc.py " \
                  + findout + \
                  " -o " + identityout + \
                  " -p " + prefix
    if verbose:
        identityout += " -v"
    os.system(identitycom)
logging.info("###Program batch.py finish###")
