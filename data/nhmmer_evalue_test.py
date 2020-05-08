#!coding:utf-8

import os
import sys
import logging

from src.Functions import logger
from src.ParseArgs import ParseCommand
from src.Functions import platform_info
from src.Functions import test_nhmmer_evalue

"""
This script must move to script folder, and runing.
"""

if len(sys.argv) == 1:
    os.system("python nhmmer_evalue_test.py -h")
    sys.exit(1)

# Parse command line
Parse = ParseCommand()
args = Parse.testparse()
# profile nhmmer need(hmm, [un]alignment file)
pofile = args.profile
# input file path
indir = args.inputdir
# run nhmmer program result
nhmmerout = args.nhmmerout
# number of parallel CPU workers to use
cpus = str(args.cpus)
# Print verbose information
verbose = args.verbose
# Sequence similarity threshold
EvalueLimit = str(args.EvalueLimit)
# test out file path
outfile = args.output

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
outlines  = []
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
    
    name, nhits = test_nhmmer_evalue(nhmmout)
    outlines.append(name + '\t' + str(nhits) + '\n')

with open(outfile, 'w') as outf:
    outf.writelines(outlines)

logging.info("###Program batch.py finish###")
