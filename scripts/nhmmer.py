#!/usr/local/bin/python
#coding:utf-8

import os
import sys
import logging

from src.Functions import logger
from src.Functions import chmod
from src.Functions import platform_info
from src.Functions import run_nhmmer
from src.ParseArgs import ParseCommand
from src.config import NHMMER


if len(sys.argv) == 1:
    os.system("python nhmmer.py -h")
    sys.exit(1)

# Parse command line
Parse = ParseCommand()
args = Parse.nhmmerparse()
# output file path
output = args.output
# Genomic data file path
gefile = args.genome
# profile nhmmer need(hmm, [un]alignment file)
pofile = args.profile
# number of parallel CPU workers to use
cpus = args.cpus
# Print verbose information
verbose = args.verbose
# Sequence similarity threshold
EvalueLimit = args.EvalueLimit

# open logging
logger('nhmmer')

# check platform
platform_info(verbose)

# change nhmmer permission
chmod(NHMMER)

if os.path.isdir(output):
    output = os.path.join(output, 'nhmmer_out.tblout')
    print("You input is folder, write data to {}".format(output))

# Run the nhmmer program
if verbose:
    print("\n\033[1;32mNHMMER program is running... \033[0m\n")
code = run_nhmmer(EvalueLimit, cpus, output, pofile, gefile, verbose)
if (code == 0) and verbose:
    print("\n\033[1;32mNHMMER program is complete.\033[0m\n")
elif verbose:
    print("\n\033[1;31mNHMMER program ran incorrectly.\033[0m\n")
logging.info("###Program nhmmer.py finish###")

