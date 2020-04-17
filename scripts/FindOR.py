#!coding:utf-8

import os
import sys
import logging

from src.Functions import logger
from src.Functions import find_cds
from src.Functions import extract_cds
from src.Functions import writer2file
from src.Functions import platform_info
from src.Functions import proc_nhmmer_out
from src.ParseArgs import ParseCommand


if len(sys.argv) == 1:
    os.system("python FindOR.py -h")
    sys.exit(1)

# Parse command line
Parse = ParseCommand()
args = Parse.parsecommand()
# nhmmer output file path, as input
infile = args.input
# Genomic data file path
gefile = args.genome
# Print verbose information
verbose = args.verbose
# Result save directory
outputdir = args.outputdir
# output file prefix
prefix = args.prefix
# Sequence similarity threshold
EvalueLimit = args.EvalueLimit
# An artificially set OR's sequence length threshold
SeqLengthLimit = args.SeqLengthLimit


# open logging
logger('FindOR')

# check platform
platform_info(verbose)

if verbose:
    print("\033[1;32mProcess nhmmer output file...\033[0m")
hmmout, hitnames = proc_nhmmer_out(infile, EvalueLimit)

if verbose:
    print("\033[1;32mExtract cds from genomic file...\033[0m")
hmmout_seq = extract_cds(hmmout, gefile)

if verbose:
    print("\033[1;32mFind ATG and STOP codons for each sequence...\033[0m")
functional, outliers = find_cds(hmmout, hmmout_seq, SeqLengthLimit)

if verbose:
    print("\033[1;32mWrite data to file...\033[0m")
writer2file(outputdir, prefix, functional, outliers, hmmout_seq, hmmout)

logging.info("###Program FindOR.py finish###")
