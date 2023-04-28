#!coding:utf-8

import os
import sys
import logging

from src.config import PYTHON
from src.config import BATCH_MESSAGE
from src.Functions import logger
from src.Functions import getCpus
from src.Functions import excuteCommand
from src.Functions import existsDirectory
from src.ParseArgs import ParseCommand
from src.Functions import platform_info
from src.main import Nhmmer
from src.main import FindOR
from src.main import Identify


sys.stdout.write(BATCH_MESSAGE)

# correct input error
if len(sys.argv) == 1:
    sys.argv.append("-h")

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
identifyout = args.identifyout
# number of parallel CPU workers to use
cpus = args.cpus
keepfile = str(args.keepfile)
# Print verbose information
verbose = args.verbose
# Sequence similarity threshold
EvalueLimit = str(args.EvalueLimit)
# An artificially set OR's sequence length threshold
SeqLengthLimit = str(args.SeqLengthLimit)

# open logging
logger('batch annotated')

# check platform
platform_info(verbose)
    
# check CPU cores
cpus = getCpus(cpus)
cpus = str(cpus)

# Check if the directory path exists,
# and create the directory if it doesn't exist.
existsDirectory(nhmmerout)
existsDirectory(findoutdir)
existsDirectory(identifyout)

files = os.listdir(indir)
# number of genomic in inputdir Folder
nums = len(files)
inpaths = [os.path.join(indir, fi) for fi in files]
npaths = [os.path.join(nhmmerout, fi + ".tblout") for fi in files]

# batch run nhmmer and FindOR program
for infile in files:
    logging.info(
        "\033[0;32mAnnotating genome {}\033[0m".format(infile)
    )
    prefix = os.path.splitext(infile)[0]
    genome = os.path.join(indir, infile)
    nhmmout = os.path.join(nhmmerout, prefix + ".tblout")
    # run nhmmer program
    nhmmer = Nhmmer(
        output = nhmmout,
        genome = genome,
        profile = pofile,
        cpus = cpus,
        verbose = verbose,
        EvalueLimit = EvalueLimit,
    )
    nhmmer.run()

    # run FindOR program
    findor = FindOR(
        input = nhmmout,
        genome = genome,
        verbose = verbose,
        outputdir= findoutdir,
        prefix = prefix,
        EvalueLimit = EvalueLimit,
        SeqLengthLimit = SeqLengthLimit,
    )
    findor.run()

    # run IdentifyFunc program
    profile = prefix + '_Pre-ORs_pro.fa'
    dnafile = prefix + "_Pre-ORs_dna.fa"
    profile = os.path.join(findoutdir, profile)
    dnafile = os.path.join(findoutdir, dnafile)
    identify = Identify(
        hitPROfile = profile,
        hitDNAfile = dnafile,
        cpus = cpus,
        outputdir = identifyout,
        prefix = prefix,
        keepfile = keepfile,
        verbose = verbose,
    )
    identify.run()

logging.info("###The batch.py has completed.###\n")



