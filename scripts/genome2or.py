#!coding:utf-8

import os
import sys
import logging

from src.config import PYTHON
from src.config import NHMMER
from src.config import CDHIT
from src.config import MAFFT
from src.config import HMMBUILD
from src.config import GENOME2OR_MESSAGE
from src.ParseArgs import ParseCommand
from src.Functions import logger
from src.Functions import chmod
from src.Functions import getCpus
from src.Functions import platform_info
from src.Functions import sequence_align
from src.Functions import excuteCommand
from src.Functions import existsDirectory
from src.main import Nhmmer
from src.main import FindOR
from src.main import Identify


def CommandParse():
    """
    Parse command line
    """
    Parse = ParseCommand()
    args = Parse.genome2or()
    # profile nhmmer need(hmm, [un]alignment file)
    pofile = args.profile
    # output dir
    outdir = args.outputdir
    # genomic data file
    genome = args.genome
    # number of parallel CPU workers to use
    cpus = args.cpus
    keepfile = str(args.keepfile)
    # output file prefix
    prefix = args.prefix
    # Print verbose information
    verbose = args.verbose
    # Sequence similarity threshold
    EvalueLimit = str(args.EvalueLimit)
    # An artificially set OR's sequence length threshold
    SeqLengthLimit = str(args.SeqLengthLimit)
    return (
        pofile, outdir, genome,
        cpus, prefix, verbose, 
        EvalueLimit, keepfile,
        SeqLengthLimit
    )


def CommanRun(
        pofile, outdir, prefix,
        genome, cpus, verbose,
        EvalueLimit, SeqLengthLimit,
        keepfile):
    nhmmout = os.path.join(outdir, prefix + ".tblout")
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
        outputdir= outdir,
        prefix = prefix,
        EvalueLimit = EvalueLimit,
        SeqLengthLimit = SeqLengthLimit,
    )
    findor.run()

    # run IdentifyFunc program
    profile = prefix + '_Pre-ORs_pro.fa'
    dnafile = prefix + "_Pre-ORs_dna.fa"
    profile = os.path.join(outdir, profile)
    dnafile = os.path.join(outdir, dnafile)
    identify = Identify(
        hitPROfile = profile,
        hitDNAfile = dnafile,
        cpus = cpus,
        outputdir = outdir,
        prefix = prefix,
        keepfile = keepfile,
        verbose = verbose,
    )
    identify.run()
    return dnafile


def main():

    sys.stdout.write(GENOME2OR_MESSAGE)
    
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    
    # get comandlines
    (
        pofile, outdir, genome,
        cpus, prefix, verbose, 
        EvalueLimit, keepfile,
        SeqLengthLimit
    ) = CommandParse()

    # open logging
    logger('genome2or')

    # change tool permission
    chmod(HMMBUILD)
    chmod(NHMMER)
    chmod(MAFFT)
    chmod(CDHIT)

    # check platform
    platform_info(verbose)

    # check CPU cores
    cpus = getCpus(cpus)
    cpus = str(cpus)

    # Check if the directory path exists,
    # and create the directory if it doesn't exist.
    existsDirectory(outdir)

    # main runing
    dnafile = CommanRun(
        pofile,
        outdir,
        prefix,
        genome,
        cpus,
        verbose,
        EvalueLimit,
        SeqLengthLimit,
        keepfile
    )
    logging.info("###The Iteration.py has completed.###\n")


if __name__ == "__main__":
    main()


