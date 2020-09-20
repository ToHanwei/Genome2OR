#!coding:utf-8

import os
import sys
import logging

from src.config import HMMBUILD
from src.Functions import logger
from src.Functions import chmod
from src.ParseArgs import ParseCommand
from src.Functions import platform_info
from src.Functions import sequence_align


def CommandParse():
    """
    Parse command line
    """
    Parse = ParseCommand()
    args = Parse.Iteraparse()
    # profile nhmmer need(hmm, [un]alignment file)
    pofile = args.profile
    # output dir
    outdir = args.outputdir
    # genomic data file
    genome = args.genome
    # number of iteration
    itera = args.iteration
    # number of parallel CPU workers to use
    cpus = str(args.cpus)
    # output file prefix
    prefix = args.prefix
    # Print verbose information
    verbose = args.verbose
    # Sequence similarity threshold
    EvalueLimit = str(args.EvalueLimit)
    # An artificially set OR's sequence length threshold
    SeqLengthLimit = str(args.SeqLengthLimit)
    return pofile, outdir, genome, itera, cpus, prefix, verbose, EvalueLimit, SeqLengthLimit


def CommanRun(pofile, outdir, prefix, genome, cpus, verbose, EvalueLimit, SeqLengthLimit):
    nhmmout = os.path.join(outdir, prefix + ".tblout")
    # run nhmmer program
    hmmcom = ("python nhmmer.py "
              + pofile + " "
              + genome + " "
              + nhmmout
              + " -e " + EvalueLimit
              + " -c " + cpus
              )
    if verbose:
        hmmcom += " -v"
    os.system(hmmcom)

    # run FindOR program
    findcom = ("time python FindOR.py "
               + nhmmout + " "
               + genome
               + " -o " + outdir
               + " -p " + prefix
               + " -e " + EvalueLimit
               + " -l " + SeqLengthLimit
               )
    if verbose:
        findcom += " -v"
    os.system(findcom)

    # run IdentityFunc program
    profile = prefix + '_Pre-ORs_pro.fa'
    dnafile = prefix + "_Pre-ORs_dna.fa"
    profile = os.path.join(outdir, profile)
    dnafile = os.path.join(outdir, dnafile)
    identitycom = ("time python IdentifyFunc.py "
                   + profile + " "
                   + dnafile
                   + " -o " + outdir
                   + " -c " + cpus
                   + " -p " + prefix
                   )
    if verbose:
        identityout += " -v"
    os.system(identitycom)
    return dnafile


def constru_hmm(dnafile):
    prefix = dnafile.split('.')[0]
    alignf = prefix + '.msa'
    hmmfile = prefix + '.hmm'
    sequence_align(dnafile, alignf)
    command = HMMBUILD + ' ' + hmmfile + ' ' + alignf
    os.system(command)
    return hmmfile


def main():

    Iterdoc = """
    ============================================
                                              
                  Iteration {}                
                                              
    ============================================
    """

    # correct input error
    if len(sys.argv) == 1:
        os.system("python Iteration.py -h")
        sys.exit(1)
    # open logging
    logger('Iteration')
    # change tool permission
    chmod(HMMBUILD)
    # get comandlines
    pofile, outdir, genome, itera, cpus, prefix, verbose, EvalueLimit, SeqLengthLimit = CommandParse()
    # check platform
    platform_info(verbose)
    # iteration runing
    hmmfile = pofile
    for _iter in range(1, itera+1):
        print(Iterdoc.format(_iter))
        prefix_i = prefix + '_itera' + str(_iter)
        dnafile = CommanRun(
            hmmfile,
            outdir,
            prefix_i,
            genome,
            cpus,
            verbose,
            EvalueLimit,
            SeqLengthLimit
        )
        if _iter < itera:
            logging.info('Iteration {} construct HMM profile'.format(_iter))
            hmmfile = constru_hmm(dnafile)
    logging.info("###Program Iteration.py finish###")


if __name__ == "__main__":
    main()


