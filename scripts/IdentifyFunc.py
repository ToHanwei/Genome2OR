#!coding:utf-8

import os
import sys
import logging

from src.config import *
from src.Functions import logger
from src.Functions import chmod
from src.Functions import refact_hitfile
from src.Functions import identify_filter
from src.Functions import identity_writer
from src.Functions import refact_list
from src.Functions import tm_prepare
from src.Functions import unredundant
from src.Functions import ReadSampleFasta
from src.Functions import platform_info
from src.Functions import Nterm_length
from src.Functions import tm_gaps_filter
from src.ParseArgs import ParseCommand


if len(sys.argv) == 1:
    os.system('python IdentifyFunc.py -h')
    sys.exit(1)

# Parse command line
Parse = ParseCommand()
args = Parse.identifyparse()
# hit protein file input file
hitpro = args.hitPROfile
# hit DNA file input file
hitdna = args.hitDNAfile
# number of parallel 
cpus = args.cpus
# Result save directory
outputdir = args.outputdir
# output file prefix
prefix = args.prefix
keepfile = args.keepfile
# Print verbose information(or not)
verbose = args.verbose

# open logging
logger('IdentityFunc')

# check platform
platform_info(verbose)

# change nhmmer permission
chmod(MAFFT)
chmod(CDHIT)

if verbose:
    print("\033[1;32mProcess hit sequence file...\033[0m")
template, hit_dict = refact_hitfile(hitpro)

if verbose:
    print("\033[1;32mRecombination hit list...\033[0m")
hit_list = refact_list(template, hit_dict, cpus)

if verbose:
    print("\033[1;32mIdentity functions\\pseudogenes...\033[0m")
pseus, funcs = [], []
num_of_paeu = {"NTERM": 0, "GAP": 0}

pre_dnas_dict = dict(ReadSampleFasta(hitdna))

for seq_list in hit_list:
    assert len(seq_list) > 1, print("Input file is FindOR.py output file. Check place!")
    # Prepare tm_list
    nterms, tm_list = tm_prepare(seq_list)
    seq_dict = dict(seq_list)

    # identity function ORs
    func, pseu, ptype = identify_filter(seq_dict, tm_list, nterms)
    if func:
        func_seq = ">" + func + seq_dict[func].replace('-', '')
        funcs.append(func_seq)
    elif pseu:
        pname = pseu.strip() + '_' + ptype + '\n'
        pseu_seq = ">" + pname + pre_dnas_dict[pseu].replace('-', '')
        pseus.append(pseu_seq)
        num_of_paeu[ptype] += 1
    else:
        print("Function either or pseudogene")

if verbose:
    print("\033[1;32mWrite data to file...\033[0m")
files = identity_writer(
    hitpro, outputdir,
    prefix, funcs,
    pseus, num_of_paeu
    )
func_file, pseu_file, prepseu = files

# Unredundance sequences and write DNA sequences to file.
if verbose:
    print("\033[1;32mCD-HIT filter redundant sequence...\033[0m")
unredundant(outputdir, prefix, func_file, pseu_file, hitdna)

# delect file, if keepfile=False
if not keepfile:
   os.remove(hitpro)
   os.remove(hitdna)
   os.remove(func_file)
   os.remove(pseu_file)
   os.remove(prepseu)
   os.system("rm *clstr")

logging.info("###Program IdentityFunc.py finish###")
