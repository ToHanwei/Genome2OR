#!coding:utf-8

import os
import sys
import logging

from src.config import *
from src.Functions import logger
from src.Functions import refact_hitfile
from src.Functions import refact_list
from src.Functions import tm_pattern
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
# hit_file input file
hitfile = args.hitfile
# Result save directory
outputdir = args.outputdir
# output file prefix
prefix = args.prefix
# Print verbose information(or not)
verbose = args.verbose

# open logging
logger('IdentityFunc')

# check platform
platform_info(verbose)


if verbose:
    print("\033[1;32mProcess hit sequence file...\033[0m")
template, hit_dict = refact_hitfile(hitfile)

if verbose:
    print("\033[1;32mRecombination hit list...\033[0m")
hit_list = refact_list(template, hit_dict)


if verbose:
    print("\033[1;32mIdentity functions\\pseudogenes...\033[0m")
pseus, funcs = [], []
for seq_list in hit_list:
    # Some pseudogenes were filtered out by pattern matching
    pseu1, nterms, tm_list = tm_pattern(seq_list)
    # Some pseudogenos were filtered out by N-term length
    func1, pseu2 = Nterm_length(nterms)
    tm_list = [(n, seq) for n, seq in tm_list if n in func1]
    # Some pseudogenoes were filtered out by TM gaps
    func2, pseu3 = tm_gaps_filter(tm_list)
    seq_dict = dict(seq_list)
    for func in func2:
        seq_func = ">" + func + seq_dict[func].replace('-', '')
        funcs.append(seq_func)
    for pseu in pseu1 + pseu2 + pseu3:
        seq_pseu = ">" + pseu + seq_dict[pseu].replace('-', '')
        pseus.append(seq_pseu)

if verbose:
    print("\033[1;32mWrite data to file...\033[0m")
func_file = os.path.join(outputdir, prefix+"_redundant_func_ORs.fasta")
pseu_file = os.path.join(outputdir, prefix+"_redundant_pseu_ORs.fasta")
with open(func_file, 'w') as funcf:
    funcf.writelines(funcs)
with open(pseu_file, 'w') as pseuf:
    pseuf.writelines(pseus)

if verbose:
    print("\033[1;32mCD-HIT filter redundant sequence...\033[0m")
nonredun_func = os.path.join(outputdir, prefix+"_func_ORs.fasta")
funcomm = CDHIT + " -i " + func_file + " -c 1.0 -o " + nonredun_func
os.system(funcomm)
nonredun_pseu = os.path.join(outputdir, prefix+"_pseu_ORs.fasta")
pseucomm = CDHIT + " -i " + pseu_file + " -c 1.0 -o " + nonredun_pseu
os.system(pseucomm)

logging.info("###Program IdentityFunc.py finish###")