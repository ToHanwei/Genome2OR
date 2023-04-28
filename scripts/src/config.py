#!coding:utf-8

import os
import sys

from shutil import which
from os.path import dirname
from os.path import abspath

# The name of this tool
TOOLNAME = "Genome2OR"

# The version
VERSION = "Genome2OR 1.2.3"

# Complementary base
CompBase = {
    'A': 'T', 'T': 'A',
    'C': 'G', 'G': 'C',
    'N': 'N',
    }

# TM boundary, refer to OR5AN1(UniprotKB:Q8NGI8)
TM_boundary = [
    24, 53, 60, 87, 94,
    131, 139, 164, 197,
    229, 235, 263, 268,
    293
]


# Olfactory receptor protein pattern
TM1 = r'[YF].....[GAESW]N'
TM2_1 = r'[MK][YF].[FL][LIV]'
TM2_2 = r'[LF]...[DEN]........P'
TM3_1 = r'C..Q..............[LFYI]..[ML]..[DN][RHCQLW]'
TM3_2 = r'A[IV]..PL'
TM5 = r'[ST]Y..[IVL]'
TM6 = r'[KR]...T[CL]..H'
TM7 = r'[NSY]P.[IVL][YF]'
Nterm = r'[FLIV].[LFIM].[GAS]'
H8 = r'[RKQ]...[VIMLF]...[LMVIFA]'
ICL1 = r'L..P'
ICL2 = r'Y...[MLVI]'
ECL2 = r'[CYFRHS].........C.........[CSYA]'

PATTERNS = [
    TM1, TM2_1, TM2_2,
    TM3_1, TM3_2, TM5,
    TM6, TM7, ICL1,
    ICL2, ECL2, Nterm, H8
]

# N-terminal length interval partition
C1 = 12
C2 = 48
B1 = 16
B2 = 36
A1 = 19
A2 = 26
PEAK = 22

DIST_PEAK = 23

# The region of transmembrane helix parameter
TM_GAPS_TOTAL = 5

# The extend length.
# If CDS shorted than EXTEND_LENGTH will be extended to that length
EXTEND_LENGTH = 1200

# the limit of pseudogene and low-quality sequences
LENGTHLIMIT = 100

# Codon mapping table
CODON_TABLE = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R',
    'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S',
    'AGC': 'S', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'TTA': 'L',
    'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GTT': 'V',
    'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'AAT': 'N', 'AAC': 'N', 'GAT': 'D', 'GAC': 'D',
    'TGT': 'C', 'TGC': 'C', 'CAA': 'Q', 'CAG': 'Q', 'GAA': 'E',
    'GAG': 'E', 'CAT': 'H', 'CAC': 'H', 'AAA': 'K', 'AAG': 'K',
    'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y', 'ATG': 'M',
    'TGG': 'W', 'TAG': '*', 'TGA': '*', 'TAA': '*'
}



class DependToolError(Exception):
    def __init__(self, tool):
        self.tool = "Dependency tool {} not found.".format(tool)

    def __str__(self):
        return self.tool


def getAbsPath():
    """
    Get absolute path to resource.
    """
    #return abs_path
    absPaht = ""
    if getattr(sys, 'frozen', False):
        absPath = dirname(abspath(sys.executable))
    elif __file__:
        absPath = dirname(abspath(__file__))
    
    SCRIPT_PATH = dirname(absPath)
    TEMPLATE_PATH = os.path.join(dirname(SCRIPT_PATH), "template")
    return SCRIPT_PATH, TEMPLATE_PATH

SCRIPT_PATH, TEMPLATE_PATH = getAbsPath()


# define nhmmer tool path 
PYTHON = which("python")
if not PYTHON:
    raise DependToolError("python")

NHMMER = which("nhmmer")
if not NHMMER:
    raise DependToolError("nhmmer")

# define mafft tool path
MAFFT = which("mafft")
if not MAFFT:
    raise DependToolError("mafft")

# define cd-hit tool path
CDHIT = which("cd-hit")
if not CDHIT:
    raise DependToolError("cd-hit")

# define hmmbuild tool path
HMMBUILD = which('hmmbuild')
if not HMMBUILD:
    raise DependToolError("hmmbuild")



GENOME2OR_MESSAGE = """\033[0;32m
    ____ _____ _   _  ___  __  __ _____ ____   ___  ____  
   / ___| ____| \ | |/ _ \|  \/  | ____|___ \ / _ \|  _ \ 
  | |  _|  _| |  \| | | | | |\/| |  _|   __) | | | | |_) |
  | |_| | |___| |\  | |_| | |  | | |___ / __/| |_| |  _ < 
   \____|_____|_| \_|\___/|_|  |_|_____|_____|\___/|_| \_|
                                                          
\033[0m
"""

BATCH_MESSAGE = """\033[0;32m
                ____    _  _____ ____ _   _ 
               | __ )  / \|_   _/ ___| | | |
               |  _ \ / _ \ | || |   | |_| |
               | |_) / ___ \| || |___|  _  |
               |____/_/   \_\_| \____|_| |_|
                                            
\033[0m
"""

FINDOR_MESSAGE = """\033[0;32m
             _____ ___ _   _ ____   ___  ____  
            |  ___|_ _| \ | |  _ \ / _ \|  _ \ 
            | |_   | ||  \| | | | | | | | |_) |
            |  _|  | || |\  | |_| | |_| |  _ < 
            |_|   |___|_| \_|____/ \___/|_| \_|
                                               
\033[0m
"""

IDENTIFYFUNC_MESSAGE = """\033[0;32m
      ___ ____  _____ _   _ _____ ___ _______   _______ _   _ _   _  ____ 
     |_ _|  _ \| ____| \ | |_   _|_ _|  ___\ \ / /  ___| | | | \ | |/ ___|
      | || | | |  _| |  \| | | |  | || |_   \ V /| |_  | | | |  \| | |    
      | || |_| | |___| |\  | | |  | ||  _|   | | |  _| | |_| | |\  | |___ 
     |___|____/|_____|_| \_| |_| |___|_|     |_| |_|    \___/|_| \_|\____|
                                                                          
\033[0m
"""

ITERATION_MESSAGE = """\033[0;32m
     ___ _____ _____ ____      _  _____ ___ ___  _   _ 
    |_ _|_   _| ____|  _ \    / \|_   _|_ _/ _ \| \ | |
     | |  | | |  _| | |_) |  / _ \ | |  | | | | |  \| |
     | |  | | | |___|  _ <  / ___ \| |  | | |_| | |\  |
    |___| |_| |_____|_| \_\/_/   \_\_| |___\___/|_| \_|
                                                       
\033[0m
"""

ITERBATH_MESSAGE = """\033[0;32m
     ___ _____ _____ ____  ____    _  _____ ____ _   _ 
    |_ _|_   _| ____|  _ \| __ )  / \|_   _/ ___| | | |
     | |  | | |  _| | |_) |  _ \ / _ \ | || |   | |_| |
     | |  | | | |___|  _ <| |_) / ___ \| || |___|  _  |
    |___| |_| |_____|_| \_\____/_/   \_\_| \____|_| |_|
                                                       
\033[0m
"""

NHMMER_MESSAGE = """\033[0;32m
           _   _ _   _ __  __ __  __ _____ ____  
          | \ | | | | |  \/  |  \/  | ____|  _ \ 
          |  \| | |_| | |\/| | |\/| |  _| | |_) |
          | |\  |  _  | |  | | |  | | |___|  _ < 
          |_| \_|_| |_|_|  |_|_|  |_|_____|_| \_|
                                                 
\033[0m
"""


ITERDOC = """\033[0;32m
       =============================================
                     ITERATION {}                  
       =============================================\033[0m
"""
