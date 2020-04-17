#!coding:utf-8

# The name of this tool
TOOLNAME = "Genome2OR"

# The version
VERSION = "Genome2OR 1.0"

# Complementary base
CompBase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
            'N': 'N', }

# Olfactory receptor protein pattern
PATTERNS = [
    r'[YFC].....[GAS]N..[ILMV]',  # TM1
    r'L..PMY..[LI]',  # TM2
    r'[IL]....C..Q',  # TM3-1
    r'M..DR..A',  # TM3-2
    r'A[IV]..PL.Y',  # TM3-3
    r'C.........C.........C',  # ECL2
    r'SY..[IVL]',  # TM5
    r'[KR]...T[CL]..H',  # TM6
    r'P..NP..[YF]'  # TM7
]

TM_boundary = [
    24, 53, 60, 87, 94,
    131, 139, 164, 197,
    229, 235, 263, 268,
    293
]


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

patterns = [
    TM1, TM2_1, TM2_2,
    TM3_1, TM3_2, TM5,
    TM6, TM7, ICL1,
    ICL2, ECL2, Nterm, H8
]

# The standard location of the pattern
GOLD_POS = [34, 54, 91, 117, 124, 168, 216, 235, 282]

# The threshold for matching the number of patterns
PATTERN_THRESHOLD = 10

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

# Olfactory receptor template sequence
ORfile = "../template/template.fasta"
# temporary directory name
TEMP = "../temp"
# temporary OR sequence file name
TEMP_OR = "temp.fasta"
# temporary OR alignment file name
TEMP_ALIGN = TEMP + "/" + "temp.fas"

# nhmmer tool path 
NHMMER = "../tools/nhmmer"
# mafft tool path
MAFFT = "../tools/mafft"
# cd-hit tool path
CDHIT = "../tools/cd-hit"

