#coding:utf-8

import re
from collections import Counter
import matplotlib.pyplot as plt

def ReadSampleFasta(seqfile):
    """
    Function:
        ReadSampleFasta
        read small FASTA format file.
        Note,large files consume a lot of memory.
    Parameter:
        seqfile, input file name, FASTA format file
    Return:
        return type -> List
        seq_list, [(name1, seq1), (name2, seq2), ...]
    """

    seq_list = []
    with open(seqfile) as seqf:
        text = seqf.read().replace('\r', '')
        seqs = text.split('>')[1:]
        for seq in seqs:
            lines = seq.split('\n')
            name = lines[0] + "\n"
            aads = ''.join(lines[1:]) + "\n"
            aads = aads.replace('*', '')
            seq_list.append((name, aads))
    return seq_list


def tm_cut(seq_list):
    """
    Function:
        tm_cut
    Parameter:
        seq_list, sequence alignment list
    Return:
        return type -> dict
        nterm_dict, key   -> sequence name
                    value ->  N terminal sequence
        return type -> dict
        cut_dict,  key   -> sequence name
                   value -> TM sequence
    """
    iaad = 0
    tm1, tm2, tm3, tm4, ecl2, tm5, tm6, tm7 = [],[],[],[],[],[],[],[]
    for i in range(len(seq_list[0][1])):
        aad = seq_list[0][1][i]
        if aad != "-":
            iaad += 1
        else:
            continue
        if iaad == 24: tm1.append(i)
        if iaad == 53: tm1.append(i+1)
        if iaad == 60: tm2.append(i)
        if iaad == 87: tm2.append(i+1)
        if iaad == 94: tm3.append(i)
        if iaad == 131: tm3.append(i+1)
        if iaad == 139: tm4.append(i)
        if iaad == 164: tm4.append(i+1)
        if iaad == 197: tm5.append(i)
        if iaad == 229: tm5.append(i+1)
        if iaad == 235: tm6.append(i)
        if iaad == 263: tm6.append(i+1)
        if iaad == 268: tm7.append(i)
        if iaad == 293: tm7.append(i+1)
    tm_list, nterm, h8 = [], [], []
    ecl1, ecl2, ecl3 = [], [], []
    icl1, icl2, icl3 = [], [], []
    for name, seq in seq_list:
        # append n-term sequence to list
        nterm.append(seq[:tm1[0]])
        h8.append(seq[tm7[1]:])
        tms = []
        for tm in [tm1, tm2, tm3, tm4, tm5, tm6, tm7]:
            tm_seq = seq[tm[0]:tm[1]]
            tms.append(tm_seq)
        tm_list.append(tms)
        icl1.append(seq[tm1[1]:tm2[0]])
        icl2.append(seq[tm3[1]:tm4[0]])
        icl3.append(seq[tm5[1]:tm6[0]])
        ecl1.append(seq[tm2[1]:tm3[0]])
        ecl2.append(seq[tm4[1]:tm5[0]])
        ecl3.append(seq[tm6[1]:tm7[0]])
    icl = [icl1, icl2, icl3]
    ecl = [ecl1, ecl2, ecl3]
    return nterm, icl, ecl, tm_list, h8


def tm_pattern_statistic(tm_list, pattern):

    leng = len(tm_list)
    tm_list = (tm.replace('-', '') for tm in tm_list)
    n = sum(1 for tm in tm_list if re.search(pattern, tm))
    return round((n*100.0)/leng, 2)


def Ploter(match):
    plt.plot(match)
    plt.plot(match, 'o')
    plt.hlines(90, -0.5, len(match)-0.5, colors='r', linestyles='--')
    plt.title("OR Pattern matching percentage")
    plt.ylabel("Percentage(%)")
    ticks = ('N-term', "TM1", "ICL1", "TM2(1)", "TM2(2)",
             "ICL2", "TM3(1)", "TM3(2)", "ECL2",
             "TM5", "TM6", "TM7", "H8"
             )
    plt.xticks(range(len(match)), ticks, rotation=45)
    plt.ylim(ymin=70, ymax=100)
    plt.savefig('OR_pattern_match.png', format='png')
    plt.show()
    with open('OR_pattern_match.tab', 'w') as outf:
        lines = ['pattern\tpercentage\n']
        d = [x+'\t'+str(y)+'\n' for x,y in zip(ticks, match)]
        lines += d
        outf.writelines(lines)


def pattern_distribution(patterns, seqs):
    match_list = []
    for seq in zip(*seqs):
        seq = [tm.replace('-', '') for tm in seq]
        m = sum(1 for pat, tm in zip(patterns, seq) if re.search(pat, tm))
        match_list.append(m)

    count = Counter(match_list)
    total = sum(count.values())
    count = sorted(count.items(), key=lambda x: x[1], reverse=True)
    X, Y = list(zip(*count))
    Y_hat = [y*1.0/total for y in Y]

    plt.bar(X, Y_hat, linewidth=1.5, width=0.9)
    plt.title('ORs Pattern Match Distribution')
    plt.xlabel("# of pattern match")
    plt.ylabel("Percentage(%)")
    plt.xticks(range(len(X)), X[::-1])
    plt.savefig("ORs_pattern_match_distribution.png", format='png')
    plt.show()
    with open("ORs_pattern_match_distribution.tab", 'w') as outf:
        lines = ['match_number\tmatch_percentage\tpattern_match\n']
        d = [str(x)+'\t'+str(y)+'\t'+str(z)+'\n' for x,y,z in zip(Y, Y_hat, X)]
        lines += d
        outf.writelines(lines)


TM1 = r'[YF].....[GAESW]N'
TM2_1 = r'[MK][YF].[FL][LIV]'
TM2_2 = r'[LF]...[DEN]........P'
#TM3_1 = r'C..Q'
# TM3_2 = r'[LFY]..M..[DN][RH]'
TM3_1 = r'C..Q..............[LFYI]..[ML]..[DN][RHCQLW]'
TM3_2 = r'A[IV]..PL'
TM5 = r'[ST]Y..[IVL]'
TM6 = r'[KR]...T[CL]..H'
# TM7 = r'[PSA]..NP..[YF]'
TM7 = r'[NSY]P.[IVL][YF]'
# Nterm = r'[FL].[LFI].[GA]'
# Nterm = r'[FLIV].[LFIM].[GA]'
# Nterm = r'[FLIV].[LFIM].[GASECPR][LFIVM]'
Nterm = r'[FLIV].[LFIM].[GAS]'
# H8 = r'[RK][NTS][KRQ][EDQK].[KRQ]'
H8 = r'[RKQ]...[VIMLF]...[LMVIFA]'
ICL1 = r'L..P'
ICL2 = r'Y...[MLVI]'
# ECL2 = r'[CYF].........C.........[CSY]'
ECL2 = r'[CYFRHS].........C.........[CSYA]'


infile = "HuOR_E60SecondE60Res_OR_sequence-lengfilter280.fas"

seq_list = ReadSampleFasta(infile)
nterm, icl, ecl, tms, h8 = tm_cut(seq_list)
tms = list(zip(*tms))

# --------------Start Pattern match percentage--------------

m_tm1 = tm_pattern_statistic(tms[0], TM1)
m_tm2_1 = tm_pattern_statistic(tms[1], TM2_1)
m_tm2_2 = tm_pattern_statistic(tms[1], TM2_2)
m_tm3_1 = tm_pattern_statistic(tms[2], TM3_1)
m_tm3_2 = tm_pattern_statistic(tms[2], TM3_2)
m_tm5 = tm_pattern_statistic(tms[4], TM5)
m_tm6 = tm_pattern_statistic(tms[5], TM6)
m_tm7 = tm_pattern_statistic(tms[6], TM7)
m_nterm = tm_pattern_statistic(nterm, Nterm)
m_h8 = tm_pattern_statistic(h8, H8)
m_icl1 = tm_pattern_statistic(icl[0], ICL1)
m_icl2 = tm_pattern_statistic(icl[1], ICL2)
m_ecl2 = tm_pattern_statistic(ecl[1], ECL2)

match = [
    m_nterm, m_tm1, m_icl1,
    m_tm2_1, m_tm2_2, m_icl2,
    m_tm3_1, m_tm3_2,
    m_ecl2, m_tm5, m_tm6,
    m_tm7, m_h8
]

Ploter(match)

# --------------End Pattern match percentage--------------

# ----------Start ORs Pattern match distribution----------
patterns = [
    TM1, TM2_1, TM2_2,
    TM3_1, TM3_2, TM5,
    TM6, TM7, ICL1,
    ICL2, ECL2, Nterm, H8
    ]
seqs = [tms[0], tms[1], tms[1],
        tms[2], tms[2], tms[4],
        tms[5], tms[6], icl[0],
        icl[1], ecl[1], nterm, h8,
        ]

pattern_distribution(patterns, seqs)
# ----------End ORs Pattern match distribution----------
