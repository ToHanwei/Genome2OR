#coding:utf-8

from collections import Counter

infile = "HuOR_E60SecondE60Res_OR_sequence-lengfilter280.fas"


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
    icl = (icl1, icl2, icl3)
    ecl = (ecl1, ecl2, ecl3)
    return nterm, icl, ecl, tm_list, h8


def statistic_tm(tm, outfile):
    """
    Function:
        statistic_tm
        The residual evaluation rate of each site in TM region was calculated
    Parameter:
        tm, tm list
        outfile, output file name
    Return: None
    """

    seq_num = len(tm)
    outlist = []
    for cols in zip(*tm):
        # Count site aads
        Count = Counter(cols)
        Count = sorted(Count.items(), key=lambda x:x[1], reverse=True)
        if Count[0][0] == '-': continue
        Count = [(k, round((v*100.0)/seq_num, 1)) for k, v in Count]
        line = '\t'.join([a[0]+'('+str(a[1])+')' for a in Count])
        line += '\n'
        outlist.append(line)
    with open(outfile, 'w') as outf:
        outf.writelines(outlist)


seq_list = ReadSampleFasta(infile)
nterm, icl, ecl, tms, h8 = tm_cut(seq_list)
tms = list(zip(*tms))

statistic_tm(h8, 'statistic_H8.tab')
statistic_tm(nterm, 'statistic_N-term.tab')
statistic_tm(tms[0], 'statistic_tm1.tab')
statistic_tm(tms[1], 'statistic_tm2.tab')
statistic_tm(tms[2], 'statistic_tm3.tab')
statistic_tm(tms[3], 'statistic_tm4.tab')
statistic_tm(tms[4], 'statistic_tm5.tab')
statistic_tm(tms[5], 'statistic_tm6.tab')
statistic_tm(tms[6], 'statistic_tm7.tab')
statistic_tm(icl[0], 'statistic_icl1.tab')
statistic_tm(icl[1], 'statistic_icl2.tab')
statistic_tm(icl[2], 'statistic_icl3.tab')
statistic_tm(ecl[0], 'statistic_ecl1.tab')
statistic_tm(ecl[1], 'statistic_ecl2.tab')
statistic_tm(ecl[2], 'statistic_ecl3.tab')

