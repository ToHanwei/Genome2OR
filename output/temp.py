import os
import sys
from collections import defaultdict

ORfile = sys.argv[1]
infile = sys.argv[2]


def sequence_align(seqf, alignf):
    """
    Function:
        sequence_align
        call MAFFT(linsi) multiple sequence alignment
    Parameter:
        seqf, sequence file name (FASTA format)
        alignf, Multiple sequence alignment results
    Return: None
    """
    command = "mafft" \
              + " --localpair " \
              + "--maxiterate 1000 " \
              + "--quiet " \
              + seqf \
              + " > " \
              + alignf
    os.system(command)


def drop_gap(tm_list):
    """
    Function:
        drop_gap
        delect the common (TM) gap columns
    Parameter:
        tm_list, list of TM
    Return:
        return type -> List
        out_list, clear TM list
    """

    gap_num = len(tm_list)
    gap_col = tuple(['-'] * gap_num)
    tm_cols = list(zip(*tm_list))
    while gap_col in tm_cols:
        tm_cols.remove(gap_col)
    out_list = [''.join(tm) for tm in zip(*tm_cols)]
    return out_list


def Writer(outdir, filename, datalist):
    """
    Function:
        Writer
        write list to file
    Parameter:
        outdir, the directory where the output file is saved
        filename, out file name
        datalist, data list, each element is a row
    Return:
        None
    """
    outfile = os.path.join(outdir, filename)
    with open(outfile, 'w') as outf:
        outf.writelines(datalist)


def drop_gap(tm_list):
    """
    Function:
        drop_gap
        delect the common (TM) gap columns
    Parameter:
        tm_list, list of TM
    Return:
        return type -> List
        out_list, clear TM list
    """

    gap_num = len(tm_list)
    gap_col = tuple(['-'] * gap_num)
    tm_cols = list(zip(*tm_list))
    while gap_col in tm_cols:
        tm_cols.remove(gap_col)
    out_list = [''.join(tm) for tm in zip(*tm_cols)]
    return out_list


with open(ORfile) as ORf:
    template = ORf.readlines()

hit_dict = defaultdict(list)
with open(infile) as hitf:
    line = hitf.readline()
    while line:
        if line[0] == ">":
            header = line[1:]
            headers = header.split('_')
            key = headers[0]
        else:
            line = line.replace('*', '')
            hit_dict[key].append((header, line))
        line = hitf.readline()

for hit, seqs in hit_dict.items():
    templist = []
    for head, seq in seqs:
        head = ">" + head
        templist.extend([head, seq])
    aligns = template + templist
    Writer('../temp', 'temp.fasta', aligns)
    sequence_align('../temp/temp.fasta', '../temp/temp.fas')

    seq_list = []
    with open('../temp/temp.fas') as temf:
        text = temf.read()
        text = text.replace('\r', '')
        seqs = text.split('>')[1:]
        for seq in seqs:
            lines = seq.split('\n')
            name = lines[0]
            aads = "".join(lines[1:])
            seq_list.append((name, aads))
    new_list = seq_list[0:1] + seq_list[4:]

#    iaad, bound = 0, 0
#    for i in range(len(new_list[0][1])):
#        aad = new_list[0][1][i]
#        if aad != "-": iaad += 1
#        if iaad == 24:
#            bound = i
#            break
#    nterms = []
#    for name, seq in new_list[1:]:
#        nterm = seq[:bound].replace('-', '')
#        len_nterm = len(nterm)
#        nterms.append((name, len_nterm))
#    area_a = [(n, l) for n, l in nterms if l <= 18]
#    area_b = [(n, l) for n, l in nterms if 28 > l > 18]
#    area_c = [(n, l) for n, l in nterms if l >= 28]
#    if len(area_a) != 0:
#        print("Area A: ", len(area_a))
#        print(area_a)
#    elif len(area_b) != 0:
#        print("Area B: ", len(area_b))
#        print(area_b)
#    elif len(area_c) !=0:
#        print("Area_C: ", len(area_c))
#        print(area_c)

#    iaad = 0
#    tm1, tm2, tm3, tm4, tm5, tm6, tm7, ecl2 = [],[],[],[],[],[],[],[]
#    for i in range(len(new_list[0][1])):
#        aad = new_list[0][1][i]
#        if aad != "-": iaad += 1
#        if iaad == 24: tm1.append(i)
#        if iaad == 53: tm1.append(i+1)
#        if iaad == 60: tm2.append(i)
#        if iaad == 87: tm2.append(i+1)
#        if iaad == 94: tm3.append(i)
#        if iaad == 131: tm3.append(i+1)
#        if iaad == 139: tm4.append(i)
#        if iaad == 164: tm4.append(i+1)
#        if iaad == 165: ecl2.append(i)
#        if iaad == 196: ecl2.append(i+1)
#        if iaad == 197: tm5.append(i)
#        if iaad == 229: tm5.append(i+1)
#        if iaad == 235: tm6.append(i)
#        if iaad == 263: tm6.append(i+1)
#        if iaad == 268: tm7.append(i)
#        if iaad == 293: tm7.append(i+1)
#
#    for tm in [tm1, tm2, tm3, tm4, ecl2, tm5, tm6, tm7]:
#        name_list = []
#        tm_list = []
#        n_term = []
#        for elem in new_list:
#            name = elem[0]
#            seq = elem[1]
#            name_list.append(name)
#            tm_seq = seq[tm[0]:tm[1]]
#            tm_list.append(tm_seq)
#        tm_list = drop_gap(tm_list)
#        count_dict = {n: c.count('-') for n, c in zip(name_list, tm_list)}
#        for k, v in count_dict.items():
#            print(k, ': ', v)
#        print("#"*15)
#    for name, seq in new_list:
#        nterm = seq[:tm1[0]].replace('-', '')
#        print(name, ': ', len(nterm))
#        print('*'*15)
