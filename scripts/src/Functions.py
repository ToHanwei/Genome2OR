#!coding:utf-8

import re
import os
import time
import logging
import platform
from collections import Counter
from collections import defaultdict
from typing import List, Any

from src.config import *
from src.CodeMessages import NtermError
from src.CodeMessages import StrandError
from src.CodeMessages import LengthError
from src.CodeMessages import PlatformError
from src.CodeMessages import FastaFormatError
from src.CodeMessages import VersionWarning


def logfun(func):
    """Function run time modifier"""

    def logtimer(*args, **kwargs):
        """Record the running time of the function"""
        fname = func.__name__
        start = time.time()
        logging.info('function %s() is running' % fname)
        temp = func(*args, **kwargs)
        end = time.time()
        runtime = end - start
        logging.info('{} use {:.3f} seconds'.format(fname, runtime))
        return temp

    return logtimer


def logger(script):
    logfile = "RunRecord.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
        datefmt='%Y-%m-%d %A %H:%M:%S',
        filename=logfile,
        filemode='a'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s '
                                  + '%(filename)s:'
                                  + '%(levelname)s '
                                  + '%(message)s'
                                  )
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)
    logging.info("###" + script + " program starts running###")


@logfun
def run_nhmmer(evalue, cpus, output, profile, genome, verbose):
    """
    Function:
        run_nhmmer
        Run the nhmmer program
    Parameter:
        evalue, report sequences <= this E-value threshold in output
        cpus,  number of parallel CPU workers to use
        output, nhmmer output file
        profile, profile nhmmer need(hmm, [un]alignment file)
        genome, genomic datath
        verbose, show version message and exit [true]
    Return:
        return type -> int
        generate an output file and return exit status
    """

    commandline = "time " \
                  + NHMMER \
                  + " -E " \
                  + str(evalue) \
                  + " --dna --cpu " \
                  + str(cpus) \
                  + " --tblout " \
                  + output \
                  + " " + profile \
                  + " " + genome
    if not verbose:
        commandline += " > /dev/null"
    code = os.system(commandline)
    return code


@logfun
def platform_info(verbose):
    """
    Function:
        platform_info
        Determine the type of operating system, python version.
    Return:
        return type -> String
        system_type, type of operating system
        version, python version
    """
    system_type = platform.system()
    if verbose:
        if system_type == 'Windows':
            print("The system you use is Windows.")
        elif system_type == 'Linux':
            print("The system you use is Linux.")
        elif system_type == 'Darwin':
            print("The system you use is Darwin.")
        elif system_type == 'Java':
            print("The system you use is Java.")
        else:
            print("Other unknown systems ")

    python_version = platform.python_version()
    version = python_version.split(".")[0]
    if verbose:
        if version == "3":
            print("Python3 is used.")
        elif version == "2":
            print("Python2 is used.")
    # check platform
    if system_type != "Linux":
        logging.error("Platform error, use Linux only")
        raise PlatformError(system_type)
    if version != '3':
        logging.warning("This program is best run on python3")
        print(VersionWarning(version))


def reverse_complement(string):
    """
    Function:
        reverse_complement
        get reverse complement of a DNA string
    Parameter:
        string, a DNA sequence
    Return:
        return type -> String
        Reversed DNA sequence
    """

    str_reverse = string.strip()[::-1]
    try:
        str_comp = ''.join([CompBase[base] for base in str_reverse])
    except KeyError:
        str_replace = ''
        for s in string:
            if s in CompBase.keys():
                str_replace += s
            else:
                str_replace += 'N'
        str_comp = reverse_complement(str_replace)
        logging.error('Illegal letters appear in the genome!')
    return str_comp


def find_all(substring, string):
    """
    Function:
        find_all
        find all indexes of a substring in a string.
        Overlapping is considered.
    Parameter:
        substring,
        string,
    Return:
        return type -> List
        substring starting position
    """

    indexes = [m.start() for m in re.finditer(substring, string)]
    return indexes


def dna_translation(dna_seq):
    """
    Function:
        dna_translation
        Translate DNA sequence to protein sequence
    Parameter:
        dna_seq, DNA sequence without header
    Return:
        return type -> String
        Translated protein sequence
    """

    protein = ''
    i = 0
    len_dna = len(dna_seq)
    if len_dna % 3 != 0:
        logging.error("The length of CDS({0}) is not divisible by 3."
                      .format(len_dna))
        raise LengthError(len_dna)
    plen = len(dna_seq) / 3
    while i < plen:
        n = dna_seq[3 * i: 3 * (i + 1)]
        if 'N' in n:
            r = '*'
        else:
            r = CODON_TABLE[n]
        i += 1
        protein += r
    return protein


def pattern_search(string):
    """
    Function:
        pattern_search
        Search "OR" pattern from sequence
    Parameter:
        string, protein sequence
    Pattern:
    '[YFC].....[GAS]N..[ILMV]',    # TM1
    'L..PMY..[LI]',                # TM2
    '[IL]....C..Q',                # TM3-1
    'M..DR..A',                    # TM3-2
    'A[IV]..PL.Y',                 # TM3-3
    'C.........C.........C',       # ECL2
    'SY..[IVL]',                   # TM5
    '[KR]...T[CL]..H',             # TM6
    'P..NP..[YF]'                  # TM7
    Return:
        return type -> List
        positions, patterns match(or dismatch) position list
        return type -> Integer
        nmatch, number of pattern match
    """

    positions = []
    for i in range(len(PATTERNS)):
        pattern = PATTERNS[i]
        pos = GOLD_POS[i]
        match = re.search(pattern, string)
        if match:
            start = match.start()
            if abs(start - pos) > 24:
                positions.append(-2)
            else:
                positions.append(start)
        else:
            positions.append(-1)
        nmatch = len(PATTERNS) - positions.count(-1) - positions.count(-2)
    return positions, nmatch


def protein_aa_comp(fasta):
    """
    Function:
        protein_aa_comp
        Count the number of residues in the sequence
    Parameter:
        fasta, sequence without header
    Return:
        return type -> Dict
        key   -> residue(base)
        value -> Ratio of residue(base)
    """

    count_dict = Counter(fasta)
    prolen = len(fasta)
    return {k: v * 100.0 / prolen for k, v in count_dict.items()}


def extract_cds_match(scaf, hmmout, dna):
    """
    Function:
        extract_cds_match
        Serves ectract_cds function
    Parameter:
        scaf, scaffold name
        hmmout, funtion 'proc_nhmmer_out' output
        dna, scaffold (DNA) squence
    Return:
        return type -> str
        hit, hit gene name from genome
        return type -> str
        cutseq, cut sequence from scaffold
    """

    outdict = {}
    for hit in hmmout:
        if hmmout[hit][0] == scaf:
            new_fr = hmmout[hit][3]
            new_to = hmmout[hit][4]
            sign = hmmout[hit][5]
            cutseq = dna[new_fr:new_to].upper()
            seq_replace = ''
            for s in cutseq:
                if s in ['A', 'T', 'C', 'G', 'N']:
                    seq_replace += s
                else:
                    seq_replace += 'N'
            if cutseq != seq_replace:
                cutseq = seq_replace
                logging.error('Illegal letters appear in the genome!')
            if sign == '-':
                cutseq = reverse_complement(cutseq)
            outdict[hit] = cutseq
    return outdict


@logfun
def extract_cds(hmmout, gefile):
    """
    Function:
         extract_cds
         Extract cds from genomic file
    Parameter:
         hmmout, funtion 'proc_nhmmer_out' output
         hitnames, hit names
    Return:
        return type -> Dict
        key   -> hit gene name
        value -> CDS sequence
    """
    hmmout_seq = {}
    fin = open(gefile)

    ss, seq_line, flag = '', '', False
    line = fin.readline()
    if line[0] != '>':
        logging.error("Your file start with '{0}'"
                      .format(line[0]))
        raise FastaFormatError(line[0])
    while line:
        if line[0:1] == '>':
            if flag:
                d = extract_cds_match(ss, hmmout, seq_line)
                hmmout_seq.update(d)
            seq_line = ''
            ss = line[1:].split()[0]
        else:
            flag = True
            seq_line += line.strip()
        line = fin.readline()
    fin.close()
    return hmmout_seq


@logfun
def proc_nhmmer_out(file, EvalueLimit):
    """
    Function:
        proc_nhmmer_out
        process nhmmer .out file, get sca, fr, to, sign,  & new_fr, new_to
    Parameter:
        file, nhmmer outfile,header like following form
        "target name", "accession", "query name", "accession",
        "hmmfrom", "hmm to", "alifrom", "ali to", "envfrom",
        "env to", "sq len", "strand", "E-value", "score", "bias",
        "description of target"
        EvalueLimit, Sequence similarity threshold
    Return:
        return type (hmmout) -> Dictionary
        key   -> gene_name, hit gene name
        value -> a out list
        return type (hitname) -> List
        hit name (genome sequence name)
    """

    with open(file) as f:
        linelist = [line for line in f if line[0:1] != '#']

    hmmout = {}
    hitname, hmmfr, hmmto, sqlen, evalue = [], [], [], [], []
    for line in linelist:
        temp = line.split()
        sca = temp[0]
        hmmfr = int(temp[4])
        hmmto = int(temp[5])
        hmmlen = hmmto - hmmfr
        envfr = int(temp[8])
        envto = int(temp[9])
        slen = int(temp[10])
        sign = temp[11]
        evalue = float(temp[12])
        gene_name = sca + '_' + str(envfr) + '_' + str(envto) + '_' + sign
        extend = 1000 - abs(envto - envfr)
        # sequence similarity and length up to grade
        if (evalue < EvalueLimit) and (hmmlen > 599):
            hitname.append(sca)
            if sign == '+':
                new_to = min(envto + extend, slen)
                if envfr < extend:
                    new_fr = 0
                else:
                    new_fr = envfr - extend - 1
                hmmout[gene_name] = [
                    sca, envfr, envto, new_fr,
                    new_to, sign, evalue, slen, hmmlen
                ]
            elif sign == '-':
                new_fr = max(0, envto - extend - 1)
                if (slen - envfr) < extend + 1:
                    new_to = slen
                else:
                    new_to = envfr + extend
                hmmout[gene_name] = [
                    sca, envto, envfr, new_fr,
                    new_to, sign, evalue, slen, hmmlen
                ]
            else:
                # sign(strand) must be '+' or '-'
                logging.error("NHMMER tool output 'strand' "
                              + "column only '+' or '-'")
                raise StrandError(sign)

    return hmmout, hitname


def find_stop_codons(seq, SeqLengthLimit):
    """
    Function:
        find_stop_codons
        Find stop codons in seq, no matter what reading frame is
    Parameter:
        seq, DNA sequence without header
    Return:
        return type -> List
        All 'stop codon' of seq(input DNA sequence)
    """

    stop_tag = find_all('(?=(TAG))', seq)
    stop_tga = find_all('(?=(TGA))', seq)
    stop_taa = find_all('(?=(TAA))', seq)
    stops = stop_tag + stop_tga + stop_taa
    # Filter stop long enough
    stops = sorted(i for i in stops if i > SeqLengthLimit)
    return stops


@logfun
def find_cds(hmmout, hmmout_seq, SeqLengthLimit):
    """
    Function:
        find_cds
        try find ATG and STOP codons for each seq,
        put good cds info into fun = {},
        put others into outliers = {} for later processing
    Parameter:
        hmmout, funtion 'proc_nhmmer_out' output
        hmmout_seq, function 'extract_cds' outpuy
        SeqLengthLimit, artificially set OR's sequence length threshold
    Return:
    return type(fun) -> Dictory
    key   -> hit gene name
    value -> function OR CDS
    return type(outliers) -> Dictory
    key   -> hit gene name
    value -> outlier OR CDS
    """

    hit = 1
    fun = {}
    outliers = {}
    for gname in hmmout_seq:
        processed = 0
        raw_cds = hmmout_seq[gname]
        len_cds = len(raw_cds)
        # The CDS length is too short
        if len_cds < SeqLengthLimit:
            outliers[gname] = raw_cds
            continue
        [sca, ori_fr, ori_to, fr, to, strand, evalue, slen, hmmlen] = hmmout[gname]
        atgs = find_all('(?=(ATG))', raw_cds)
        starts = sorted(i for i in atgs if i < len_cds - SeqLengthLimit)
        stops = find_stop_codons(raw_cds, SeqLengthLimit)
        iso = 1
        for iatg in starts:
            for istop in stops:
                klen = istop - iatg
                if (klen > SeqLengthLimit) and (klen % 3 == 0):
                    fine_cds = raw_cds[iatg:istop]
                    a_list = re.findall('...', fine_cds)
                    # interrupting stop codon
                    if ('TAG' in a_list) \
                            or ('TGA' in a_list) \
                            or ('TAA' in a_list):
                        continue
                    processed = 1
                    if strand == '+':
                        real_fr = int(fr) + iatg
                        real_to = int(fr) + istop + 3
                    elif strand == '-':
                        real_fr = int(to) - istop - 3
                        real_to = int(to) - iatg
                    if 'N' in raw_cds[iatg:istop]:
                        outliers[gname] = raw_cds
                    else:
                        hits = "hit" + str(hit)
                        new_name = sca + '_' \
                                   + str(real_fr) + '_' \
                                   + str(real_to) + '_' \
                                   + strand + '_iso' \
                                   + str(iso)
                        final_cds = raw_cds[iatg:(istop + 3)]
                        fun[new_name] = [
                            hits, gname, iso, final_cds,
                            real_fr, real_to, strand
                        ]
                        iso += 1
                        break
        if processed == 0:
            outliers[gname] = raw_cds
        hit += 1
    return fun, outliers


def sequence_align(indir, seqf, alignf):
    """
    Function:
        sequence_align
        call MAFFT(linsi) multiple sequence alignment
    Parameter:
        seqf, sequence file name (FASTA format)
        alignf, Multiple sequence alignment results
    Return: None
    """

    seqf = os.path.join(indir, seqf)
    command = MAFFT \
              + " --localpair " \
              + "--maxiterate 1000 " \
              + "--quiet " \
              + seqf \
              + " > " \
              + alignf
    os.system(command)


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


@logfun
def refact_hitfile(hitfile):
    """
    Function:
        refact_hitfile
        Refactoring hit file
    Parameter:
        hitfile, OR sequence annotated from genome
    Return:
        return type -> Dict
        hit_dict, output refactoring dict
    """
    with open(ORfile) as ORf:
        template = ORf.readlines()

    hit_dict = defaultdict(list)
    seqlist = ReadSampleFasta(hitfile)
    for header, seq in seqlist:
        key = header.split('_')[0]
        hit_dict[key].append((header, seq))
    return template, hit_dict


@logfun
def refact_list(template, hit_dict):
    """
    Function:
        refact_list
        Reorganization of the sequence list.
        Note: Return a generator.
    Parameter:
        template, OR template sequence list
        hit_dict, hit OR sequence list
    Return:
        return type -> list
        seq_list, sequence alignment list
    """

    for hit, seqs in hit_dict.items():
        templist = []
        for head, seq in seqs:
            head = ">" + head
            templist.extend([head, seq])
        aligns = template + templist
        Writer(TEMP, TEMP_OR, aligns)
        sequence_align(TEMP, TEMP_OR, TEMP_ALIGN)
        seq_list = ReadSampleFasta(TEMP_ALIGN)
        # drop some template OR sequence, retain OR5AN1 only
        seq_list = seq_list[0:1] + seq_list[4:]
        yield seq_list


def tm_cut(seq_list):
    """
    Function:
        tm_cut
        cut sequence to TMs(and ICL, ECL, N-term. H8)
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
        if iaad == TM_boundary[0]: tm1.append(i)
        if iaad == TM_boundary[1]: tm1.append(i+1)
        if iaad == TM_boundary[2]: tm2.append(i)
        if iaad == TM_boundary[3]: tm2.append(i+1)
        if iaad == TM_boundary[4]: tm3.append(i)
        if iaad == TM_boundary[5]: tm3.append(i+1)
        if iaad == TM_boundary[6]: tm4.append(i)
        if iaad == TM_boundary[7]: tm4.append(i+1)
        if iaad == TM_boundary[8]: tm5.append(i)
        if iaad == TM_boundary[9]: tm5.append(i+1)
        if iaad == TM_boundary[10]: tm6.append(i)
        if iaad == TM_boundary[11]: tm6.append(i+1)
        if iaad == TM_boundary[12]: tm7.append(i)
        if iaad == TM_boundary[13]: tm7.append(i+1)
    icl1, icl2, icl3 = [], [], []
    ecl1, ecl2, ecl3 = [], [], []
    name_list, cut_list, nterm_list = [], [], []
    h8_list = []
    hit_list = seq_list[1:]
    for name, seq in hit_list:
        # append n-term sequence to list
        nterm_list.append(seq[:tm1[0]].replace('-', ''))
        h8_list.append(seq[tm7[1]:].rstrip())
        tm_list = []
        for tm in [tm1, tm2, tm3, tm4, tm5, tm6, tm7]:
            tm_seq = seq[tm[0]:tm[1]]
            tm_list.append(tm_seq)
        cut_list.append(tm_list)
        name_list.append(name)
        icl1.append(seq[tm1[1]:tm2[0]])
        icl2.append(seq[tm3[1]:tm4[0]])
        ecl2.append(seq[tm4[1]:tm5[0]])
    nongap_list = [drop_gap(tms) for tms in zip(*cut_list)]
    nongap_cuts = list(zip(*nongap_list))
    cut_dict = dict(zip(name_list, nongap_cuts))
    nterm_dict = dict(zip(name_list, nterm_list))
    icl1 = dict(zip(name_list, icl1))
    icl2 = dict(zip(name_list, icl2))
    ecl2 = dict(zip(name_list, ecl2))
    h8 = dict(zip(name_list, h8_list))
    return nterm_dict, cut_dict, icl1, icl2, ecl2, h8


def tm_pattern(seq_list):
    """
    Function:
        tm_pattern
        match pattern from TMs count
    Parameter:
        seq_list, sequence alignment list
    Return:
        return type -> list
        pseu, pseudogene OR gene names list
        nterms, N-term list [(seq_name, N-term), ...]
        tm_lists, TM list [(seq_name, (tm1, tm2, ..., tm7), ...]
    """

    pseu = []
    nterms = []
    tm_lists = []
    nterm, tm_dict, icl1, icl2, ecl2, h8 = tm_cut(seq_list)
    for k, tms in tm_dict.items():
        seqs = [tms[0], tms[1], tms[1],
                tms[2], tms[2], tms[4],
                tms[5], tms[6], icl1[k],
                icl2[k], ecl2[k],
                nterm[k], h8[k],
                ]
        seqs = [seq.replace('-', '') for seq in seqs]
        n = sum(1 for pat, seq in zip(patterns, seqs) if re.search(pat, seq))
        if n < PATTERN_THRESHOLD:
            pseu.append(k)
        else:
            nterms.append((k, nterm[k]))
            tm_lists.append((k, tms))
    return pseu, nterms, tm_lists


def Nterm_length(nterms):
    """
    Function:
        Nterm_length
        choose function OR gene through N-term length
    Parameter:
        nterm_list, Every element is a (name, N-term)
    Return:
        return type -> list
        funcs, function OR gene names list
        pseu, pseudogene OR gene names list
    """

    funcs, pseus = [], []
    other, area_a = [], []
    area_b1, area_b2 = [], []
    area_c1, area_c2 = [], []
    for name, nterm in nterms:
        length = len(nterm)
        if 13 <= length <= 16:
            area_c1.append((name, length))
        elif 17 <= length <= 19:
            area_b1.append((name, length))
        elif 20 <= length <=25:
            area_a.append((name, length))
        elif 26 <= length <= 35:
            area_b2.append((name, length))
        elif 36 <= length <= 47:
            area_c2.append((name, length))
        else:
            other.append(name)

    # Processing area A
    if len(area_a) != 0:
        a_shift = [(n, abs(l-22)) for n, l in area_a]
        a_name = sorted(a_shift, key=lambda x: x[1])[0][0]
        funcs.append(a_name)
    # Processing area B
    elif len(area_b1) or len(area_b2):
        if (len(area_b1) > 0) and (len(area_b2) > 0):
            b2 = sorted(area_b2, key=lambda x: x[1])[0]
            b1 = sorted(area_b1, key=lambda x: x[1], reverse=True)[0]
            if abs(b2[1] - 25) <= abs(b1[1] - 20):
                funcs.append(b2[0])
            else:
                funcs.append(b1[0])
        elif (len(area_b1) > 0) and (len(area_b2) == 0):
            b1 = sorted(area_b1, key=lambda x: x[1], reverse=True)[0]
            funcs.append(b1[0])
        else:
            b2 = sorted(area_b2, key=lambda x: x[1])[0]
            funcs.append(b2[0])
    # Processing area C
    elif len(area_c1) or len(area_c2):
        if (len(area_c1) > 0) and (len(area_c2) > 0):
            c2 = sorted(area_c2, key=lambda x: x[1])[0]
            c1 = sorted(area_c1, key=lambda x: x[1], reverse=True)[0]
            if abs(c2[1] - 35) <= abs(c1[1] - 17):
                funcs.append(c2[0])
            else:
                funcs.append(c1[0])
        elif (len(area_c1) > 0) and (len(area_c2) == 0):
            c1 = sorted(area_c1, key=lambda x: x[1], reverse=True)[0]
            funcs.append(c1[0])
        else:
            c2 = sorted(area_c2, key=lambda x: x[1])[0]
            funcs.append(c2[0])
    elif len(other):
        other = str(other)
        logging.info("{0}'s N-term out of range!".format(other))
    else:
        logging.info("process an empty list!")

    if len(funcs) == 0:
        pseus = [name for name, _ in nterms]
    return funcs, pseus


def tm_gaps_filter(seq_list):
    """
    Function:
        tm_gaps_filter
        The sequence is filtered by the number of gaps in region TM
    Parameter:
        seq_list, tms list (seq_name, tms), 'tms' is a tuple save s tms
    Return:
        return type -> list
        funcs, function OR gene names list
        pseu, pseudogene OR gene names list
    """

    func, pseu = [], []
    for name, tms in seq_list:
        gaps = [tm.count('-') for tm in tms if tm.count('-') > 0]
        gap_tms = len(gaps)
        gap_total = sum(gaps)
        # no more than 2 tm contain gao, and total gap less than 5
        if (gap_tms <= 2) and (gap_total <= 5):
            func.append(name)
    if len(func) == 0:
        pseu = [name for name, tms in seq_list]
    return func, pseu


def calculate_offset(sum_d, nmatch):
    """
    Function:
       calculate_offset
       calculate average pattern offset
    Parameter:
        sum_d, pattern total offset
        nmatch, number of pattern match
    Return:
        return type -> Float
        avg_d, average pattern offset
    """

    npattern = len(PATTERNS)
    try:
        # it's mean nmatch=0, no pattern match
        avg = sum_d / nmatch
    except ZeroDivisionError:
        avg = float('inf')
    # mismatch pattern, fill with average
    avg_d = (sum_d + avg * (npattern - nmatch)) / npattern
    avg_d = round(avg_d, 2)
    return avg_d


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


def writer2file(outdir, prefix, funcs, outliers, hmmout_seq, hmmout):
    """
    Function:
        Writer
        write output data to file.
    Parameter:
        outdir, the directory where the output file is saved
        prefix, output file prefix.
        funcs, find_cds() return, function OR dict
        outliers, find_cds() return
        hmmout_seq, extract_cds() return
        hmmout, proc_nhmmer_out() return
    Return:
        None
    """

    count = 0
    dna_list, pro_list, summary_list = [], [], []
    s_header = 'count\tsca\tfr\tto\treal_fr\treal_to\tstrand\tevalue\thmmlen\tlen_cds\n'
    summary_list.append(s_header)
    for index in funcs:
        [hit, old_i, iso, cds, real_fr, real_to, _a] = funcs[index]
        [sca, fr, to, _b, _c, sign, evalue, _d, hmmlen] = hmmout[old_i]
        header = [hit, sca, str(real_fr), str(real_to), sign, 'iso' + str(iso)]
        header = "_".join(header)
        seq = dna_translation(cds)
        dna_list.append('>' + header + '\n' + cds + '\n')
        pro_list.append('>' + header + '\n' + seq + '\n')
        count += 1
        summary_line = [
            count, sca,
            fr, to,
            real_fr, real_to,
            sign, evalue,
            hmmlen, len(cds)
        ]
        summary_line = "\t".join([str(line) for line in summary_line])
        summary_list.append(summary_line + "\n")
    outlier_list = ['>' + name + '\n' + line + '\n' for name, line in outliers.items()]

    Writer(outdir, prefix + "_ORs_DNA.fa", dna_list)
    Writer(outdir, prefix + "_ORs_pro.fa", pro_list)
    Writer(outdir, prefix + "_outliers_dna.fa", outlier_list)
    Writer(outdir, prefix + "_summary_cds.txt", summary_list)

    # result write to log file
    num_cds = len(hmmout_seq)
    num_out = len(outliers)
    num_fun_cds = len(hmmout_seq) - len(outliers)
    num_fun_OR = len(funcs)
    num_iso = num_fun_OR - num_fun_cds
    cds = "{} OR CDS found by nhmmer. Among them." \
        .format(num_cds)
    outlie = "{} outliers cds were found." \
        .format(num_out)
    funcds = "{} functional cds were found;" \
        .format(num_fun_cds)
    funor = "generating {} functional ORs;" \
        .format(num_fun_OR)
    isoform = "including {} isoforms".format(num_iso)
    logging.info("###The result as follows###")
    logging.info(prefix + " processing completed")
    logging.info(cds)
    logging.info(outlie)
    logging.info(funcds)
    logging.info(funor)
    logging.info(isoform)
    logging.info("###The program finish###")


if __name__ == "__main__":
    platform_info(True)
