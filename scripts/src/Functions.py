#!coding:utf-8

import re
import os
import stat
import sys
import time
import pickle
import logging
import tempfile
import platform
from multiprocessing import Process
from multiprocessing import Queue
from multiprocessing import Lock
from collections import Counter
from collections import defaultdict

from src.config import *
from src.CodeMessages import NtermError
from src.CodeMessages import FileNotExists
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

    commandline = ("time "
                   + NHMMER
                   + " -E "
                   + str(evalue)
                   + " --dna --cpu "
                   + str(cpus)
                   + " --noali"
                   + " --popen 0.02"  # default 0.02
                   + " --pextend 0.4"  # default 0.4
                   + " --tblout "
                   + output
                   + " " + profile
                   + " " + genome
                   )
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
    sourdict = {}
    hmmout = {k: v for k, v in hmmout.items() if v[0] == scaf}
    for hit in hmmout:
        fr, to = sorted(hmmout[hit][1:3])
        new_fr = hmmout[hit][3]
        new_to = hmmout[hit][4]
        sign = hmmout[hit][5]
        cutseq = dna[new_fr:new_to].upper()
        sourseq = dna[fr:to].upper()
        seq_replace = ''
        sour_replace = ''
        for s in cutseq:
            # Other non-standard bases are converted to N
            if s in ['A', 'T', 'C', 'G', 'N']:
                seq_replace += s
            else:
                seq_replace += 'N'
        for s in sourseq:
            # Other non-standard bases are converted to N
            if s in ['A', 'T', 'C', 'G', 'N']:
                sour_replace += s
            else:
                sour_replace += 'N'
        cutseq = seq_replace
        sourseq = sour_replace
        if sign == '-':
            cutseq = reverse_complement(cutseq)
            sourseq = reverse_complement(sourseq)
        outdict[hit] = cutseq
        sourdict[hit] = sourseq
    return outdict, sourdict


@logfun
def extract_cds(hmmout, gefile):
    """
    Function:
         extract_cds
         Extract cds from genomic file
    Parameter:
         hmmout, function 'proc_nhmmer_out' output
         gefile, genome sequence file
    Return:
        return type -> Dict
        key   -> hit gene name
        value -> CDS sequence
    """
    hmmout_seq = {}
    sour_seq = {}
    genomef = open(gefile)
    header, seq_line, flag = '', '', False
    line = genomef.readline()
    if line[0] != '>':
        logging.error("Your genome start with '{0}', FASTA format?"
                      .format(line[0]))
        raise FastaFormatError(line[0])
    while line:
        if line[0] == '>':
            if flag:
                cds, sourcds = extract_cds_match(header, hmmout, seq_line)
                hmmout_seq.update(cds)
                sour_seq.update(sourcds)
                flag = False
            seq_line = ''
            header = line[1:].split()[0]
        else:
            flag = True
            seq_line += line.strip()
        line = genomef.readline()
    # Process last sequence in genome
    cds, sourcds = extract_cds_match(header, hmmout, seq_line)
    hmmout_seq.update(cds)
    sour_seq.update(sourcds)

    genomef.close()
    return hmmout_seq, sour_seq


@logfun
def proc_nhmmer_out(file, EvalueLimit, SeqLengthLimit):
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
        SeqLengthLimit, Sequence length threshold
    Return:
        return type (hmmout) -> Dictionary
        key   -> gene_name, hit gene name
        value -> a out list
    """

    with open(file) as tbloutf:
        # Filter the header and tail
        linelist = [line for line in tbloutf if line[0:1] != '#']

    hmmout, trunc = {}, {}
    for line in linelist:
        temp = line.strip().split()
        # Extract information from NHMMER outfile
        sca = temp[0]
        hmmfr = int(temp[4])
        hmmto = int(temp[5])
        envfr = int(temp[8])
        envto = int(temp[9])
        slen = int(temp[10])
        sign = temp[11]
        evalue = float(temp[12])
        hmmlen = abs(hmmto - hmmfr) + 1
        envlen = abs(envto - envfr) + 1
        gene_name = sca + '_' + str(envfr) + '_' + str(envto) + '_' + sign
        # if CDS length less than EXTEND_LENGTH, extend it
        extend = EXTEND_LENGTH - envlen
        if extend > 0:
            stop_extend = int(0.5 * extend)
            start_extend = extend - stop_extend
        else:
            stop_extend = 0
            start_extend = 0
        # sequence similarity filter
        if evalue > EvalueLimit: continue
        if sign == '+':
            new_fr = max(envfr - start_extend, 0)
            new_to = min(envto + stop_extend, slen)
        elif sign == '-':
            # The chain of antisense needs to be reversed
            new_fr = max(envto - start_extend, 0)
            new_to = min(envfr + stop_extend, slen)
        else:
            # sign(strand) must be '+' or '-'
            logging.error("NHMMER tool output 'strand' "
                          + "column only '+' or '-'")
            raise StrandError(sign)
        fraglen = abs(new_to - new_fr) 
        # fraglen equal to EXTEND_LENGTH, otherwise truncat
        if fraglen < SeqLengthLimit:
            trunc[gene_name] = [ 
                sca, envfr, envto, new_fr,
                new_to, sign, evalue, slen, hmmlen
            ]
        else:
            hmmout[gene_name] = [
                sca, envfr, envto, new_fr,
                new_to, sign, evalue, slen, hmmlen
            ]
    truncdoc = ("\033[0;32m{}\033[0m truncated gene(s) was discovered"
                .format(len(trunc)))
    logging.info(truncdoc)
    return hmmout, trunc


def find_stop_codons(seq, SeqLengthLimit):
    """
    Function:
        find_stop_codons
        Find stop codons in seq, no matter what reading frame is
    Parameter:
        seq, DNA sequence without header
        SeqLengthLimit, artificially set OR's sequence length threshold
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


def cds_length_filter(cdslist, SeqLengthLimit):
    """
    Function:
        cds_length_filter
        Determine if all CDS's meet the length limit
    Parameter:
        cdslist, cds list from a DNA fragment
    Return:
        return type -> list
        Filtered cdslist
    """
    filterlist = [cds for cds in cdslist if len(cds[2]) >= SeqLengthLimit]
    return filterlist


def insert_filter(cdslist):
    """
    Function:
        insert_filter
        Determine if all CDS's ware inserted or deleted
    Parameter:
        cdslist, cds list from a DNA fragment
    Return:
        return type -> list
        Filtered cdslist
    """
    filterlist = [cds for cds in cdslist if len(cds[2]) % 3 == 0]
    return filterlist


def interrupt_stop_codon(cdslist):
    """
    Function:
        interrupt_stop_codon
        Determine if all CDS ware interrupting stop codon
    Parameter:
        cdslist, cds list from a DNA fragment
    Return:
        return type -> list
        Filtered cdslist
    """
    filterlist = []
    for cds in cdslist:
        codons = re.findall('...', cds[2])
        if 'TAG' in codons: continue
        if 'TGA' in codons: continue
        if 'TAA' in codons: continue
        filterlist.append(cds)
    return filterlist


@logfun
def find_cds(hmmout, hmmout_seq, sour_seq, SeqLengthLimit):
    """
    Function:
        find_cds
        try find ATG and STOP codons for each seq,
        put good cds info into fun = {},
    Parameter:
        hmmout, function 'proc_nhmmer_out' output
        hmmout_seq, function 'extract_cds' output
        SeqLengthLimit, artificially set OR's sequence length threshold
    Return:
    return type(funcdict) -> Dictory
    key   -> hit gene name
    value -> assume that OR CDS
    return type(pseudos) -> Dictory
    key   -> hit gene name
    value -> assume that OR pseudogenes
    """

    hit = 1
    pseutype = ''
    funcdict = {}
    pseudos = {}
    lengfilter = 0
    insertfilter = 0
    interrfilter = 0
    func_gnames = []
    pseunum = []
    for gname in hmmout_seq:
        iso = 1
        temp_pseu = []
        funcgene = False
        real_fr, real_to = '', ''
        raw_cds = hmmout_seq[gname]
        sour_cds = sour_seq[gname]
        len_cds = len(raw_cds)
        [sca, _a, _b, fr, to, strand, _c, _d, _e] = hmmout[gname]
        starts = find_all('(?=(ATG))', raw_cds)
        starts = sorted(i for i in starts if i < len_cds - SeqLengthLimit)
        stops = find_stop_codons(raw_cds, SeqLengthLimit)
        cdslist = [(i, j, raw_cds[i:j]) for i in starts for j in stops]
        # Determine if all CDS's meet the length limit
        cdslist = cds_length_filter(cdslist, SeqLengthLimit)
        if cdslist:
            # Insert or delete codon (pseudogene)
            cdslist = insert_filter(cdslist)
            if cdslist:
                # Interrupting stop codon (pseudogene)
                cdslist = interrupt_stop_codon(cdslist)
                if cdslist:
                    for cds in cdslist:
                        iatg, istop, cds_seq = cds
                        if strand == '+':
                            real_fr = int(fr) + iatg
                            real_to = int(fr) + istop + 3
                        elif strand == '-':
                            real_fr = int(to) - iatg
                            real_to = int(to) - istop - 3
                        if 'N' not in cds_seq:
                            hits = "hit" + str(hit)
                            new_name = (sca + '_'
                                        + str(real_fr) + '_'
                                        + str(real_to) + '_'
                                        + strand + '_iso'
                                        + str(iso)
                                        )
                            final_cds = raw_cds[iatg:(istop + 3)]
                            funcdict[new_name] = [
                                hits, gname, iso, final_cds,
                                real_fr, real_to, strand
                            ]
                            func_gnames.append(gname)
                            funcgene = True
                            iso += 1
                else:
                    # all cds fragment has interrupt stop codon
                    pseutype = "INTER"
                    interrfilter += 1
                    
            else:
                # all cds fragment can not be multiple by 3
                pseutype = "INDEL"
                insertfilter += 1
        else:
            # all cds fragment length too short
            pseutype = "SHORT"
            lengfilter += 1
        if not funcgene:
            new_gname = gname + '_' + pseutype
            pseudos[new_gname] = sour_cds
        else:
            hit += 1
    pseunum = (lengfilter, insertfilter, interrfilter)
    return funcdict, pseudos, pseunum


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

    command = (MAFFT
               + " --localpair "
               + "--maxiterate 1000 "
               + " --thread -1 "
               + "--quiet "
               + seqf
               + " > "
               + alignf
               )
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

    with open(seqfile) as seqf:
        seq_list = []
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
    with open(OR_TEMPLATE) as ORf:
        template = ORf.readlines()

    hit_dict = defaultdict(list)
    seqlist = ReadSampleFasta(hitfile)
    for header, seq in seqlist:
        key = header.split('_')[0]
        hit_dict[key].append((header, seq))
    return template, hit_dict


def tempfile_align(seqs, template, queue, lock):
    """
    Function:
        tempfile_align
        Alignment sequence with tempfile
    Parameter:
        seqs: sequences list, [(name, sequence), ...]
        template: templact sequences list
        queue: multiprocessing.Queue, save result
    """
    templist = []
    for head, seq in seqs:
        head = ">" + head
        templist.extend([head, seq])
    aligns = template + templist
    # temporary OR sequence file name
    TempOR = tempfile.NamedTemporaryFile('w+t')
    TempName = TempOR.name
    TempOR.writelines(aligns)
    TempOR.seek(0)
    # temporary OR alignment file name
    TempAlign = TempName + '.fas'
    sequence_align(TempName, TempAlign)
    seq_list = ReadSampleFasta(TempAlign)
    # drop some template OR sequence, retain OR5AN1 only
    index = int(len(template) / 2)
    seq_list = seq_list[0:1] + seq_list[index:]
    os.remove(TempAlign)
    TempOR.close()
    lock.acquire()
    queue.put(seq_list)
    lock.release()


@logfun
def refact_list(template, hit_dict, cpus):
    """
    Function:
        refact_list
        Reorganization of the sequence list.
        Note: Return a generator.
    Parameter:
        template, OR template sequence list
        hit_dict, hit OR sequence list
        cpus, number of parallel
    Return:
        return type -> list
        seq_list, sequence alignment list
    """

    seq_lists = []
    keys = list(hit_dict.keys())
    num_of_jobs = len(keys)
    for i in range(0, num_of_jobs, cpus):
        batch_jobs = keys[i:i+cpus]
        queue = Queue()
        lock = Lock()
        jobs = []
        for hit in batch_jobs:
            seqs = hit_dict[hit]
            p = Process(
                target=tempfile_align,
                args=(seqs, template, queue, lock)
            )
            p.start()
            jobs.append(p)
        for job in jobs:
            job.join()
        while not queue.empty():
            seq_lists.append(queue.get())
    return seq_lists

#@logfun
#def refact_list(template, hit_dict):
#    """
#    Function:
#        refact_list
#        Reorganization of the sequence list.
#        Note: Return a generator.
#    Parameter:
#        template, OR template sequence list
#        hit_dict, hit OR sequence list
#    Return:
#        return type -> list
#        seq_list, sequence alignment list
#    """
#
#    seq_lists = []
#    for hit, seqs in hit_dict.items():
#        templist = []
#        for head, seq in seqs:
#            head = ">" + head
#            templist.extend([head, seq])
#        aligns = template + templist
#        # temporary OR sequence file name
#        TempOR = tempfile.NamedTemporaryFile('w+t')
#        TempName = TempOR.name
#        TempOR.writelines(aligns)
#        TempOR.seek(0)
#        # temporary OR alignment file name
#        TempAlign = TempName + '.fas'
#        sequence_align(TempName, TempAlign)
#        seq_list = ReadSampleFasta(TempAlign)
#        # drop some template OR sequence, retain OR5AN1 only
#        index = int(len(template) / 2)
#        seq_list = seq_list[0:1] + seq_list[index:]
#        os.remove(TempAlign)
#        TempOR.close()
#        seq_lists.append(seq_list)
#    return seq_lists


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
    tm1, tm2, tm3, tm4, ecl2, tm5, tm6, tm7 = [], [], [], [], [], [], [], []
    for i in range(len(seq_list[0][1])):
        aad = seq_list[0][1][i]
        if aad != "-":
            iaad += 1
        else:
            continue
        if iaad == TM_boundary[0]: tm1.append(i)
        if iaad == TM_boundary[1]: tm1.append(i + 1)
        if iaad == TM_boundary[2]: tm2.append(i)
        if iaad == TM_boundary[3]: tm2.append(i + 1)
        if iaad == TM_boundary[4]: tm3.append(i)
        if iaad == TM_boundary[5]: tm3.append(i + 1)
        if iaad == TM_boundary[6]: tm4.append(i)
        if iaad == TM_boundary[7]: tm4.append(i + 1)
        if iaad == TM_boundary[8]: tm5.append(i)
        if iaad == TM_boundary[9]: tm5.append(i + 1)
        if iaad == TM_boundary[10]: tm6.append(i)
        if iaad == TM_boundary[11]: tm6.append(i + 1)
        if iaad == TM_boundary[12]: tm7.append(i)
        if iaad == TM_boundary[13]: tm7.append(i + 1)
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
    # if tm_dict none, do nothing
    if len(tm_dict) == 0:
        return pseu, nterms, tm_lists
    for k, tms in tm_dict.items():
        seqs = [tms[0], tms[1], tms[1],
                tms[2], tms[2], tms[4],
                tms[5], tms[6], icl1[k],
                icl2[k], ecl2[k],
                nterm[k], h8[k],
                ]
        seqs = [seq.replace('-', '') for seq in seqs]
        n = sum(1 for pat, seq in zip(PATTERNS, seqs) if re.search(pat, seq))
        if n < PATTERN_THRESHOLD:
            # DNA similarity but pattern not match, means missense mutation.
            pseu.append((n, k))
        else:
            nterms.append((k, nterm[k]))
            tm_lists.append((k, tms, n))
    if len(tm_lists) == 0:
        pseu = sorted(pseu, key=lambda x:x[0], reverse=True)
        pseu = pseu[0][1]
    else:
        pseu = None
    return pseu, nterms, tm_lists


def tm_prepare(seq_list):
    """
    Function:
        tm_prepare
        prepare tm list from seq_list
    Parameter:
        seq_list, sequence alignment list
    Return:
        return type -> list
        nterms, N-term list [(seq_name, N-term), ...]
        tm_lists, TM list [(seq_name, (tm1, tm2, ..., tm7), ...]
    """

    nterms = []
    tm_lists = []
    nterm, tm_dict, icl1, icl2, ecl2, h8 = tm_cut(seq_list)
    # if tm_dict none, do nothing
    if len(tm_dict) == 0:
        return nterms, tm_lists
    for k, tms in tm_dict.items():
        seqs = [tms[0], tms[1], tms[1],
                tms[2], tms[2], tms[4],
                tms[5], tms[6], icl1[k],
                icl2[k], ecl2[k],
                nterm[k], h8[k],
                ]
        seqs = [seq.replace('-', '') for seq in seqs]
        n = sum(1 for pat, seq in zip(PATTERNS, seqs) if re.search(pat, seq))
        nterms.append((k, nterm[k]))
        tm_lists.append((k, tms, n))
    return nterms, tm_lists


def Nterm_length(nterms, tm_list):
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

    func, pseu = None, None
    nterm_len = []
    stand_nterm = False
    for name, nterm in nterms:
        leng = len(nterm)
        nterm_len.append((name, leng))
        if leng > 0: stand_nterm = True
    if stand_nterm:
        nterm_len = [(name, abs(leng-23)) for name, leng in nterm_len]
        func = sorted(nterm_len, key=lambda x:x[1])[0][0]
    else:
        tm_list_sort = sorted(tm_list, key=lambda x:x[2], reverse=True)
        pseu = tm_list_sort[0][0]
    #other, area_a = [], []
    #area_b1, area_b2 = [], []
    #area_c1, area_c2 = [], []
    #for name, nterm in nterms:
    #    length = len(nterm)
    #    if 13 <= length <= 16:
    #        area_c1.append((name, length))
    #    elif 17 <= length <= 19:
    #        area_b1.append((name, length))
    #    elif 20 <= length <= 25:
    #        area_a.append((name, length))
    #    elif 26 <= length <= 35:
    #        area_b2.append((name, length))
    #    elif 36 <= length <= 47:
    #        area_c2.append((name, length))
    #    else:
    #        other.append(name)

    ## Processing area A
    #if len(area_a) != 0:
    #    a_shift = [(n, abs(l - 22)) for n, l in area_a]
    #    a_name = sorted(a_shift, key=lambda x: x[1])[0][0]
    #    funcs.append(a_name)
    ## Processing area B
    #elif len(area_b1) or len(area_b2):
    #    if (len(area_b1) > 0) and (len(area_b2) > 0):
    #        b2 = sorted(area_b2, key=lambda x: x[1])[0]
    #        b1 = sorted(area_b1, key=lambda x: x[1], reverse=True)[0]
    #        if abs(b2[1] - 25) <= abs(b1[1] - 20):
    #            funcs.append(b2[0])
    #        else:
    #            funcs.append(b1[0])
    #    elif (len(area_b1) > 0) and (len(area_b2) == 0):
    #        b1 = sorted(area_b1, key=lambda x: x[1], reverse=True)[0]
    #        funcs.append(b1[0])
    #    else:
    #        b2 = sorted(area_b2, key=lambda x: x[1])[0]
    #        funcs.append(b2[0])
    ## Processing area C
    #elif len(area_c1) or len(area_c2):
    #    if (len(area_c1) > 0) and (len(area_c2) > 0):
    #        c2 = sorted(area_c2, key=lambda x: x[1])[0]
    #        c1 = sorted(area_c1, key=lambda x: x[1], reverse=True)[0]
    #        if abs(c2[1] - 35) <= abs(c1[1] - 17):
    #            funcs.append(c2[0])
    #        else:
    #            funcs.append(c1[0])
    #    elif (len(area_c1) > 0) and (len(area_c2) == 0):
    #        c1 = sorted(area_c1, key=lambda x: x[1], reverse=True)[0]
    #        funcs.append(c1[0])
    #    else:
    #        c2 = sorted(area_c2, key=lambda x: x[1])[0]
    #        funcs.append(c2[0])
    #elif len(other):
    #    other = str(other)
    #    logging.info("{0}'s N-term out of range!".format(other))
    #else:
    #    logging.info("process an empty list!")

    #if len(funcs) == 0:
    #    pseus = [name for name, _ in nterms]
    return pseu, func


def tm_gaps_filter(seq_list, nterms):
    """
    Function:
        tm_gaps_filter
        The sequence is filtered by the number of gaps in region TM
    Parameter:
        nterms, N-term list (seq_name, n-term_sequences)
        seq_list, tms list (seq_name, tms), 'tms' is a tuple save s tms
    Return:
        return type -> list
        funcs, function OR gene names list
        pseu, pseudogene OR gene names list
    """

    funcs, pseus = [], []
    pseu = None
    for name, tms, npattern in seq_list:
        gaps = [tm.count('-') for tm in tms if tm.count('-') > 0]
        gap_tms = len(gaps)
        gap_total = sum(gaps)
        # total gap less than 5
        if gap_total <= TM_GAPS_TOTAL:
            funcs.append(name)
        else:
            pseus.append((name, npattern))
    if len(funcs) == 0:
        # a series of fragment select the most npattern as pseuogene
        pseu = sorted(pseus, key=lambda x:x[1], reverse=True)[0][0]
    nterm_filter = [(name, nterm) for name, nterm in nterms if name in funcs]
    return pseu, nterm_filter


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


def writer2file(outdir, prefix, funcs, pseudos, hmmout_seq, hmmout, trunc, pseunum):
    """
    Function:
        Writer
        write output data to file.
    Parameter:
        outdir, the directory where the output file is saved
        prefix, output file prefix.
        funcs, find_cds() return, assume function OR dict
        pseudos, assume pseudos olfactory receptor genes
        hmmout_seq, extract_cds() return
        hmmout, proc_nhmmer_out() return
        trunc, proc_nhmmer_out() return, truncated gene infomation
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
    
    count = 0
    t_header = 'count\tsca\tfr\tto\treal_fr\treal_to\tstrand\tevalue\thmmlen\n'
    trunc_list = [t_header]
    for gname in trunc:
        [sca, fr, to, newfr, newto, sign, evalue, _d, hmmlen] = trunc[gname]
        trunc_line = [
            count, sca,
            fr, to,
            newfr, newto,
            sign, evalue,
            hmmlen]
        trunc_line = "\t".join([str(line) for line in trunc_line])
        trunc_list.append(trunc_line + "\n")

    # writer sequence to file
    logging.info("Merge pseudogene fragement")
    pseudos_list, pseu_type_num = merge_pseudo(pseudos)
    Writer(outdir, prefix + "_Pre-ORs_dna.fa", dna_list)
    Writer(outdir, prefix + "_Pre-ORs_pro.fa", pro_list)
    Writer(outdir, prefix + "_Pre-pseudos_dna.fa", pseudos_list)
    Writer(outdir, prefix + "_summary_cds.txt", summary_list)
    Writer(outdir, prefix + "_truncated.txt", trunc_list)

    # result write to log file
    num_cds = len(hmmout_seq)
    num_fun_OR = len(funcs)
    num_pse_OR = len(pseudos_list)
    merge_pse_num = sum(pseunum) - num_pse_OR
    length_cause_pseu = pseu_type_num['SHORT']
    insert_cause_pseu = pseu_type_num['INDEL']
    interr_cause_pseu = pseu_type_num['INTER']
    cds = "\033[0;32m{}\033[0m OR fragments found by nhmmer.".format(num_cds)
    funcdoc = "\033[1;32m{}\033[0m functional ORs were discover.".format(num_fun_OR)
    pseudoc = "\033[1;32m{}\033[0m pseudogene ORs were discover.".format(num_pse_OR)
    lengthpseo = ("\033[0;32m{}\033[0m pseudogenes cause by too short sequence length."
                .format(length_cause_pseu))
    insertpseo = ("\033[0;32m{}\033[0m pseudogenes cause by insert or delect base."
                .format(insert_cause_pseu))
    interrpseo = ("\033[0;32m{}\033[0m pseudogenes cause by contains termination codons."
                .format(interr_cause_pseu))
    mergedpseo = ("\033[0;32m{}\033[0m pseudogenes fragment were merged"
                .format(merge_pse_num))
    logging.info("###The result as follows###")
    logging.info(prefix + " processing completed")
    logging.info(cds)
    logging.info(funcdoc)
    logging.info(mergedpseo)
    logging.info(pseudoc)
    logging.info(lengthpseo)
    logging.info(insertpseo)
    logging.info(interrpseo)
    logging.info("###The program finish###")

    # record pseudogene number
    if not os.path.exists("PseudoNumRecode.pkl"):
        fwrite = open("PseudoNumRecode.pkl", 'wb')
        pickle.dump({}, fwrite)
        fwrite.close()
    with open('PseudoNumRecode.pkl', 'rb') as pseudoNumR:
        pseudoNum = pickle.load(pseudoNumR)
    with open('PseudoNumRecode.pkl', 'wb') as pseudoNumW:
        pseudoNum[prefix] = pseu_type_num
        pickle.dump(pseudoNum, pseudoNumW)


def identify_filter(seq_dict, tm_list, nterms):
    func, pseu, ptype = None, None, None
    if tm_list:
        # Some pseudogenoes were filtered out by TM gaps
        gap_pseu, nterm_filter = tm_gaps_filter(tm_list, nterms)
        if nterm_filter:
            # Some pseudogenoes were filtered out by N-term length
            nterm_pseu, function = Nterm_length(nterm_filter, tm_list)
            if function:
                func = function
            else:
                pseu = nterm_pseu
                ptype = "NTERM"
        else:
            pseu = gap_pseu
            ptype = "GAP"
    return func, pseu, ptype


@logfun
def identity_writer(hitfile, outputdir, prefix, funcs, pseus, num_of_pseu):
    """
    Functions:
        identity_writer
        IdentityFunc.py module writer
    Parameter:
        hitfile, OR sequence annotated from genome
        outputdir, Result save directory
        prefix, output file prefix
        funcs, functions olfactory receptor list
        pseus, pseudo gene olfactory receptor list
        num_of_pseu: number of pseudogene every type
    Return: None
    """

    func_file = os.path.join(outputdir, prefix + "_redundant_func_ORs.fasta")
    pseu_file = os.path.join(outputdir, prefix + "_redundant_pseu_ORs.fasta")
    with open(func_file, 'w') as funcf:
        funcf.writelines(funcs)

    prepseu = hitfile[:-14] + "Pre-pseudos_dna.fa"
    if not os.path.exists(prepseu):
        raise FileNotExists(prepseu, outputdir)
    with open(prepseu) as prepseuf:
        prepseu_list = prepseuf.readlines()
    with open(pseu_file, 'w') as pseuf:
        prepseu_list.extend(pseus)
        pseuf.writelines(prepseu_list)
    func_OR_num = len(funcs)
    pseu_OR_num = len(prepseu_list)
    ntermpseu = num_of_pseu["NTERM"]
    gappseu = num_of_pseu["GAP"]
    funcdoc = "\033[1;32m{}\033[0m functional ORs were discover.".format(func_OR_num)
    pseudoc = "\033[1;32m{}\033[0m pseudogene ORs were discover.".format(pseu_OR_num)
    ntermdoc = ("\033[0;32m{}\033[0m pseudogenes cause by no N-term."
        .format(ntermpseu))
    gapdoc = ("\033[0;32m{}\033[0m pseudogenes cause by TM has gaps(more than {})."
        .format(gappseu, TM_GAPS_TOTAL))
    logging.info(funcdoc)
    logging.info(pseudoc)
    logging.info(ntermdoc)
    logging.info(gapdoc)
    return func_file, pseu_file, prepseu


def unredundant(outputdir, prefix, func_file, pseu_file, hitdna):
    """
    Function:
        unredundant
        sequence unredundance use CD-HIT
    Parameter:
        func_file, function olfactory file name.(Non redundant)
        pseu_file, pseudogene olfactory file name.(Non redundant)
        hitdna, Olfactory receptor DNA sequences file
    Return: None
    """

    nonredun_func = os.path.join(outputdir, prefix+"_func_ORs_pro.fasta")
    funcomm = (CDHIT
               + " -i " + func_file 
               + " -c 1.0 -T 0"
               + " -o " + nonredun_func
              )
    os.system(funcomm)

    nonredun_pseu = os.path.join(outputdir, prefix+"_pseu_ORs.fasta")
    pseucomm = (CDHIT 
                + " -i " + pseu_file 
                + " -c 1.0 -T 0"
                + " -o " + nonredun_pseu
               )
    os.system(pseucomm)

    # find and write functin ORs dna to file
    nonredun_dna = os.path.join(outputdir, prefix+"_func_ORs_dna.fasta")
    dnaseqs = ReadSampleFasta(hitdna)
    dnadict = dict(dnaseqs) 
    funcseqs = ReadSampleFasta(nonredun_func)
    funcdnas = ['>'+n+dnadict[n] for n, _ in funcseqs]
    with open(nonredun_dna, 'w') as outf:
        outf.writelines(funcdnas)


def merge_pseudo(pseudos):
    """
   Function:
       merge_pseudo
       merge pseudogene fragments to a complete sequences
       example: [100, 2000] + [1000, 2500] => [100, 2500]
       write sequence to file
    Parameter:
       pseudos: pseudogenes sequence dict.{name: seq}
    Return: None
    """
    pseudos = {'>'+k+'\n': v+'\n' for k, v in pseudos.items()}
    names = pseudos.keys()
    named = defaultdict(list)
    for name in names:
        ns = name.strip().split('_')
        head = ns[0]
        sign = ns[3]
        ptype = ns[4]
        if sign == '+':
            start = int(ns[1])
            end = int(ns[2])
        else:
            start = int(ns[2])
            end = int(ns[1])
        # save to a dict which from the same scaford
        named[head].append((start, end, name, sign, ptype))
    merge = {}
    for key in named:
        frags = named[key]
        length = len(frags)
        if length == 1:
            n = frags[0][2]
            merge[n] = pseudos[n]
            continue
        droup = []
        for i in range(length):
            s1, e1, n1, sign1, ptype1 = frags[i]
            n, seq = None, None
            for j in range(i+1, length):
                s2, e2, n2, sign2, ptype2 = frags[j]
                signx = '|'.join((sign1, sign2))
                ptypex = '|'.join((ptype1, ptype2))
                # [m, n] + [n-x, n+p] => [m, n+p]
                if (s1 < s2 < e1) and (e1 < e2):
                    seq = (
                        pseudos[n1][0:s2-s1-1]
                        + pseudos[n2][0:e2-s2] + '\n'
                    )
                    n = (
                        key + '_'
                        + str(s1) + '_'
                        + str(e2) + '_'
                        + signx + '_'
                        + ptypex + '\n'
                    )
                    droup.extend([n1, n2])
                # [m, n] + [m+x, n-p] => [m, n]
                elif (s1 <s2) and (e2 < e1):
                    seq = seq[n1]
                    n = n1
                    droup.extend([n1, n2])
                # [m-x, n-p] + [m, n] => [m-x, n]
                elif (s2 < s1 < e2) and (e2 < e1):
                    seq = (
                        pseudos[n2][0:s1-s2-1]
                        + pseudos[n1][0:e1-s1] + '\n'
                    )
                    n = (
                        key + '_'
                        + str(s2) + '_'
                        + str(e1) + '_'
                        + signx + '_'
                        + ptypex + '\n'
                    )
                    droup.extend([n1, n2])
                # [m+x, n-p] + [m, n] => [m, n]
                elif (s2 < s1) and (e1 < e2):
                    seq = pseudos[n2]
                    n = n2
                    droup.extend([n1, n2])
                else:
                    pass
            if seq:
                merge[n] = seq
            else:
                merge[n1] = pseudos[n1]
        for d in set(droup):
            try:
                del merge[d]
            except KeyError:
                pass
    # pseudos gene lines, for output
    outlines = ['>'+k+v for k, v in merge.items()]
    
    # number of every pseudos type
    type_of_pseu = {
        'INTER': 0,
        'SHORT': 0,
        'INDEL': 0,
        'INTER|INTER': 0,
        'INTER|SHORT': 0,
        'INTER|INDEL': 0,
        'SHORT|SHORT': 0,
        'SHORT|INDEL': 0,
        'SHORT|INTER': 0,
        'INDEL|INDEL': 0,
        'INDEL|SHORT': 0,
        'INDEL|INTER': 0,
    }
    ptyles = [k.split('_')[-1].strip() for k in merge.keys()]
    for pt in ptyles:
        type_of_pseu[pt] += 1
    pseu_type_num = {}
    pseu_type_num['INTER'] = (
        type_of_pseu['INTER']
        + type_of_pseu['INTER|INTER']
        + type_of_pseu['INTER|SHORT']
        + type_of_pseu['INTER|INDEL']
    )
    pseu_type_num['SHORT'] = (
        type_of_pseu['SHORT']
        + type_of_pseu['SHORT|SHORT']
        + type_of_pseu['SHORT|INDEL']
        + type_of_pseu['SHORT|INTER']
    )
    pseu_type_num['INDEL'] = (
        type_of_pseu['INDEL']
        + type_of_pseu['INDEL|INDEL']
        + type_of_pseu['INDEL|SHORT']
        + type_of_pseu['INDEL|INTER']
    )
     
    return outlines, pseu_type_num


def test_nhmmer_evalue(nhmmout):
    """
    Function:
        test_nhmmer_evalue
    Parameter:
        nhmmout: nhmmer tblout format file
    Return:
        return type -> string
        name, species name
        return type -> int
        nhits, number of nhmmer hits
    """

    name = nhmmout.split('-')[0]
    with open(nhmmout) as nhmmf:
        nhits = sum(1 for line in nhmmf if line[0] != '#')
    return name, nhits


def chmod(path):
    """
    Function:
        chmod
        Change file permission
    Parameter:
        path, file/folder path
    """
    if not os.path.exists(path):
        logging.error(path + " not exists!")
        sys.exit()
    elif os.path.isfile(path):
        os.chmod(path, stat.S_IRWXO + stat.S_IRWXG + stat.S_IRWXU)
    elif os.path.isdir(path):
        for _file in os.listdir(path):
            os.chmod(_file, stat.S_IRWXO + stat.S_IRWXG + stat.S_IRWXU)
    else:
        logging.info("chmod() nothing to do.")


if __name__ == "__main__":
    platform_info(True)
