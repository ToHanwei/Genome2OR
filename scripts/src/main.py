#!coding:utf-8

import os
import sys
import logging

from src.config import *
from src.Functions import *
from src.data.data import *


class Nhmmer(object):
    """
    runing nhmmer program
    """
    
    def __init__(self, output, genome, profile, cpus, verbose, EvalueLimit):
        self.output = output
        self.gefile = genome
        self.pofile = profile
        self.cpus = cpus
        self.verbose = verbose
        self.EvalueLimit = EvalueLimit

    def check(self):
        """
        check the environment and variables
        """
        # open logging
        logger('nhmmer')
        
        # check platform
        platform_info(self.verbose)
        
        # change nhmmer permission
        chmod(NHMMER)
        
        # check CPU cores
        self.cpus = getCpus(self.cpus)
        self.cpus = str(self.cpus)
        
        # check profile path
        #self.pofile = getProfilePath(self.pofile)
        
        
    def run(self):
        """
        excute the command
        """
        sys.stdout.write(NHMMER_MESSAGE)
        if os.path.isdir(self.output):
            self.output = os.path.join(self.output, 'nhmmer_out.tblout')
            sys.stdout.write(
                "You have entered a directory. The output file will be redirected to {}.\n".format(self.output)
            )
        
        # Run the nhmmer program
        if self.verbose:
            sys.stdout.write(
                "\n\033[1;32mThe NHMMER program is running... \033[0m\n"
            )
        code, outext = run_nhmmer(
            self.EvalueLimit, 
            self.cpus,
            self.output, 
            self.pofile,
            self.gefile, 
            self.verbose,
        )
        
        if (code == 0) and self.verbose:
            sys.stdout.write(parseNhmmerOut(outext))
            sys.stdout.write(
                "\n\033[1;32mThe NHMMER program has finished running.\033[0m\n"
            )
        elif self.verbose:
            sys.stdout.write(
                "\n\033[1;31mThe NHMMER program ran incorrectly.\033[0m\n"
            )
        logging.info("###The nhmmer.py has completed.###")
        
        
        
class FindOR(object):
    """
    extract cds from genome
    """
    def __init__(self, input, genome, verbose, outputdir, prefix, EvalueLimit, SeqLengthLimit):
        self.infile = input
        self.gefile = genome
        self.verbose = verbose
        self.outputdir = outputdir
        self.prefix = prefix
        self.EvalueLimit = EvalueLimit
        self.SeqLengthLimit = SeqLengthLimit

    def check(self):
        """
        check the environment and variables
        """
        # open logging
        logger('FindOR')
        
        # check platform
        platform_info(self.verbose)
        
        # Check if the directory path exists,
        # and create the directory if it doesn't exist.
        existsDirectory(self.outputdir)
        
        
    def run(self):
        """
        excute the command
        """
        sys.stdout.write(FINDOR_MESSAGE)
        if self.verbose:
            sys.stdout.write(
                "\033[1;32mProcess nhmmer output file...\033[0m\n"
            )
        hmmout, trunc = proc_nhmmer_out(
            self.infile,
            self.EvalueLimit,
            self.SeqLengthLimit,
        )
        
        if self.verbose:
            sys.stdout.write(
                "\033[1;32mExtract cds from genomic file...\033[0m\n"
            )
        hmmout_seq, sour_seq = extract_cds(hmmout, self.gefile)
        
        if self.verbose:
            sys.stdout.write(
                "\033[1;32mFind ATG and STOP codons for each sequence...\033[0m\n"
            )
        functional, pseudos, low_quality = find_cds(
            hmmout,
            hmmout_seq,
            sour_seq, 
            self.SeqLengthLimit,
        )
        
        if self.verbose:
            sys.stdout.write(
                "\033[1;32mWrite data to file...\033[0m\n"
            )
        writer2file(
            self.outputdir, 
            self.prefix,
            functional,
            pseudos,
            hmmout_seq,
            hmmout,
            trunc,
            low_quality
        )
        
        logging.info("###The FindOR.py has completed.###")


class Identify(object):
    """
    identify olfactory receptor genes
    """

    def __init__(self, hitPROfile, hitDNAfile, cpus, outputdir, prefix, keepfile, verbose):
        self.hitpro = hitPROfile
        self.hitdna = hitDNAfile
        self.cpus = cpus
        self.outputdir = outputdir
        self.prefix = prefix
        self.keepfile = keepfile
        self.verbose = verbose

    def check(self):
        """
        check the environment and variables
        """
        # open logging
        logger('IdentityFunc')
        
        # check platform
        platform_info(self.verbose)
        
        # change tool permission
        chmod(MAFFT)
        chmod(CDHIT)
        
        # check CPU cores
        self.cpus = getCpus(self.cpus)
        self.cpus = str(self.cpus)
        
        # Check if the directory path exists,
        # and create the directory if it doesn't exist.
        existsDirectory(self.outputdir)
    

    def run(self):
        """
        excute the command
        """
        sys.stdout.write(IDENTIFYFUNC_MESSAGE)
        
        if self.verbose:
            sys.stdout.write(
                "\033[1;32mProcess hit sequence file...\033[0m\n"
            )
        hit_dict = refact_hitfile(self.hitpro)
        
        if self.verbose:
            sys.stdout.write(
                "\033[1;32mRecombination hit list...\033[0m\n"
            )
        #hit_list = refact_list(template, hit_dict, cpus)
        hit_list = refact_list(OR_TEMPLATE, hit_dict)
        
        if self.verbose:
            sys.stdout.write(
                "\033[1;32mIdentity functions\\pseudogenes...\033[0m\n"
            )
        pseus, funcs = [], []
        pasu_of_gap = 0
        
        pre_dnas_dict = dict(ReadSampleFasta(self.hitdna))
        
        for seq_list in hit_list:
            assert len(seq_list) > 1, \
            sys.stdout.write(
                "Please ensure that the input file here "
                + "is the output of the FindOR.py.\n"
            )
            # Prepare tm_list
            nterms, tm_list = tm_prepare(seq_list)
            seq_dict = dict(seq_list)
        
            # identity function ORs
            func, pseu, ptype = identify_filter(tm_list, nterms)
            if func:
                func_seq = (
                    ">" + func
                    + seq_dict[func].replace('-', '')
                )
                funcs.append(func_seq)
            elif pseu:
                pname = pseu.strip() + '_' + ptype + '\n'
                pseu_seq = (
                    ">" + pname
                    + pre_dnas_dict[pseu].replace('-', '')
                )
                pseus.append(pseu_seq)
                pasu_of_gap += 1
            else:
                sys.stdout.write(
                    "Function either or pseudogene\n"
                )
        
        if self.verbose:
            logging.info(
                "\033[1;32mWrite data to file...\033[0m"
            )
        
        # write functional ORs to file [Protein]
        profile = self.prefix + "_redundant_func_ORs.fasta" 
        Writer(self.outputdir, profile, funcs)
        
        # merge pseugenes
        pseufile = self.prefix + "_redundant_pseu_ORs.fasta"
        pseus = merge_pseudogene(self.hitpro, self.outputdir, pseus)
        Writer(self.outputdir, pseufile, pseus)
        
        # Unredundance sequences and write DNA sequences to file.
        if self.verbose:
            sys.stdout.write(
                "\033[1;32mThe CD-HIT program filters "
                + "redundant sequence...\033[0m\n"
            )
        
        # unredundant functional sequences
        # unredundant proteins
        profile = os.path.join(self.outputdir, profile)
        pro_outfile = self.prefix + "_final_func_pro_ORs.fasta"
        pro_outfile = os.path.join(self.outputdir, pro_outfile)
        message = "Remove redundant functional proteins"
        run_CDHIT(profile, pro_outfile, message, self.verbose)
        
        # unredundant DNAs
        dna_outfile = self.prefix + "_final_func_dna_ORs.fasta"
        dnaseqs = ReadSampleFasta(self.hitdna)
        dnadict = dict(dnaseqs) 
        funcseqs = ReadSampleFasta(pro_outfile)
        funcdnas = ['>'+n+dnadict[n] for n, _ in funcseqs]
        Writer(self.outputdir, dna_outfile, funcdnas)
        
        # unredundant pseugenes
        pseu_outfile = self.prefix + "_final_pseu_ORs.fasta"
        pseufile = os.path.join(self.outputdir, pseufile)
        pseu_outfile = os.path.join(self.outputdir, pseu_outfile)
        message = "Remove redundant pseugenes"
        run_CDHIT(pseufile, pseu_outfile, message, self.verbose)
        
        # unredundant low-quality sequences
        prelow = self.hitpro[:-14] + "Pre-low-quality_dna.fa"
        low_outfile = self.prefix +  "_final_low-quality.fasta"
        low_outfile = os.path.join(self.outputdir, low_outfile)
        if not os.path.exists(prelow):
            raise FileNotExists(prelow, None)
        message = "Remove redundant low-quality sequences"
        run_CDHIT(prelow, low_outfile, message, self.verbose)
        
        
        # print summary
        truncfile = self.hitpro[:-14] + "truncated.txt"
        printSummary(pseu_outfile, pro_outfile, low_outfile, truncfile)
        
        
        # delect file, if keepfile=False
        self.keepfile = boolean_string(self.keepfile)
        if not self.keepfile:
            os.remove(self.hitpro)
            logging.info(
                "Delete the {} file.".format(self.hitpro)
            )
            os.remove(self.hitdna)
            logging.info(
                "Delete the {} file.".format(self.hitdna)
            )
            os.remove(profile)
            logging.info(
                "Delete the {} file.".format(profile)
            )
            os.remove(pseufile)
            logging.info(
                "Delete the {} file.".format(pseufile)
            )
            os.remove(prelow)
            logging.info(
                "Delete the {} file.".format(prelow)
            )
            del_path = os.path.dirname(
                os.path.abspath(self.hitpro)
            )
            removeSuffixFile(del_path, ".clstr")
            prepseu = self.hitpro[:-14] + "Pre-pseudos_dna.fa"
            os.remove(prepseu)
            logging.info(
                "Delete the {} file.".format(prepseu)
            )
        
        logging.info("###The IdentityFunc.py has completed.###")
