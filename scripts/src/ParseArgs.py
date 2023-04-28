#!coding:utf-8

import os
import sys
import argparse

from multiprocessing import cpu_count
from src.config import TOOLNAME
from src.config import VERSION
from src.Functions import boolean_string


class ParseCommand(object):
    """
    Parse command line
    """

    def parsecommand(self):
        """
        FindOR.py command line parse
        """
        #sys.stdout.write(FINDOR_MESSAGE)

        parser = argparse.ArgumentParser(
            prog=TOOLNAME,
            description="Mining for olfactory receptor genes from genome.",
            epilog="https://link.springer.com/article/10.1007/s11427-021-2081-6"
        )
        # positional arguments
        parser.add_argument(
            "input",
            type=str,
            help="String. Input file, which is also the output file of the nhmmer.py program.",
        )
        parser.add_argument(
            'genome',
            type=str,
            help="String. File path of the genome to be annotated.",
        )

        # optional arguments
        parser.add_argument(
            '-o', '--outputdir',
            type=str,
            help="String. Directory path where the output files are stored.[default:Current directory]",
            default=os.getcwd(),
            metavar='',
        )
        parser.add_argument(
            '-p', '--prefix',
            type=str,
            help="String. Output file name prefix.[default:ORannotation]",
            default="ORannotation",
            metavar='',
        )
        parser.add_argument(
            '-e', '--EvalueLimit',
            type=float,
            help="Float, The e-value threashold used for extract olfactory receptor gene fragment(s) from the genome.[default:1e-20]",
            # An empirical value
            default=1e-20,
            metavar='',
        )
        parser.add_argument(
            '-l', '--SeqLengthLimit',
            type=int,
            help="Integer. Threshold of sequence length, sequences shoter than this value will not be considered as the preferred targets for functional olfactory receptors.[default:868]",
            # The OR sequence length is approximately 310
            default=868,
            metavar='',
        )
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            help="Print detailed running messages of the program.",
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            help='Show version message and exit.',
            version=VERSION,
        )
        args = parser.parse_args()
        return args

    def nhmmerparse(self):
        """nhmmer.py command line parse"""

        #sys.stdout.write(NHMMER_MESSAGE)

        parser = argparse.ArgumentParser(
            prog="Run nhmmer",
            description="Identify the nucleotide fragments of olfactory receptor genes from the genome.",
            epilog="https://link.springer.com/article/10.1007/s11427-021-2081-6"
        )
        # positional arguments
        parser.add_argument(
            'profile',
            help="""
                Select an HMM profile for the annotated species from the following options: 
                "Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", 
                "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Crocodylia", 
                "Hyperoartia", "Lepidosauria", "Mammalia", "Myxini", "Reptiles", "Testudines".
                Alternativa, provide the path to a  HMM profile file. However, 
                we do not generally recommend doing so unless there is no corresponding 
                option for the species you need to annotate in the list we provide.
                """
        )
        parser.add_argument(
            'genome',
            help="String. File path of the genome to be annotated.",
        )
        parser.add_argument(
            "output",
            help="String. Output file path.",
        )

        # optional arguments
        parser.add_argument(
            '-e', '--EvalueLimit',
            type=float,
            help="Float. The e-value threshold for searching homologous sequences of odorant receptors in the genome.[default:1e-20]",
            # An empirical value
            default=1e-20,
            metavar='',
        )
        parser.add_argument(
            '-c', '--cpus',
            type=int,
            help="Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']",
            default=int(cpu_count() * 2 / 3),
            metavar='',
        )
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            help="Print detailed running messages of the program.",
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            help='Show version message and exit.',
            version=VERSION,
        )
        args = parser.parse_args()
        return args

    def batchparse(self):
        """
        batch.py command line parse
        """

        #sys.stdout.write(BATCH_MESSAGE)

        parser = argparse.ArgumentParser(
            prog="batch.py",
            description="Batch annotation of olfactory receptor genes in the genome.",
            epilog="https://link.springer.com/article/10.1007/s11427-021-2081-6"
        )
        # positional arguments
        parser.add_argument(
            'profile',
            help="""
                Select an HMM profile for the annotated species from the following options: 
                "Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", 
                "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Crocodylia", 
                "Hyperoartia", "Lepidosauria", "Mammalia", "Myxini", "Reptiles", "Testudines".
                Alternativa, provide the path to a  HMM profile file. However, 
                we do not generally recommend doing so unless there is no corresponding 
                option for the species you need to annotate in the list we provide.
                """
        )
        parser.add_argument(
            'inputdir',
            help="""
                String. Enter a folder path for several genomes to be annotated.
                Please note that the species in the same input directory should 
                be of the same type (belonging to the same class or closer taxa), 
                as the species in the input directory share the same HMM profile.
            """,
        )
        parser.add_argument(
            "nhmmerout",
            help="String. Output folder path for the output of the nhmmer program.",
        )
        parser.add_argument(
            "findoutdir",
            help="String. Output folder path for the output of the FindOR program.",
        )
        parser.add_argument(
            "identifyout",
            help="String. Output folder path for the output of the Identifyout program."
        )

        # optional arguments
        parser.add_argument(
            '-e', '--EvalueLimit',
            type=float,
            help="Float, The e-value threashold used for extract olfactory receptor gene fragment(s) from the genome.[default:1e-20]",
            # An empirical value
            default=1e-20,
            metavar='',
        )
        parser.add_argument(
            '-l', '--SeqLengthLimit',
            type=int,
            help="Integer. Threshold of sequence length, sequences shoter than this value will not be considered as the preferred targets for functional olfactory receptors.[default:868]",
            # The OR sequence length is approximately 310
            default=869,
            metavar='',
        )
        parser.add_argument(
            '-c', '--cpus',
            type=int,
            help="Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']",
            default=int(cpu_count() * 2 / 3),
            metavar='',
        )
        parser.add_argument(
                '-k', '--keepfile',
                help="Bool. whether to keep detailed intermediate file(True/False).[default:True]",
                default=True,
                metavar='',
                type=boolean_string,
        )
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            help="Print detailed running messages of the program.",
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            help='Show version message and exit.',
            version=VERSION,
        )
        args = parser.parse_args()
        return args


    def Iteraparse(self):
        """Iteration.py command line parse"""

        #sys.stdout.write(ITERATION_MESSAGE)

        parser = argparse.ArgumentParser(
            prog="Iteration.py",
            description="Iterative annotation of olfactory receptor genes in the genome.",
            epilog="https://link.springer.com/article/10.1007/s11427-021-2081-6"
        )
        # positional arguments
        parser.add_argument(
            'profile',
            help="""
                Select an HMM profile for the annotated species from the following options: 
                "Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", 
                "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Crocodylia", 
                "Hyperoartia", "Lepidosauria", "Mammalia", "Myxini", "Reptiles", "Testudines".
                Alternativa, provide the path to a  HMM profile file. However, 
                we do not generally recommend doing so unless there is no corresponding 
                option for the species you need to annotate in the list we provide.
                """
        )
        parser.add_argument(
            'outputdir',
            help="String. Directory path where the output files are stored.[default:Current directory]",
        )
        parser.add_argument(
            'genome',
            help="String. File path of the genome to be annotated.",
        )

        # optional arguments
        parser.add_argument(
            '-i', '--iteration',
            help="Int. Number of iterations.[default:2]",
            type=int,
            default=2,
            metavar=''
        )
        parser.add_argument(
            '-e', '--EvalueLimit',
            type=float,
            help="Float, The e-value threashold used for extract olfactory receptor gene fragment(s) from the genome.[default:1e-20]",
            # An empirical value
            default=1e-20,
            metavar='',
        )
        parser.add_argument(
            '-l', '--SeqLengthLimit',
            type=int,
            help="Integer. Threshold of sequence length, sequences shoter than this value will not be considered as the preferred targets for functional olfactory receptors.[default:868]",
            # The OR sequence length is approximately 310
            default=869,
            metavar='',
        )
        parser.add_argument(
            '-c', '--cpus',
            type=int,
            help="Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']",
            default=int(cpu_count() * 2 / 3),
            metavar='',
        )
        parser.add_argument(
            '-p', '--prefix',
            help="String. Output file name prefix.[default:ORannotation]",
            default="ORannotation",
            metavar='',
        )
        parser.add_argument(
                '-k', '--keepfile',
                help="Bool. whether to keep detailed intermediate file(True/False).[default:True]",
                default=True,
                metavar='',
                type=boolean_string,
        )
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            help="Print detailed running messages of the program.",
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            help='Show version message and exit.',
            version=VERSION,
        )
        args = parser.parse_args()
        return args


    def genome2or(self):
        """genome2or.py command line parse"""

        #sys.stdout.write(GENOME2OR_MESSAGE)

        parser = argparse.ArgumentParser(
            prog="genome2or.py",
            description="Annotating Olfactory Receptor Genes in Vertebrate Genomes in One Step.",
            epilog="https://link.springer.com/article/10.1007/s11427-021-2081-6"
        )
        # positional arguments
        parser.add_argument(
            'profile',
            help="""
                Select an HMM profile for the annotated species from the following options: 
                "Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", 
                "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Crocodylia", 
                "Hyperoartia", "Lepidosauria", "Mammalia", "Myxini", "Reptiles", "Testudines".
                Alternativa, provide the path to a  HMM profile file. However, 
                we do not generally recommend doing so unless there is no corresponding 
                option for the species you need to annotate in the list we provide.
                """
        )
        parser.add_argument(
            'outputdir',
            help="String. Directory path where the output files are stored.[default:Current directory]",
        )
        parser.add_argument(
            'genome',
            help="String. File path of the genome to be annotated.",
        )

        # optional arguments
        parser.add_argument(
            '-e', '--EvalueLimit',
            type=float,
            help="Float, The e-value threashold used for extract olfactory receptor gene fragment(s) from the genome.[default:1e-20]",
            # An empirical value
            default=1e-20,
            metavar='',
        )
        parser.add_argument(
            '-l', '--SeqLengthLimit',
            type=int,
            help="Integer. Threshold of sequence length, sequences shoter than this value will not be considered as the preferred targets for functional olfactory receptors.[default:868]",
            # The OR sequence length is approximately 310
            default=869,
            metavar='',
        )
        parser.add_argument(
            '-c', '--cpus',
            type=int,
            help="Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']",
            default=int(cpu_count() * 2 / 3),
            metavar='',
        )
        parser.add_argument(
            '-p', '--prefix',
            help="String. Output file name prefix.[default:ORannotation]",
            default="ORannotation",
            metavar='',
        )
        parser.add_argument(
                '-k', '--keepfile',
                help="Bool. whether to keep detailed intermediate file(True/False).[default:True]",
                default=True,
                metavar='',
                type=boolean_string,
        )
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            help="Print detailed running messages of the program.",
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            help='Show version message and exit.',
            version=VERSION,
        )
        args = parser.parse_args()
        return args


    def IteraBatchparse(self):
        """IterationBatch.py command line parse"""

        #sys.stdout.write(ITERBATH_MESSAGE)

        parser = argparse.ArgumentParser(
            prog="IterationBatch.py",
            description="Iteration annotated batch genome",
            epilog="https://link.springer.com/article/10.1007/s11427-021-2081-6"
        )
        # positional arguments
        parser.add_argument(
            'profile',
            help="""
                Select an HMM profile for the annotated species from the following options: 
                "Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", 
                "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Crocodylia", 
                "Hyperoartia", "Lepidosauria", "Mammalia", "Myxini", "Reptiles", "Testudines".
                Alternativa, provide the path to a  HMM profile file. However, 
                we do not generally recommend doing so unless there is no corresponding 
                option for the species you need to annotate in the list we provide.
                """
        )
        parser.add_argument(
            'outputdir',
            help="String. Directory path where the output files are stored.[default:Current directory]",
        )
        parser.add_argument(
            'genomedir',
            help="String. Enter a folder path for several genomes to be annotated.",
        )

        # optional arguments
        parser.add_argument(
            '-i', '--iteration',
            help="Int, Number of iterations.[default:2]",
            type=int,
            default=2,
            metavar=''
        )
        parser.add_argument(
            '-e', '--EvalueLimit',
            type=float,
            help="Float, The e-value threashold used for extract olfactory receptor gene fragment(s) from the genome.[default:1e-20]",
            # An empirical value
            default=1e-20,
            metavar='',
        )
        parser.add_argument(
            '-l', '--SeqLengthLimit',
            type=int,
            help="Integer. Threshold of sequence length, sequences shoter than this value will not be considered as the preferred targets for functional olfactory receptors.[default:868]",
            # The OR sequence length is approximately 310
            default=869,
            metavar='',
        )
        parser.add_argument(
            '-c', '--cpus',
            type=int,
            help="Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']",
            metavar='',
        )
        parser.add_argument(
                '-k', '--keepfile',
                help="Bool. whether to keep detailed intermediate file(True/False).[default:True]",
                default=True,
                metavar='',
                type=boolean_string,
        )
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            help="Print detailed running messages of the program.",
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            help='Show version message and exit.',
            version=VERSION,
        )
        args = parser.parse_args()
        return args


    def identifyparse(self):
        """
        IdentifyFunc.py command line parse
        """

        #sys.stdout.write(IDENTIFYFUNC_MESSAGE)

        parser = argparse.ArgumentParser(
            prog="IdentifyFunc.py",
            description="Identifying olfactory receptor gene(s) in the genome",
            epilog="https://link.springer.com/article/10.1007/s11427-021-2081-6"
        )
        # positional arguments
        parser.add_argument(
            'hitPROfile',
            help="IdentifyFunc.py script output file(protein sequence file). \
            hit sequence from genome"
        )
        parser.add_argument(
            'hitDNAfile',
            help="IdentifyFunc.py script output file(DNA sequence file). \
            hit sequence from genome"
        )

        # optional arguments
        parser.add_argument(
            '-o', '--outputdir',
            help="String. Directory path where the output files are stored.[default:Current directory]",
            default=".",
            metavar='',
        )
        parser.add_argument(
            '-p', '--prefix',
            help="String. Output file name prefix.[default:ORannotation]",
            default="ORannotation",
            metavar='',
        )
        parser.add_argument(
            '-c', '--cpus',
            type=int,
            help="Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']",
            default=int(cpu_count() * 2 / 3),
            metavar='',
        )
        parser.add_argument(
                '-k', '--keepfile',
                help="Bool. whether to keep detailed intermediate file(True/False).[default:True]",
                default=True,
                metavar='',
                type=boolean_string,
        )
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            help="Print detailed running messages of the program.",
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            help='Show version message and exit.',
            version=VERSION,
        )
        args = parser.parse_args()
        return args

    def testparse(self):
        """nhmmer_evalue_test.py command line parse"""

        parser = argparse.ArgumentParser(
            prog="nhmmer_evalue_test.py",
            description="Test the number of nhmmer hits",
            epilog="https://link.springer.com/article/10.1007/s11427-021-2081-6"
        )
        # positional arguments
        parser.add_argument(
            'profile',
            help="""
                Select an HMM profile for the annotated species from the following options: 
                "Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", 
                "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Crocodylia", 
                "Hyperoartia", "Lepidosauria", "Mammalia", "Myxini", "Reptiles", "Testudines".
                Alternativa, provide the path to a  HMM profile file. However, 
                we do not generally recommend doing so unless there is no corresponding 
                option for the species you need to annotate in the list we provide.
                """
        )
        parser.add_argument(
            'inputdir',
            help="""
                String. Enter a folder path for several genomes to be annotated.
                Please note that the species in the same input directory should 
                be of the same type (belonging to the same class or closer taxa), 
                as the species in the input directory share the same HMM profile.
            """,
        )
        parser.add_argument(
            "nhmmerout",
            help="String. Output folder path for the output of the nhmmer program.",
        )

        # optional arguments
        parser.add_argument(
            '-e', '--EvalueLimit',
            type=float,
            help="Float, The e-value threashold used for extract olfactory receptor gene fragment(s) from the genome.[default:1e-20]",
            # An empirical value
            default=1e-20,
            metavar='',
        )
        parser.add_argument(
            '-o', '--output',
            type=str,
            help="String. Output file path.",
            # The OR sequence length is approximately 310
            default="nhmmer_hits_test.csv",
            metavar='',
        )
        parser.add_argument(
            '-c', '--cpus',
            type=int,
            help="Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']",
            default=int(cpu_count() * 2 / 3),
            metavar='',
        )
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            help="Print detailed running messages of the program.",
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            help='Show version message and exit.',
            version=VERSION,
        )
        args = parser.parse_args()
        return args


if __name__ == "__main__":
    parse = ParseCommand()
    args = parse.nhmmerparse()
    print(args)
