#!coding:utf-8

import argparse
from multiprocessing import cpu_count

from src.config import TOOLNAME
from src.config import VERSION


class ParseCommand(object):
	"""Parse command line"""

	def parsecommand(self):
		"""FindOR.py command line parse"""

		parser = argparse.ArgumentParser(prog=TOOLNAME,
										 description="Olfactory receptor annotation",
										 epilog="http://zhaolab.shanghaitech.edu.cn/"
										 )
		# positional arguments
		parser.add_argument("input",
							help="String, nhmmer output file path.",
							)
		parser.add_argument('genome',
							help="String, Genomic data file path.",
							)

		# optional arguments
		parser.add_argument('-o', '--outputdir',
							help="String, Result save directory.(default:../output)",
							default="../output",
							metavar='',
							)
		parser.add_argument('-p', '--prefix',
							help="String, output file prefix.(default:ORannotation)",
							default="ORannotation",
							metavar='',
							)
		parser.add_argument('-e', '--EvalueLimit',
							type=float,
							help="Float, Sequence similarity threshold.\
								 (default:1e-60)",
							# An empirical value
							default=1e-60,
							metavar='',
							)
		parser.add_argument('-l', '--SeqLengthLimit',
							type=int,
							help="Int, An artificially set OR's sequence \
								 length threshold.(default:868)",
							# The OR sequence length is approximately 310
							default=869,
							metavar='',
							)
		parser.add_argument('-v', '--verbose',
							action='count',
							help="Print verbose information.",
							)
		parser.add_argument('-V', '--version',
							action='version',
							help='Show version message and exit.',
							version=VERSION,
							)
		args = parser.parse_args()
		return args

	def nhmmerparse(self):
		"""nhmmer.py command line parse"""

		parser = argparse.ArgumentParser(prog="run_nhmmer",
										 description="Autorun nhmmer",
										 epilog="http://zhaolab.shanghaitech.edu.cn/"
										 )
		# positional arguments
		parser.add_argument('profile',
							help="String, profile nhmmer need(hmm, [un]alignment file)"
							)
		parser.add_argument('genome',
							help="String, genomic data file path.",
							)
		parser.add_argument("output",
							help="String, save nhmmer output file path.",
							)

		# optional arguments
		parser.add_argument('-e', '--EvalueLimit',
							type=float,
							help="Float, Sequence similarity threshold.\
								 (default:1e-10)",
							# An empirical value
							default=1e-10,
							metavar='',
							)
		parser.add_argument('-c', '--cpus',
							type=int,
							help="number of parallel CPU workers to use. \
								 (default='2/3 of all cores')",
							default=int(cpu_count() * 2 / 3),
							metavar='',
							)
		parser.add_argument('-v', '--verbose',
							action='count',
							help="Print verbose information.",
							)
		parser.add_argument('-V', '--version',
							action='version',
							help='Show version message and exit.',
							version=VERSION,
							)
		args = parser.parse_args()
		return args

	def batchparse(self):
		"""batch.py command line parse"""

		parser = argparse.ArgumentParser(prog="BatchProcess.py",
										 description="Batch annotated genome",
										 epilog="http://zhaolab.shanghaitech.edu.cn/"
										 )
		# positional arguments
		parser.add_argument('profile',
							help="String, profile nhmmer need(hmm, [un]alignment file).\
							Notice: A group of genomes share a profile."
							)
		parser.add_argument('inputdir',
							help="String, genomic directory.",
							)
		parser.add_argument("nhmmerout",
							help="String, run nhmmer program result.",
							)
		parser.add_argument("findoutdir",
							help="String, processing results directory, \
								 run FindOR.py result.",
							)
		parser.add_argument("identityout",
							help="String, processing results directory, \
								 run IdentityFunc.py result.",
							)

		# optional arguments
		parser.add_argument('-e', '--EvalueLimit',
							type=float,
							help="Float, Sequence similarity threshold.\
								 (default:1e-60)",
							# An empirical value
							default=1e-60,
							metavar='',
							)
		parser.add_argument('-l', '--SeqLengthLimit',
							type=int,
							help="Int, An artificially set OR's sequence \
								 length threshold.(default:868)",
							# The OR sequence length is approximately 310
							default=869,
							metavar='',
							)
		parser.add_argument('-c', '--cpus',
							type=int,
							help="number of parallel CPU workers to use. \
								 (default='2/3 of all cores')",
							default=int(cpu_count() * 2 / 3),
							metavar='',
							)
		parser.add_argument('-v', '--verbose',
							action='count',
							help="Print verbose information.",
							)
		parser.add_argument('-V', '--version',
							action='version',
							help='Show version message and exit.',
							version=VERSION,
							)
		args = parser.parse_args()
		return args

	def identifyparse(self):
		"""IdentityFunc.py command line parse"""

		parser = argparse.ArgumentParser(prog="IdentityFunc.py",
										 description="Idntity Function OR",
										 epilog="http://zhaolab.shanghaitech.edu.cn/"
										 )
		# positional arguments
		parser.add_argument('hitfile',
							help="IdentityFunc.py script output file(protein sequence file). \
								 hit sequence from genome"
							)

		# optional arguments
		parser.add_argument('-o', '--outputdir',
							help="String, Result save directory.(default:../output)",
							default="../output",
							metavar='',
							)
		parser.add_argument('-p', '--prefix',
							help="String, output file prefix.(default:Identity)",
							default="Identity",
							metavar='',
							)
		parser.add_argument('-v', '--verbose',
							action='count',
							help="Print verbose information.",
							)
		parser.add_argument('-V', '--version',
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

