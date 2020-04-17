#!coding:utf-8


standmes = """
\033[1;31mSymbol Error! \033[0m
The NHMMER tool output 'strand' column only '+' or '-' symbol.
But '{0}' has exists.
"""

lengthmes = """
\033[1;31mShort Error! \033[0m
The length of CDS({0}) is not divisible by 3.
It cannot be translated into a correct protein sequence
"""

platformmes = """
\033[1;31mPlatform Error! \033[0m
This program only runs on the Linux platform.
Your platform: '{0}'
"""

versionmes = """
\033[0;31mVersion Warning! \033[0m
This program is best run on python3.
Your Python version: Python{0}
"""

fastaerror = """
\033[1;31mFasta Format Error\033[0m
Fasta format must start with '>'
Your file start '{0}'.
"""

ntermerror = """
\033[1;31mN-term length error!\033[0m
Your N-term list is empty. 
"""


class StrandError(Exception):

	def __init__(self, sing):
		self.sing = sing

	def __str__(self):
		return standmes.format(self.sing)


class LengthError(Exception):

	def __init__(self, length):
		self.length = length

	def __str__(self):
		return lengthmes.format(self.length)


class PlatformError(Exception):

	def __init__(self, system):
		self.system = system

	def __str__(self):
		return platformmes.format(self.system)


class VersionWarning(Exception):

	def __init__(self, version):
		self.version = str(version)

	def __str__(self):
		return versionmes.format(self.version)


class FastaFormatError(Exception):

	def __init__(self, word):
		self.word = word

	def __str__(self):
		return fastaerror.format(self.word)


class NtermError(Exception):

	def __str__(self):
		return ntermerror
