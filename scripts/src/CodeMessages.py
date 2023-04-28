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

file_error = """
\033[1;31mFile Not ExistsError\033[0m
There is no {0} directory {1}.
Have you changed the {0} file name?
Or, you move it to other directory. 
"""

file_stype_error = """
\033[1;31mFile Type Error\033[0m
Input file this step is FindOR.py output file.
It's mean that {0} is 'Protein' sequence FASTA file,
and {1} is 'DNA' sequence FASTA file.
Place check!
"""

profile_error = """
\033[1;31mProfile Error\033[0m
The {0} option is incorrect. You can choose one of "Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Crocodylia", "Hyperoartia", "Lepidosauria", "Mammalia", "Myxini", "Reptiles", "Testudines", or provide the path to a better HMM profile file. However, we do not generally recommend doing so unless there is no corresponding option for the species you need to annotate in the list we provide. If you are sure that you have entered a file path, please check if the path is correct.
"""


class StrandError(Exception):

    def __init__(self, sing):
        self.sing = str(sing) + "\n"

    def __str__(self):
        return standmes.format(self.sing)


class LengthError(Exception):

    def __init__(self, length):
        self.length = str(length)+"\n"

    def __str__(self):
        return lengthmes.format(self.length)


class PlatformError(Exception):

    def __init__(self, system):
        self.system = str(system) + "\n"

    def __str__(self):
        return platformmes.format(self.system)


class VersionWarning(Exception):

    def __init__(self, version):
        self.version = str(version) + "\n"

    def __str__(self):
        return versionmes.format(self.version)


class FastaFormatError(Exception):

    def __init__(self, word):
        self.word = str(word) + "\n"

    def __str__(self):
        return fastaerror.format(self.word)


class NtermError(Exception):

    def __str__(self):
        return ntermerror


class FileNotExists(Exception):

    def __init__(self, filen, dirname):
        self.filen = filen
        self.dirname = dirname

    def __str__(self):
        if self.dirname is None:
            self.dirname = "Working"
        return file_error.format(self.filen, self.dirname)

class FileStypeError(Exception):

    def __init__(self, hitpro, hitdna):
        self.hitpro = hitpro
        self.hitdna = hitdna

    def __str__(self):
        return file_stype_error.format(self.hitpro, self.hitdna)

class ProfileError(Exception):

    def __init__(self, path):
        self.path = path

    def __str__(self):
        return profile_error.format(self.path)


