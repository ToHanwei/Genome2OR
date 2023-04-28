# Genome2OR
Annotate Olfactory receptor CDS from genome
==============================================


## Contents Index
* [Cite](#cite)
* [Abstract](#abstract)
* [Installation](#installation)
* [Contact us](#contact-us)

## Cite
* Cite us </br>
[Welcome to our paper](https://link.springer.com/article/10.1007/s11427-021-2081-6). </br>
Han W, Wu Y, Zeng L, Zhao S. Building the Chordata Olfactory Receptor Database using more than 400,000 receptors annotated by Genome2OR. Sci China Life Sci. 2022 Dec;65(12):2539-2551. doi: 10.1007/s11427-021-2081-6. Epub 2022 Jun 10. PMID: 35696018. </br>
[Article](./data/s11427-021-2081-6.pdf) </br>
[Support information](./data/11427_2021_2081_MOESM1_ESM.pdf) </br>

* Resource </br>
Welcome to "[Chordata Olfactory Receptor Database(CORD)](https://cord.ihuman.shanghaitech.edu.cn/#/home)". You will be able to obtain a large amount of data on olfactory receptors in chordata.

## Abstract
[Genome2OR](https://github.com/ToHanwei/Genome2OR.git) is a genetic annotation tool based on [HMMER](http://hmmer.org), [MAFFT](https://mafft.cbrc.jp/alignment/software/) and [CD-HIT](http://weizhongli-lab.org/cd-hit/).  
[HMMER](http://hmmer.org) searches biological sequence databases for homologous sequences, using either single sequences or multiple sequence alignments as queries. HMMER implements a technology called "profile hidden Markov models (profile HMMs)".  

[MAFFT](https://mafft.cbrc.jp/alignment/software/) is a multiple sequence alignment program for unix-like operating systems. It offers a range of multiple alignment methods, L-INS-i (accurate; for alignment of <∼200 sequences),  FFT-NS-2 (fast; for alignment of <∼30,000 sequences), etc.

[CD-HIT](http://weizhongli-lab.org/cd-hit/) is very fast and can handle extremely large databases. CD-HIT helps to significantly reduce the computational and manual efforts in many sequence analysis tasks and aids in understanding the data structure and correct the bias within a dataset.

## Installation

Use anaconda for installation.

```bash
conda install -c whanwei genome2or
```

After installation you can use the following two commands to annotate the olfactory receptor genes in the genome.

### Basic annotation method

```bash
genome2or Actinopteri outputdir genome.fasta -e 1e-10 -c 4 -p prefix
```

For more information.

```bash
genome2or --help

    ____ _____ _   _  ___  __  __ _____ ____   ___  ____
   / ___| ____| \ | |/ _ \|  \/  | ____|___ \ / _ \|  _ \
  | |  _|  _| |  \| | | | | |\/| |  _|   __) | | | | |_) |
  | |_| | |___| |\  | |_| | |  | | |___ / __/| |_| |  _ <
   \____|_____|_| \_|\___/|_|  |_|_____|_____|\___/|_| \_|


usage: genome2or.py [-h] [-e] [-l] [-c] [-p] [-k] [-v] [-V] profile outputdir genome

Annotating Olfactory Receptor Genes in Vertebrate Genomes in One Step.

positional arguments:
  profile               Select an HMM profile for the annotated species from the following options: 		   						 "Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", "Chondrichthyes", 							"Cladistia", "Coelacanthimorpha", "Crocodylia", "Hyperoartia", "Lepidosauria", 								"Mammalia", "Myxini", "Reptiles", "Testudines".
                        Alternativa, provide the path to a HMM profile file. However, we do not generally 							recommend doing so unless there is no corresponding option for the species you need 						to annotate in the list we provide.
  outputdir             String. Directory path where the output files are stored.[default:Current directory]
  genome                String. File path of the genome to be annotated.

optional arguments:
  -h, --help            show this help message and exit
  -e , --EvalueLimit    Float, The e-value threashold used for extract olfactory receptor gene fragment(s) from the genome.[default:1e-20]
  -l , --SeqLengthLimit
                        Integer. Threshold of sequence length, sequences shoter than this value will not be considered as the preferred targets for functional olfactory
                        receptors.[default:868]
  -c , --cpus           Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']
  -p , --prefix         String. Output file name prefix.[default:ORannotation]
  -k , --keepfile       Bool. whether to keep detailed intermediate file(True/False).[default:True]
  -v, --verbose         Print detailed running messages of the program.
  -V, --version         Show version message and exit.

https://link.springer.com/article/10.1007/s11427-021-2081-6
```

### Iterative annotation of the genome

```bash
Iteration Actinopteri outputdir genome.fasta -i 3 -e 1e-10 -c 4 -p prefix
```

For more information.

```bash
Iteration --help

     ___ _____ _____ ____      _  _____ ___ ___  _   _
    |_ _|_   _| ____|  _ \    / \|_   _|_ _/ _ \| \ | |
     | |  | | |  _| | |_) |  / _ \ | |  | | | | |  \| |
     | |  | | | |___|  _ <  / ___ \| |  | | |_| | |\  |
    |___| |_| |_____|_| \_\/_/   \_\_| |___\___/|_| \_|


usage: Iteration.py [-h] [-i] [-e] [-l] [-c] [-p] [-k] [-v] [-V] profile outputdir genome

Iterative annotation of olfactory receptor genes in the genome.

positional arguments:
  profile               Select an HMM profile for the annotated species from the following options: 								"Actinopteri", "Amphibia", "Aves", "Branchiostoma_floridae", "Chondrichthyes", 								"Cladistia", "Coelacanthimorpha", "Crocodylia", "Hyperoartia", "Lepidosauria", 								"Mammalia", "Myxini", "Reptiles", "Testudines".
                        Alternativa, provide the path to a HMM profile file. However, we do not generally 							recommend doing so unless there is no corresponding option for the species you need 						to annotate in the list we provide.
  outputdir             String. Directory path where the output files are stored.[default:Current directory]
  genome                String. File path of the genome to be annotated.

optional arguments:
  -h, --help            show this help message and exit
  -i , --iteration      Int. Number of iterations.[default:2]
  -e , --EvalueLimit    Float, The e-value threashold used for extract olfactory receptor gene fragment(s) from the genome.[default:1e-20]
  -l , --SeqLengthLimit
                        Integer. Threshold of sequence length, sequences shoter than this value will not be considered as the preferred targets for functional olfactory
                        receptors.[default:868]
  -c , --cpus           Integer. number of parallel, with 0, all CPUs will be used.[default='2/3 of all cores']
  -p , --prefix         String. Output file name prefix.[default:ORannotation]
  -k , --keepfile       Bool. whether to keep detailed intermediate file(True/False).[default:True]
  -v, --verbose         Print detailed running messages of the program.
  -V, --version         Show version message and exit.

https://link.springer.com/article/10.1007/s11427-021-2081-6
```

### Batch annotation

For batch annotation, you can use the following simple shell script to achieve it.

```bash
GenomeDir="Path to the directory where you store your genome"
for genome in `ls $GenomeDir`; do
	genome2or Actinopteri outputdir $genome -e 1e-10 -c 4 -p ${genome%.*}
done
```

Here we assume that you are annotating species of the "Actinopteri" and that the genomes stored in your catalog all belong to the "Actinopteri". For annotation of species of other orders, please select the corresponding HMM profile with the "profile" parameter.

Or you want to use iterations for batch annotation.

```bash
GenomeDir="Path to the directory where you store your genome"
for genome in `ls $GenomeDir`; do
	Iteration Actinopteri outputdir $genome -i 3 -e 1e-10 -c 4 -p ${genome%.*}
done
```



## Contact us

* zhaosw@shanghaitech.edu.cn, Suwen Zhao
* hanwei@shanghaitech.edu.cn, Wei Han
