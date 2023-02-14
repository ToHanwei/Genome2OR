# Genome2OR
Annotate Olfactory receptor CDS from genome
==============================================


## Contents Index
* [Cite](#cite)
* [Abstract](#abstract)
* [Quick start](#quick-start)
* [Download script](#download-script)
* [Usage](#usage)
* [Example](#example)
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
[Genome2OR](https://github.com/ToHanwei/Genome2OR.git) is a genetic annotation tool based on 
[HMMER](http://hmmer.org), [MAFFT](https://mafft.cbrc.jp/alignment/software/) and 
[CD-HIT](http://weizhongli-lab.org/cd-hit/).  
[HMMER](http://hmmer.org) searches biological sequence databases for
homologous sequences, using either single sequences or multiple
sequence alignments as queries. HMMER implements a technology called
"profile hidden Markov models (profile HMMs)".  

[MAFFT](https://mafft.cbrc.jp/alignment/software/) is a multiple sequence alignment program 
for unix-like operating systems. 
It offers a range of multiple alignment methods, L-INS-i (accurate; for alignment of <∼200 sequences), 
FFT-NS-2 (fast; for alignment of <∼30,000 sequences), etc.

[CD-HIT](http://weizhongli-lab.org/cd-hit/) is very fast and can handle extremely large databases. 
CD-HIT helps to significantly reduce the computational and manual efforts in many sequence analysis 
tasks and aids in understanding the data structure and correct the bias within a dataset.

## Quick start
1. STEP 1: Execute nhmmer.py
    ```
   cd YouDir/Genome2OR/scripts
   python nhmmer.py ../template/Mammalia.hmm genome.fasta nhmmer_out.tblout -v
    ```
    * The "Mammalia.hmm" file is the "HMM profile" for mammalia. This file provided by us. Of course, you can also use the NHMMER tool to construct it yourself according to your needs. Generally speaking, we do not recommend doing this unless you konw why you need to construct this file anew.
    * You can find more species of "HMM profile" in "the directory YouDir/Genome2OR/template" if needed. To ensure the quality of annotation, please select the appropriate "HMM profile" for the species you need to annotate.
    * The "genome.fasta" file is the genome file that you need to  annotate. This file is in the FASTA file format.
    * The "nhmmer_out.tblout" file is the output file. You will find it in the current working directory.

2. STEP 2: Execute FindOR.py
    ```
   python FindOR.py nhmmer_out.tblout genome.fasta -o ../output -v 
    ```
    * The "nhmmer_out.tblout" file here is the output file of the "STEP 1".
    * The "genome.fasta" file is the genome file that you need to  annotate. This file is the same as the file with the same name in the "STEP 1".
    * The "../output" is the directory where the output files are saved.
    * After the program is finished running, you will find five output files with the prefix "ORannotation" in the "../output" directory.
   
3. STEP 3: Execute IdentityFunc.py
    ```
   python IdentifyFunc.py ../output/ORannotation_Pre-ORs_pro.fa ../output/ORannotation_Pre-ORs_dna.fa -o ../output -p Identity -v
    ```
    * The "../output/ORannotation_Pre-ORs_pro.fa" file here is the output file of the "STEP 2".
    * The "../output/ORannotation_Pre-ORs_dna.fa" file here is the output file of the "STEP 2".
    * The "../output" is the directory where the output files are saved.
    * After the program is finished running, you will find seven output files with the prefix "Identity" in the "../output" directory.
    * In the output directory, the file "Identity_redundant_func_ORs.fasta" contains annotated functional olfactory receptor protein sequences, the file "Identity_redundant_pseu_ORs.fasta" contatins annotated pseudogene sequences, and the file "ORannotation_truncated.txt" records truncated genes.


## Download script
1. Download script from github
    ```
    cd YourDir
    git clone https://github.com/ToHanwei/Genome2OR.git
    ```
2. Enter the directory
    ```
    cd ./Genome2OR/scripts 
    ```

## Usage
### Main modules
* [nhmmer.py](#nhmmer): Simplify running program nhmmer program.
* [FindOR.py](#findor): Extract olfactory receptor cds from genome.
* [IdentifyFunc.py](#identity): Recognition function OR gene
* [Iteration.py](#iterate): Iterration annotated a species
* [batch.py](#batch): Batch annotated genome.

#### <span id='nhmmer'>nhmmer.py Usage</span>
```
python nhmmer.py -h[--help, None]
```
```
usage: run_nhmmer [-h] [-e] [-c] [-v] [-V] profile genome output

Autorun nhmmer

positional arguments:
    profile              String, profile nhmmer need(hmm, [un]alignment file)
    genome               String, genomic data file path.
    output               String, save nhmmer output file path.

optional arguments:
    -h, --help           show this help message and exit
    -e , --EvalueLimit   Float, Sequence similarity threshold. (default:1e-10)
    -c , --cpus          number of parallel CPU workers to use. (default='2/3 of
                         all cores')
    -v, --verbose        Print verbose information.
    -V, --version        Show version message and exit.

http://zhaolab.shanghaitech.edu.cn/

```

#### <span id="findor">FindOR.py Usage</span>
```
python FindOR.py -h[--help, None]
```
```
usage: FindOR [-h] [-o] [-p] [-e] [-l] [-v] [-V] input genome

Olfactory receptor annotation

positional arguments:
    input                 String, nhmmer output file path.
    genome                String, Genomic data file path.

optional arguments:
    -h, --help            show this help message and exit
    -o , --outputdir      String, Result save directory.(default:../output)
    -p , --prefix         String, output file prefix.(default:ORannotation)
    -e , --EvalueLimit    Float, Sequence similarity threshold. (default:1e-60)
    -l , --SeqLengthLimit 
                          Int, An artificially set OR's sequence length
                          threshold.(default:868)
    -v, --verbose         Print verbose information.
    -V, --version         Show version message and exit.

http://zhaolab.shanghaitech.edu.cn/
```

#### <span id='identity'>IdentifyFunc.py Usage</span>
```
python IdentityFunc.py -h[--help, None]
```
```
usage: IdentityFunc.py [-h] [-o] [-p] [-c] [-k] [-v] [-V]
                       hitPROfile hitDNAfile

Idntity Function OR

positional arguments:
    hitPROfile         IdentityFunc.py script output file(protein sequence
                       file). hit sequence from genome
    hitDNAfile         IdentityFunc.py script output file(DNA sequence file).
                       hit sequence from genome

optional arguments:
    -h, --help         show this help message and exit
    -o , --outputdir   String, Result save directory.(default:../output)
    -p , --prefix      String, output file prefix.(default:Identity)
    -c , --cpus        number of parallel CPU workers to use. (default='2/3 of
                       all cores')
    -k , --keepfile    Bool, whether to keep intermediate file.(default:True)
    -v, --verbose      Print verbose information.
    -V, --version      Show version message and exit.

http://zhaolab.shanghaitech.edu.cn/

```

#### <span id='iterate'>Iteration.py Usage</span>
```
python Iteration.py -h[--help, None]
```
```
usage: Iteration.py [-h] [-i] [-e] [-l] [-c] [-p] [-v] [-V]
                    profile outputdir genome

Iteration annotated a genome

positional arguments:
    profile               String, profile nhmmer need(hmm, [un]alignment file).
                          Notice: A group of genomes share a profile.
    outputdir             String, output directory.
    genome                String, genomic data file.
  
optional arguments:
    -h, --help            show this help message and exit
    -i , --iteration      Int, Number of iterations.default=2
    -e , --EvalueLimit    Float, Sequence similarity threshold. (default:1e-20)
    -l , --SeqLengthLimit 
                          Int, An artificially set OR's sequence length
                          threshold.(default:868)
    -c , --cpus           number of parallel (default='2/3 of all cores')
    -p , --prefix         String, output file prefix.(default:Identity)
    -v, --verbose         Print verbose information.
    -V, --version         Show version message and exit.

http://zhaolab.shanghaitech.edu.cn/

```

#### <span id='batch'>batch.py Usage</span>
```
python batch.py -h[--help, None]
```
```
sage: BatchProcess.py [-h] [-e] [-l] [-c] [-v] [-V]
                   profile inputdir nhmmerout outputdir

Batch annotated genome

positional arguments:
    profile               String, profile nhmmer need(hmm, [un]alignment file).
                          Notice: A group of genomes share a profile.
    inputdir              String, genomic directory.
    nhmmerout             String, run nhmmer program result.
    outputdir             String, processing results directory, run FindOR.py
                          result.
  
optional arguments:
    -h, --help            show this help message and exit
    -e , --EvalueLimit    Float, Sequence similarity threshold. (default:1e-60)
    -l , --SeqLengthLimit 
                          Int, An artificially set OR's sequence length
                          threshold.(default:868)
    -c , --cpus           number of parallel CPU workers to use. (default='2/3
                          of all cores')
    -v, --verbose         Print verbose information.
    -V, --version         Show version message and exit.
  
http://zhaolab.shanghaitech.edu.cn/
 
```

### Other modules
* Genome2OR/data
    * statistic_nterm.py
    * statistic_pattern_match.py
    * statistic_conserved_site.py
* Genome2OR/scripts/src
    * Functions: Functions module.
    * ParseArgs: Command line parsing module.
    * CodeMessages: Error message module.
    * config: Configuration module.
* Genome2OR/template
    * '*.hmm' file, HMM profiles
    * 'template.fasta' is template OR sequence file.You can replace it with your template file, but the first sequence(OR5AN1) cannot be change.Note that too many squences in 'template.fasta' will cause the program to run slowly, typically no more than five sequenses. 

## Example


*[Welcome to our laboratory website](http://zhaolab.shanghaitech.edu.cn/)*


## Contact us
* zhaosw@shanghaitech.edu.cn, Suwen Zhao
* hanwei@shanghaitech.edu.cn, Wei Han
