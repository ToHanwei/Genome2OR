# Genome2OR
Annotate Olfactory receptor CDS from genome
==============================================


## Contents Index
* [Abstract](#abstract)
* [Quick start](#quick-start)
* [Download script](#download-script)
* [Usage](#usage)
* [Example](#example)


## Abstract
[Genome2OR](https://github.com/ToHanwei/Genome2OR.git) is a genetic annotation tool based on 
[HMMER](http://hmmer.org), [MAFFT](https://mafft.cbrc.jp/alignment/software/) and 
[CD-HIT](http://weizhongli-lab.org/cd-hit/).  
[HMMER](http://hmmer.org) searches biological sequence databases for
homologous sequences, using either single sequences or multiple
sequence alignments as queries. HMMER implements a technology called
"profile hidden Markov models" (profile HMMs).  

[MAFFT](https://mafft.cbrc.jp/alignment/software/) is a multiple sequence alignment program 
for unix-like operating systems. 
It offers a range of multiple alignment methods, L-INS-i (accurate; for alignment of <∼200 sequences), 
FFT-NS-2 (fast; for alignment of <∼30,000 sequences), etc.

[CD-HIT](http://weizhongli-lab.org/cd-hit/) is very fast and can handle extremely large databases. 
CD-HIT helps to significantly reduce the computational and manual efforts in many sequence analysis 
tasks and aids in understanding the data structure and correct the bias within a dataset.

## Quick start
1. Execute nhmmer.py
    ```
   cd YouDir/Genome2OR/scripts
   python nhmmer.py ../template/Mammalia.hmm genome.fasta nhmmer_out.tblout -v
    ```
   Mammalia.hmm is the HMM profile for mammalia. You can find more species of HMM profile in YouDir/Genome2OR/template directory if needed.
   genome.fasta is genome(DNA) that needs annotation.

2. Execute FindOR.py
    ```
   cd YourDir/Genome2OR/scripts
   python FindOR.py nhmmer_out.tblout genome.fasta -o ../output -o ORannotation -v 
    ```
   
3. Execute IdentityFunc.py
    ```
   cd YourDir/Genome2OR/scripts
   python IdentityFunc.py ../output/ORannotation_ORs_pro.fa -o ../output -p Identity -v
    ```
   After running 1,2,3, you will find file Identity_func_ORs.fasta in directory 
   ../output, which is the Olfactory receptor we finally found
   
4. Execute batch.py
    ```
   cd YourDir/Genome2OR/scripts
   python batch.py profile.hmm inputdir nhmmeroutdir outputdir -v
    ```
   For batch annotation of genomes.

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
usage: IdentityFunc.py [-h] [-o] [-p] [-v] [-V] hitfile

Idntity Function OR

positional arguments:
    hitfile            FindOR.py script output file(protein sequence file). hit
                       sequence from genome

optional arguments:
    -h, --help         show this help message and exit
    -o , --outputdir   String, Result save directory.(default:../output)
    -p , --prefix      String, output file prefix.(default:Identity)
    -v, --verbose      Print verbose information.
    -V, --version      Show version message and exit.

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
