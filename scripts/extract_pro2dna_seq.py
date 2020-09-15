#!:coding:utf-8

import os
import sys


def proc_seq(infile):
    seqdict = {}
    with open(infile) as seqf:
        seqs = seqf.read().split('>')[1:]
    for seq in seqs:
        lines = seq.split('\n')
        name = lines[0]
        aads = ''.join(lines[1:])
        seqdict[name] = aads
    return(seqdict)


def main():
    prodir = sys.argv[1]
    dnadir = sys.argv[2]
    outdir = sys.argv[3]
    profiles = os.listdir(prodir)
    for profile in profiles:
        outfile = profile.replace('_pro.', '_dna.')
        dnafile = profile.split('_func_')[0] + '_func_ORs_dna.fasta'
        profile = os.path.join(prodir, profile)
        dnafile = os.path.join(dnadir, dnafile)
        outfile = os.path.join(outdir, outfile)
        pros = proc_seq(profile)
        dnas = proc_seq(dnafile)
        with open(outfile, 'w') as outf:
            for name in pros.keys():
                line = '>' + name + '\n' + dnas[name] + '\n'
                outf.write(line)
        
    

if __name__ == "__main__":
    main()
