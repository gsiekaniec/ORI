#!/usr/bin/env python3
# coding: utf-8

import os
from itertools import groupby

def iter_fasta (file: str) -> str :
    '''Iteration on the genomes'''
    f = open(file, 'rb')
    faiter = (x[1] for x in groupby(f, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield (seq)
    f.close()

def get_length(genomes,out):
    with open(out,'w') as o:
        for genome in os.listdir(genomes):
            length = 0
            nb = 0
            if genome.endswith(".fasta") or genome.endswith(".fa") or genome.endswith(".fna"):        
                for it in iter_fasta (genomes+'/'+genome):
                    nb+=1
                    for seq in it:
                        length = length+len(seq)
                o.write(f"{'.'.join(genome.split('.')[:-1])}\t{length}\t{nb}\n")

def main(args):
    get_length(args.genomes,args.out)
    

if __name__ == "__main__":
    main()
    