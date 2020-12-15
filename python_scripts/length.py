#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
from itertools import groupby

def iter_fasta (file: str) -> str :
    '''Iteration on the genomes'''
    f = open(file, 'rb')
    faiter = (x[1] for x in groupby(f, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield (seq)
    f.close()

def get_length(genomes,out,seed_size,false_positif_rate):
    maximum_length = 0
    with open(out,'w') as o:
        for genome in os.listdir(genomes):
            length = 0
            nb = 0
            if genome.endswith(".fasta") or genome.endswith(".fa") or genome.endswith(".fna"):        
                for it in iter_fasta (genomes+'/'+genome):
                    nb+=1
                    for seq in it:
                        length = length+len(seq)
                if seed_size != None:
                    if maximum_length < (int(length)-(int(seed_size)-1)):
                        maximum_length = int(length)-(int(seed_size)-1)
                o.write(f"{'.'.join(genome.split('.')[:-1])}\t{length}\t{nb}\n")
    if seed_size != None:
        with open('bf_min_size.txt' ,'w') as out:
            print(f'Maximum genome lenght = {maximum_length} and false positif rate = {false_positif_rate}')
            len_bf = int(round(-maximum_length/(np.log(1-false_positif_rate))))
            out.write(f'{len_bf}')
            print(f'{len_bf}')

def main(args):
    get_length(args.genomes,args.out,args.seed_size,float(args.false_positif_rate))
    

if __name__ == "__main__":
    main()
    