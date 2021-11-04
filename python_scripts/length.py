#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
from itertools import groupby
from pathlib import Path

def iter_fasta (file: str) -> str :
    '''Iteration on the genomes'''
    f = open(file, 'rb')
    faiter = (x[1] for x in groupby(f, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield (seq)
    f.close()

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def get_length(genomes,out,seed_size,false_positif_rate):
    ####Tests
    try:
        if not Path(out).resolve().parent.exists() or Path(out).resolve().parent.is_file():
            raise FileNotFoundError(f'Directory of the output file \'{out}\' doesn\'t exist ! Please create the directory first!')
    except FileNotFoundError:
        raise
    output = Path(out)
    output.touch(exist_ok=True) 
    try:
        if not output.is_file():
            raise IsADirectoryError(f'\'{out}\' is a directory not an output file !')
    except IsADirectoryError:
        raise
    try:
        if not Path(genomes).exists():
            raise FileNotFoundError(f'\'{genomes}\' doesn\'t exists !')
        elif not Path(genomes).is_dir():
            raise NotADirectoryError(f'\'{genomes}\' isn\'t a directory !')
    except NotADirectoryError:
        raise
    except FileNotFoundError:
        raise
    if seed_size != None:
        try :
            if not seed_size.isnumeric():
                raise ValueError(f'Seed_size \'{seed_size}\' is not a valid number ! (must be a positive integer)')
        except ValueError:
            raise
    if false_positif_rate != None:
        try:
            if not isfloat(false_positif_rate):
                raise ValueError(f'False positif rate \'{false_positif_rate}\' is not a float !')
            elif float(false_positif_rate) < 0 or float(false_positif_rate) > 1:
                raise ValueError(f'False positif rate \'{false_positif_rate}\' is not valid ! (must be between 0 and 1)')
            false_positif_rate = float(false_positif_rate)
        except ValueError:
            raise
    ####
    
    maximum_length = 0
    genomes_length = {}
    genomes_order = []
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
                genomes_order.append('.'.join(genome.split('.')[:-1])) 
                genomes_length['.'.join(genome.split('.')[:-1])] = (length,nb)
        for g in sorted(genomes_order):
            o.write(f"{g}\t{genomes_length[g][0]}\t{genomes_length[g][1]}\n")
    if seed_size != None:
        with open('bf_min_size.txt' ,'w') as out:
            print(f'Maximum genome lenght = {maximum_length} and false positif rate = {false_positif_rate}')
            len_bf = int(round(-maximum_length/(np.log(1-false_positif_rate))))
            out.write(f'{len_bf}')
            print(f'{len_bf}')

def main(args):
    get_length(args.genomes,args.out,args.seed_size,args.false_positif_rate)
    

if __name__ == "__main__":
    main()
    