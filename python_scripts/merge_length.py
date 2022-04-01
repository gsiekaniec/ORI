#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from pathlib import Path

def merge_length(bf,length,correspondance,out):
    
    ####Tests
    try:
        if not Path(bf).is_file():
            raise FileNotFoundError(f'File \'{bf}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
        
    try:
        if not Path(length).is_file():
            raise FileNotFoundError(f'File \'{length}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
    
    try:
        if not Path(correspondance).is_file():
            raise FileNotFoundError(f'File \'{correspondance}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
    
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
    ####
    
    with open(out,'w') as o:

        #Get the number id corresponding to each strain
        corresp = {}
        try:
            with open(correspondance,'r') as c:
                for line in c:
                    if len(line.strip().split('\t')) != 2:
                        raise ValueError(f'Incorrect \'{correspondance}\' correspondances file !')
                    line = line.strip().split('\t')
                    if '.fna' in line[1]:
                        strain = line[1].split('.fna')[0]
                    else:
                        strain = line[1].split('.fasta')[0]
                    num = int(line[0])
                    corresp[num]=strain
        except ValueError:
            raise

        #Get the length of each strain
        leng = {}   
        try:             
            with open(length,'r') as l:
                for line in l:
                    if len(line.strip().split('\t')) != 3:
                        raise ValueError(f'Incorrect \'{length}\' lengths file !')
                    line = line.strip().split('\t')
                    strain = line[0]
                    longueur = line[1]
                    leng[strain]= int(longueur)
        except ValueError:
            raise
        
        #For each cluster of strains get the mean length of the cluster
        try:
            with open(bf,'r') as b:
                for line in b:
                    longueurs = []
                    if len(" ".join(line.strip().split('\t')).split(' ')) > 1:
                         raise ValueError(f'Incorrect \'{bf}\' names file !')
                    strain = line.strip()
                    
                    strain = strain.split('.bf')[0]
                    if '.fna' in strain:
                        strain = strain.split('.fna')[0]
                    elif '.fasta' in strain:
                        strain = strain.split('.fasta')[0]
                        
                    if strain in leng.keys():
                        longueurs = [leng[strain]]
                    else:
                        for num in strain.split('_'):
                            if num != '':
                                try:
                                    longueurs.append(leng[corresp[int(num)]])
                                except ValueError:
                                    raise ValueError(f'Names in \'{bf}\' names file doesn\'t correspond to names in \'{correspondance}\'!')
                    #Write the length in the output file
                    o.write(f'{strain}\t{int(np.mean(longueurs))}\n')
        except ValueError:
            raise
        

def main(args):
    
    merge_length(args.bf,args.length,args.correspondance,args.out)
    

if __name__ == "__main__":
    main()
    
