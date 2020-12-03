#!/usr/bin/env python3
# coding: utf-8

import numpy as np 

def iter_file(file: str, listbf: list, grid) :
    '''Read the output file of HowDeSBT and create the {Strains x Reads} matrix '''
    col=-1
    with open(file, 'r') as f:
        for line in f:
            if line[0] == '*':
                col +=1
            else:
                l = line.strip()
                if '.fna' in l:
                    if len(l.split('.fna')) > 1:
                        strain = str(l.split('.bf')[0].split('.fna')[0])
                    else:
                        strain = str(l.split(' ')[0].split('.bf')[0].split('.fna')[0])
                elif '.fasta' in l:
                    if len(l.split('.fasta')) > 1:
                        strain = str(l.split('.bf')[0].split('.fasta')[0])
                    else:
                        strain = str(l.split(' ')[0].split('.bf')[0].split('.fasta')[0])
                elif '.bf' in l:
                    if len(l.split('.bf')) > 1:
                        strain = str(l.split('.bf')[0])
                    else:
                        strain = str(l.split(' ')[0].split('.bf')[0])
                else:
                    strain = l.split(' ')[0]
                number = float(l.split(' ')[2])
                
                li = listbf[strain]
                grid[col,li] = number
    return grid

def empty_grid(N: int, M: int) -> np.zeros:
    '''Create an empty matrix'''
    grid = np.zeros(shape=(N,M))
    return grid

def main(args):

    #Count the number of reads
    with open(args.file, 'r') as f:
        nb_reads = 0
        for line in f: 
            if line[0] == '*': nb_reads+=1
            
    #Get the file name
    with open(args.list_name,'r') as f:
        listbfdict = {}
        i = 0
        for line in f:
            line = line.strip()
            if '.bf' in line:
                line = line.split('.bf')[0]
            if '.fasta' in line:
                line = line.split('.fasta')[0]
            elif '.fna' in line:
                line = line.split('.fna')[0]
            listbfdict[line]=i
            i+=1
            
    #Creation of an empty matrix
    grid = empty_grid(int(nb_reads),int(i))
    
    #Read the output file of HowDeSBT
    grid = iter_file (args.file, listbfdict, grid)
    
    #Write the {Strains x Reads} matrix
    with open(args.out,'w') as o:
        for i in grid:
            for j in i:
                o.write(str(j)+' ')
            o.write('\n')







