#!/usr/bin/env python3
# coding: utf-8

import os
import statistics
import math
import gzip
from pathlib import Path

def iterFastq (file: str, threshold: int, length:int, gz: bool, nb_read_tot: int, general_length: int) -> (tuple):
    '''Read a fastq file'''
    if gz:
        f = gzip.open(file, 'rb')
    else :
        f = open(file,'rb')
        
    for line in f:
        if str(line,'utf-8')[0] == '@':
            header = str(line,'utf-8')
            nb_read_tot +=1
            seq = str(f.readline(),'utf-8')
            f.readline()
            qual = []
            seqQual = str(f.readline(),'utf-8')
            general_length += len(seq)
            if len(seq) > length:
                for char in seqQual:
                    q = ord(char)-33
                    error_rate = 10**(-(q/10))
                    qual.append(error_rate)
                read_quality = -10*(math.log10(statistics.mean(qual)))
                if read_quality >= threshold :
                    yield (header, seq, seqQual)
            yield (nb_read_tot,general_length)
            
    f.close()  

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def run (file: str, threshold: float, length:int, gz: bool):

    
    ####Tests
    try:
        if not Path(file).is_file():
            raise FileNotFoundError(f'File \'{file}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
        
    
    if threshold != None:
        try:
            if not isfloat(threshold):
                raise ValueError(f'Quality threshold \'{threshold}\' is not an integer !')
            elif float(threshold) < 0:
                raise ValueError(f'Quality threshold \'{threshold}\' is not valid ! (must be positive)')
            threshold = float(threshold)
        except ValueError:
            raise
    
    if length != None:
        try :
            if not isinstance(length, (int)):
                if not length.isnumeric():
                    raise ValueError(f'Seed_size \'{length}\' is not a valid number ! (must be a positive integer)')
            length = int(length)
        except ValueError:
            raise
    
    ####
    
    out = file_to_save (file, threshold)
    print (f"Write reads in : {out}")
    o = open(out,'w')
    nb_read_tot = 0
    nb_read_end = 0
    general_length = 0
    length_sup = 0
    for i in iterFastq(file,threshold,length,gz,nb_read_tot,general_length):
        if len(i) == 2:
            nb_read_tot = i[0]
            general_length = i[1]
        elif len(i) == 3 :
            nb_read_end+=1
            length_sup += len(i[1])
            o.write(f'{i[0]}\n{i[1]}\n+\n{i[2]}\n')
        
    print (f"Reads number : {nb_read_tot}\nReads number after filter : {nb_read_end}")
    print (f"Mean reads length : {int(general_length/nb_read_tot)}\nMean reads length after filter : {int(length_sup/nb_read_end)}")
    o.close()
    
def file_to_save (file: str, threshold: int) -> str:
    out = file.split('.fastq')[0]+"_better_than_"+str(int(threshold))+".fastq"
    return out
    
def main(args):
    run(args.fastq,args.qualityMin,args.lengthMin,args.gzip)
