#!/usr/bin/env python3
# coding: utf-8

import numpy as np

def merge_length(bf,length,correspondance,out):
    with open(out,'w') as o:

        #Get the number id corresponding to each strain
        corresp = {}
        with open(correspondance,'r') as c:
            for line in c:
                line = line.strip().split('\t')
                if '.fna' in line[1]:
                    strain = line[1].split('.fna')[0]
                else:
                    strain = line[1].split('.fasta')[0]
                num = int(line[0])
                corresp[num]=strain

        #Get the length of each strain
        leng = {}                
        with open(length,'r') as l:
            for line in l:
                line = line.strip().split('\t')
                strain = line[0]
                longueur = line[1]
                leng[strain]= int(longueur)
        
        #For each cluster of strains get the mean length of the cluster
        with open(bf,'r') as b:
            for line in b:
                longueurs = []
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
                            longueurs.append(leng[corresp[int(num)]])
                #Write the length in the output file
                o.write(f'{strain}\t{int(np.mean(longueurs))}\n')
        

def main(args):
    
    merge_length(args.bf,args.length,args.correspondance,args.out)
    

if __name__ == "__main__":
    main()
    