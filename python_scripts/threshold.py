#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt


def threshold_determination_help (matrix,threshold):
    
    threshold = float(threshold)
    distances = []
    with open (matrix,'r') as f:
        for line in f:
            line = line.strip().split(' ')
            for number in line:
                number=float(number)*10000.0
                distances.append(round(number))
    
    #maximum = float(max(distances)+0.0005)*10000
    maximum = int(max(distances)+1)
    #pas = int(maximum/0.00001)
    pas = maximum*2*2
    
    ax = plt.subplot()
    ax.hist(distances, bins = pas, color='blue', density=True) 

    #ax.legend((['1e-05']), loc='upper right');
    plt.xlabel('Hamming distances (1e-05)')
    plt.ylabel('Density')
    plt.title('Hamming distance between pairs of strains.')
    
    #plt.xlim(-0.0001,maximum)
    plt.xlim(-1,maximum)
    
    plt.axvline(x=threshold*10000, color='r', linestyle='--',lw=0.5)
    
    plt.tight_layout()
    
    plt.savefig('threshold.png', format='png')    

def main(args):
    
    threshold_determination_help(args.matrix, args.threshold)
    
      
    
    
    
    
    
