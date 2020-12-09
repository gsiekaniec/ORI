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
                distances.append(float(number))
    
    maximum = float(max(distances)+0.0005)
    pas = int(maximum/0.00005)
    
    ax = plt.subplot()
    ax.hist(distances, bins = pas, density=True) 
    
    plt.xlabel('Hamming distances')
    plt.ylabel('Density')
    plt.title('Hamming distance between pairs of strains.')
    
    plt.xlim(-0.0001,maximum)
    
    plt.axvline(x=threshold, color='r', linestyle='--',lw=0.5)
    
    plt.tight_layout()
    
    plt.savefig('threshold.png', format='png')    

def main(args):
    
    threshold_determination_help(args.matrix, args.threshold)
    
      
    
    
    
    
    
