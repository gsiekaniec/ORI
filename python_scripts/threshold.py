#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.cluster.hierarchy import ClusterWarning
from warnings import simplefilter
simplefilter("ignore", ClusterWarning)
simplefilter("ignore", UserWarning)
import numpy as np
import pandas as pd


def order_names(size : int, names : [],new_order : []) -> []:
    '''Gets the names of the strains in the new order after the hierarchical clustering'''
    new_names = ['']*size
    for i in range(size):
        new_names[i]=names[new_order[i]]
    return new_names
            
def new_list (data : [], names : [], threshold : float):
    '''Returns a dictionary containing in key the strains that have a sibling strain and in value their number of sibling strains'''
    data = np.array(data)
    name_to_save = {}
    for i in range(len(data)):
        for j in range(len(data)):
            if i != j:
                if data[i,j] <= threshold:
                    if names[i] in name_to_save.keys():
                        name_to_save[names[i]] +=1
                    else:
                        name_to_save[names[i]] =1
    return name_to_save

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def threshold_determination_help (mat,threshold, name, output):
    
    ####Tests
    try:
        if not Path(output).resolve().parent.exists() or Path(output).resolve().parent.is_file():
            raise FileNotFoundError(f'Directory of the output file \'{output}\' doesn\'t exist ! Please create the directory first!')
    except FileNotFoundError:
        raise
    out = Path(output)
    out.touch(exist_ok=True) 
    try:
        if not out.is_file():
            raise IsADirectoryError(f'\'{output}\' is a directory not an output file !')
    except IsADirectoryError:
        raise
    try:
        if not isfloat(threshold):
            raise ValueError(f'Threshold \'{threshold}\' is not a number !')
        elif float(threshold) < 0 or float(threshold) > 1:
            raise ValueError(f'Threshold \'{threshold}\' is not a valid number ! (must be between 0 and 1)')
    except ValueError:
        raise
    try:
        if not Path(mat).is_file():
            raise FileNotFoundError(f'File \'{mat}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
    try:
        if not Path(name).is_file():
            raise FileNotFoundError(f'File \'{name}\' doesn\'`t exist !')
    except FileNotFoundError:
        raise            
    ####
    
    threshold = float(threshold)
    print(f'Threshold used: {threshold}\n')
    
    try:
        matrix = np.loadtxt(mat, comments='#', skiprows=0, usecols=None, unpack=False, ndmin=0) #Load matrix
        matrix_length = len(matrix)
    except ValueError:
        raise ValueError(f'File \'{mat}\' containing something other than Hamming distances')
    try:
        names = [i.strip().split('.')[0] for i in open(name,'r')] #Load names
    except UnicodeDecodeError:
        print(f'Wrong names file \'{name}\' !')
        raise
    try:
        if not matrix_length == len(names):
            raise ValueError(f'File \'{name}\' doesn\'t contain the same number of strains as the Hamming distance matrix')
    except ValueError:
        raise
    #print(2.5/np.sqrt(len(names)))
    
    #Get distances to print histogram
    df = pd.DataFrame(matrix, columns = names, index = names) #Create a dataframe from the matrix
    distances = [round(float(number)*10000.0) for sublist in df.values.tolist() for number in sublist]

    

    #######Hamming distances histogram########
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
    plt.close()    
    ######################
    
    sns.set(font_scale=(2.5/np.sqrt(len(names)))) #Set the labels font size
    
    #######Clustermap########
    clust = sns.clustermap(df, cmap = 'viridis', method = 'ward',xticklabels = 1, yticklabels = 1, cbar_kws={'label': 'Hamming distance'})
    plt.savefig('hc_heatmap.svg', format="svg")
    plt.close()
    ######################
    
    #######Heatmap with threshold filter########
    data = clust.data2d #get the clustering matrix
    new_order = clust.dendrogram_col.reordered_ind
    new_names = order_names(len(matrix), names, new_order)
    df = pd.DataFrame(data, columns = new_names, index = new_names) #new dataframe with the clustering matrix
    sns.heatmap(df, cbar = True, cmap = 'viridis', vmax = threshold+(threshold/10), xticklabels = 1, yticklabels = 1, cbar_kws={'label': 'Hamming distance'}) 
    plt.tight_layout()
    plt.savefig('heatmap_threshold.svg', format="svg")
    ######################
    
    #######Write list the names of strains that have at least one sibling strain########
    d = new_list (data, new_names, threshold)
    with open(out,'w') as o:
        for name in new_names:
            if name in d.keys():
                o.write(f'{name}\t{d[name]}\n')
            else:
                o.write(f'#\n')
    
    ######################

def main(args):
    
    threshold_determination_help(args.matrix, args.threshold, args.names, args.output)
    
      
    
    
    
    
    
