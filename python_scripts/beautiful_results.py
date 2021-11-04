#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def improve_results(file,number_name_list):
    
    #get the number for each strains
    numbername = {}
    try:
        if number_name_list != None:
            with open(number_name_list,'r') as f:
                for line in f:
                    if line != '':
                        if len(line.strip().split('\t')) != 2:
                            raise ValueError(f'Wrong \'{number_name_list}\' correspondance file !')
                        line = line.strip()
                        line = line.split('\t')
                        
                        number = line[0]
                        strain = line[1]
                        numbername[int(number)] = strain
    except ValueError:
        raise
    
    with open(file,'r') as f:
        try:
            strains_list = []
            percents_list = []
            
            for line in f :
                if len(line.strip().split(' : ')) != 2:
                                raise ValueError(f'Wrong \'{file}\' result file !')
                line = line.strip()
                line = line.split(' : ')
    
                strain = line[0]
                number = float(line[1])
                
                if ''.join(strain.split('_')).isnumeric():
                    all_strains = ''
                    percents_list.append(int(number*100))
                    strains = strain.split('_')
                    for num in strains:
                        if num != '':
                            if all_strains == '':
                                all_strains=numbername[int(num)]
                            else:
                                all_strains = all_strains+' , '+numbername[int(num)]
                    strains_list.append(all_strains)
                else:
                    percents_list.append(int(number*100))
                    strains_list.append(strain)
            
            list_complete = []
            for i,val in enumerate(percents_list):
                list_complete.append((strains_list[i],val))    
            list_complete.sort(key=lambda x: x[1], reverse=True)
        except ValueError:
            raise
        return list_complete

def func(pct, allvalues): 
    absolute = int(pct / 100.*np.sum(allvalues)) 
    return "{:.1f}%".format(pct, absolute) 

def create_graph(complete_list):
 
    colors=['bisque','orange','lightblue','lightgreen','lime','brown','mediumturquoise','ivory']
        
    fig, ax = plt.subplots(figsize =(13, 13)) 
    
    wedges, texts, autotexts = ax.pie([i[1] for i in complete_list],     
                              labels = [i[0] for i in complete_list],
                              autopct = lambda pct: func(pct, [i[1] for i in complete_list]),
                              colors=colors,
                              textprops = dict(color ="black")) 

    ax.legend(wedges, [i[0] for i in complete_list],  
          loc ="center left", 
          bbox_to_anchor =(0.9, 0.45, 0, 1),
          prop={'size': 10}) 
    
    plt.setp(autotexts, size = 12, weight ="bold") 

    plt.savefig('identification.png', format='png')
    
def create_output_file(complete_list,output):
    with open(output,'w') as o:
        o.write(f'***** ORI *****\n')
        print(f'***** ORI *****\n')
        for i, j in complete_list:
            s = i
            n = j
            print(f'{n}%\t{s}\n')
            o.write(f'{n}%\t{s}\n')
        print(f'***************\n')
        o.write(f'***************\n')
    
def main(args):
    
    ####Tests
    try:
        if not Path(args.file).is_file():
            raise FileNotFoundError(f'File \'{args.file}\' doesn\'t exist !')
    except FileNotFoundError:
        raise

    try:
        if not Path(args.number_name_list).is_file():
            raise FileNotFoundError(f'File \'{args.number_name_list}\' doesn\'t exist !')
    except FileNotFoundError:
        raise    
    
    try:
        if not Path(args.output).resolve().parent.exists() or Path(args.output).resolve().parent.is_file():
            raise FileNotFoundError(f'Directory of the output file \'{args.output}\' doesn\'t exist ! Please create the directory first!')
    except FileNotFoundError:
        raise
    output = Path(args.output)
    output.touch(exist_ok=True) 
    try:
        if not output.is_file():
            raise IsADirectoryError(f'\'{args.output}\' is a directory not an output file !')
    except IsADirectoryError:
        raise
    ####
    
    complete_list = improve_results(args.file,args.number_name_list)
    
    if args.pie_chart:
        create_graph(complete_list)
    
    create_output_file(complete_list,args.output)
    
    
    
    
    
    
    
    
    
    
