#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def improve_results(file,number_name_list):
    
    #get the number for each strains
    numbername = {}
    with open(number_name_list,'r') as f:
        for line in f:
            if line != '':
                line = line.strip()
                line = line.split('\t')
                
                number = line[0]
                strain = line[1]
                numbername[int(number)] = strain
    
    with open(file,'r') as f:
        strains_list = []
        percents_list = []
        label = []
        
        for line in f :
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
                            all_strains = all_strains+'\n'+numbername[int(num)]
                strains_list.append(all_strains)
            else:
                percents_list.append(int(number*100))
                strains_list.append(strain)
            label.append(strain)
        
        return percents_list, strains_list, label

def func(pct, allvalues): 
    absolute = int(pct / 100.*np.sum(allvalues)) 
    return "{:.1f}%".format(pct, absolute) 

def create_graph(percents_list, strains_list, label):
 
    colors=['bisque','orange','lightblue','lightgreen','lime','brown','mediumturquoise','ivory']
        
    fig, ax = plt.subplots(figsize =(13, 13)) 
    
    wedges, texts, autotexts = ax.pie(percents_list,     
                              labels = label,
                              autopct = lambda pct: func(pct, percents_list),
                              colors=colors,
                              textprops = dict(color ="black")) 

    ax.legend(wedges, strains_list,  
          loc ="center left", 
          bbox_to_anchor =(0.9, 0.45, 0, 1),
          prop={'size': 10}) 
    
    plt.setp(autotexts, size = 12, weight ="bold") 

    plt.savefig('identification.png', format='png')
    
def create_output_file(percents_list, strains_list,output):
    with open(output,'w') as o:
        o.write(f'***** ORI *****\n')
        print(f'***** ORI *****\n')
        for i,n in enumerate(percents_list):
            s=''
            if '\n' in strains_list[i]:                   
                s = ' , '.join(strains_list[i].split('\n'))
            print(f'{n}%\t{s}\n')
            o.write(f'{n}%\t{s}\n')
        print(f'***************\n')
        o.write(f'***************\n')
    
def main(args):
    
    percents_list, strains_list, label = improve_results(args.file,args.number_name_list)
    
    if args.pie_chart:
        create_graph(percents_list, strains_list, label)
    
    create_output_file(percents_list, strains_list,args.output)
    
    
    
    
    
    
    
    
    
    