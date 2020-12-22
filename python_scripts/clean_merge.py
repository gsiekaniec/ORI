#!/usr/bin/env python3
# coding: utf-8

import os

def get_names(name_file : str) -> list:
    names = []
    with open(name_file,'r') as f:
        for line in f :
            names.append(line.strip())
    return names
    
def bad_files_suppression_list (repository : str, names : list) -> list:
    to_suppr = []
    for element in os.listdir(repository):
        if element.endswith('.bf'):
            element = element.strip()
            if element not in names:
                if (''.join(element.split('_')[:-1])).isnumeric():
                    for number in element.split('_')[:-1]:
                        number = int(number)
                        to_suppr.append(names[number])
    return to_suppr
    
def main(args):
    
    names = get_names(args.names)
    
    #creation of the correspondances {number/file} file 
    with open(args.repository+'/'+args.out,'w') as o:    
        for num,name in enumerate(names):
            o.write(f'{num}\t{".".join(name.split(".")[:-2])}\n')
    
    to_suppr = bad_files_suppression_list(args.repository, names)
    
    # suppression of the bf file already merged
    for bf_file in to_suppr:
        if bf_file in names:
            if os.path.exists(args.repository+'/'+bf_file):
                print(f'Suppresion of {args.repository}/{bf_file}')
                os.remove(args.repository+'/'+bf_file) 
            

if __name__ == "__main__":
    main()
