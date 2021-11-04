#!/usr/bin/env python3
# coding: utf-8

import os
from pathlib import Path

def get_names(name_file : str) -> list:
    names = []
    try:
        with open(name_file,'r') as f:
            for line in f :
                if len(line.strip().split(' ')) > 1:
                    raise ValueError(f'Incorrect \'{name_file}\' names file !')
                names.append(line.strip())
    except ValueError:
        raise
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
    
    ####Tests
    try:
        if not Path(args.names).is_file():
            raise FileNotFoundError(f'File \'{args.names}\' doesn\'t exist !')
    except FileNotFoundError:
        raise

    try:
        if not Path(args.repository).exists():
            raise FileNotFoundError(f'\'{args.repository}\' doesn\'t exists !')
        elif not Path(args.repository).is_dir():
            raise NotADirectoryError(f'\'{args.repository}\' isn\'t a directory !')
    except NotADirectoryError:
        raise
    except FileNotFoundError:
        raise
        
    try:
        if not Path(args.out).resolve().parent.exists() or Path(args.out).resolve().parent.is_file():
            raise FileNotFoundError(f'Directory of the output file \'{args.out}\' doesn\'t exist ! Please create the directory first!')
    except FileNotFoundError:
        raise
    output = Path(args.out)
    output.touch(exist_ok=True) 
    try:
        if not output.is_file():
            raise IsADirectoryError(f'\'{args.out}\' is a directory not an output file !')
    except IsADirectoryError:
        raise
    ####
    
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
