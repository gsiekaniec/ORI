#!/usr/bin/env python3
# coding: utf-8

import argparse 
import python_scripts.matrix
import python_scripts.identification
import python_scripts.length
import python_scripts.merge_length
import python_scripts.clean_merge
import python_scripts.suppr_bad_quality_reads

__version__ = '0.0.1'

if __name__ == '__main__':
    
    #Global parser####
    
    parser = argparse.ArgumentParser()
    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global arguments'
    parser.add_argument('--version', action='version',
        version=f'ORI v{__version__}',
        help='Display ORI version'
    )
    
    #Subparser
    
    subparsers = parser.add_subparsers(help='Functions')
    
    #Matrix parser
    
    parser_matrix = subparsers.add_parser('matrix',help='Create the {genes X strains} matrix using the results from HowDeSBT')
    
    required_matrix = parser_matrix.add_argument_group('required arguments')
    optional_matrix = parser_matrix.add_argument_group('optional arguments')

    required_matrix.add_argument('--file', '-f', metavar='IN_FILE', type=str,
    required=True, 
    help='Results file from HowDe output.'
    )
    required_matrix.add_argument('--list_name', '-l', metavar='IN_FILE2', type=str, 
    required=True, 
    help='Name of strains in the same order than in matrix.'
    )
    optional_matrix.add_argument('--output', '-out', dest='out', metavar='OUT_FILE' , 
    default='matrix.tsv',
    help='Output matrix file.')

    parser_matrix.set_defaults(parser_matrix=True, parser_identification=False,parser_length=False,parser_merge_length=False, parser_clean=False, parser_suppr_reads=True)

    #Identification parser

    parser_identification = subparsers.add_parser('identification', help='Identification of strains from a reads file')
    
    required_identification = parser_identification.add_argument_group('Required argument')
    optional_identification = parser_identification.add_argument_group('Optional argument')
    
    required_identification.add_argument('--matrix', '-m', dest='matrix', metavar='IN_FILE', 
    required=True, 
    help='Matrix file strains X reads.'
    )
    required_identification.add_argument('--clingo_path', '-clingo', dest='clingo_path', metavar='CLINGO_PASS', 
    required=True, 
    help='Clingo path.'
    )
    required_identification.add_argument('--file', '-f', dest='file',  metavar='IN_FILE2', 
    required=True,  
    help='Results file from Howde output.'
    )
    required_identification.add_argument('--length', '-le', dest='length', metavar='IN_FILE3',  
    required=True,  
    help='File with one genome and is length per line.'
    )
    required_identification.add_argument('--listname', '-l', dest='listname', metavar='IN_FILE4', 
    required=True, 
    help='Name of strains in the same order than in matrix.'
    )
    optional_identification.add_argument('--output', '-o', dest='output', metavar='OUTPUT_FILE', 
    default='out.txt', 
    help='Output file.'
    )
    optional_identification.add_argument('--threshold', '-t', dest='threshold',  metavar='THRESHOLD', 
    default=50, 
    help='Minimum percent value in the matrix for association between reads and species.'
    ) 
    optional_identification.add_argument('--nbchoices', '-n', dest='nbchoices', metavar='NB_CHOICE',
    default=12, 
    help='Only the nbchoices maximum values of a row are considered. Warning, must be less or equal to the number of species.'
    )    

    parser_identification.set_defaults(parser_identification=True,parser_matrix=False,parser_length=False,parser_merge_length=False, parser_clean=False, parser_suppr_reads=True)

    #Length parser

    parser_length = subparsers.add_parser('length', help='Get the length of genomes from fasta files')
    
    required_length = parser_length.add_argument_group('Required argument')
    optional_length = parser_length.add_argument_group('Optional argument')
    
    required_length.add_argument('--genomes', '-g', dest='genomes', 
    required=True, 
    help='Genomes repository containing the fasta files.'
    )
    
    optional_length.add_argument('--outfile', '-out', dest='out', 
    default='length.txt', 
    help='Out file containing one genome name and length per line.'
    )


    parser_length.set_defaults(parser_identification=False,parser_matrix=False,parser_length=True,parser_merge_length=False, parser_clean=False, parser_suppr_reads=True)

    #Merge length parser

    parser_merge_length = subparsers.add_parser('merge_length', help='Merge length for cluster of genomes')
    
    required_merge_length = parser_merge_length.add_argument_group('Required argument')
    optional_merge_length = parser_merge_length.add_argument_group('Optional argument')
    
    required_merge_length.add_argument('--bflist', '-b', dest='bf', 
    required=True,                                    
    help='File list of bf.') 
    
    required_merge_length.add_argument('--lengthfile', '-l', dest='length',
    required=True, 
    help='File containing one genome name and length per line.'
    )
    
    required_merge_length.add_argument('--correspondance', '-c', dest='correspondance',
    required=True, 
    help='File containing correspondance between numbers and genomes.'
    )
    
    optional_merge_length.add_argument('--outfile', '-out', dest='out', 
    default='length_merge.txt', 
    help='Out file containing one genome name and length per line.'
    )
    
    parser_merge_length.set_defaults(parser_identification=False,parser_matrix=False,parser_length=False,parser_merge_length=True, parser_clean=False, parser_suppr_reads=True)

    #Clean bf and create the number x bf list
    
    parser_clean = subparsers.add_parser('clean_merge',help='Suppression of merge bloom filters and creation of the number x strain list')
    
    required_clean = parser_clean.add_argument_group('required arguments')
    optional_clean = parser_clean.add_argument_group('optional arguments')

    required_clean.add_argument('--names', '-n', dest='names', 
    required=True,                        
    help='bf files names.'
    )
    
    required_clean.add_argument('--repository', '-r', dest='repository', 
    required=True,
    help='bf files repository.'
    )
    
    required_clean.add_argument('--outfile', '-out', dest='out', 
    default='list_number_file.txt', 
    help='Out file containing one id number and one genome name per line.'
    )

    parser_clean.set_defaults(parser_matrix=False, parser_identification=False,parser_length=False,parser_merge_length=False, parser_clean=True, parser_suppr_reads=True)

    #Supression of poor quality reads
    
    parser_suppr_reads = subparsers.add_parser('suppr_bad_reads',help='Suppression of poor quality reads from fastq')
    
    required_suppr_reads = parser_suppr_reads.add_argument_group('required arguments')
    optional_suppr_reads = parser_suppr_reads.add_argument_group('optional arguments')

    required_suppr_reads.add_argument("--fastq", "-fq", 
    required=True, 
    help="Fastq file"
    )
    
    optional_suppr_reads.add_argument("--qualityMin", "-q", 
    default=9, 
    help="Minimum quality (Phred score) threshold to save a read. Default = 9."
    )
    
    optional_suppr_reads.add_argument("--lengthMin", "-l", 
    default=2000, 
    help="Minimum length threshold to save a read. Default = 2000."
    )
    
    optional_suppr_reads.add_argument("--gzip", 
    required=False, 
    action="store_true", 
    help="fastq.gz"
    )


    parser_suppr_reads.set_defaults(parser_matrix=False, parser_identification=False,parser_length=False,parser_merge_length=False, parser_clean=False, parser_suppr_reads=True)


    #End parser#######

    args = parser.parse_args()
    
    if args != argparse.Namespace(): 
        if args.parser_matrix:
            python_scripts.matrix.main(args)
        elif args.parser_identification:
            python_scripts.identification.main(args)
        elif args.parser_length:
            python_scripts.length.main(args)
        elif args.parser_merge_length:
            python_scripts.merge_length.main(args)
        elif args.parser_clean:
            python_scripts.clean_merge.main(args)
        elif args.parser_suppr_reads:
            python_scripts.suppr_bad_quality_reads.main(args)
    else :
        parser.print_help()

