#!/usr/bin/env python3
# coding: utf-8

import decimal
import clyngor
from clyngor import solve
from operator import itemgetter
from python_scripts.EMalgo import vraissemblance, expectation, maximization
from pathlib import Path

ASP="""\

parameter(threshold(threshold)).  parameter(nbchoices(nbchoices)).
parameter(margin(margin)).        parameter(highpass(highpass)).
parameter(parsimony(parsimony)).  parameter(rank0(rank0)).

data(J,F,I,W,R):-  readsstrainsrank(F,I,J,W,R); strain(J,S).
%strainread(Strain name,Read number,Rank of the strain in the read)
strainread(S,I,R):- data(J,_,I,_,R), strain(J,S).
strainread(S,I):- strainread(S,I,_).

% Selected subsets (strain names) and  elements of the subset (reads)
% A strain is candidate for identification 
% if it is present with a minimum of highpass qgram proportion 
% in at least one best read (with low ambiguity)
subset(S):- data(J,"best",_,W,_), strain(J,S), W>highpass.
element(X) :- strainread(S,X), subset(S).

nb(subset(N)):- N={subset(X)}.
nb(element(N)):- N={element(X)}.

totrank(N):- N=#sum{R,S,X:strainread(S,X,R), subset(S)}.

% Compute the threshold for the size of marginal subsets 
%(under this threshold, subsets are considered marginal): 
% it's 1% of elements (reads)
onepercentelement((N+1)/100):- nb(element(N)).
% Compute the threshold for the sum of ranks of marginal subsets 
% it's 2% of the total sum of ranks of elements (read ranks)
%twopercenttotrank((N+1)/50):- totrank(N).


% A strain is marginal if it is associated to less than 1% of the reads
marginal(S):- onepercentelement(N), subset(S), {strainread(S,X)}N-1.

marginal_element(X):- marginal(S), strainread(S,X).

remainingsubset(S):- subset(S), not marginal(S).
nb(remainingsubset(N)):- N={remainingsubset(X)}.
remainingelement(X):- subset(S), not marginal(S), strainread(S,X).
nb(remainingelement(N)):- N={remainingelement(X)}.

%Sum of ranks of elements covered by a remaining subset
totrank(S,N):- remainingsubset(S); N=#sum{R,X:strainread(S,X,R)}.


% Choose the subsets, at least one
% The solution is made of non marginal selected strains
%---------------------------------------------------------
1{ choice(S): remainingsubset(S)}.

ranksupport(X,N):- remainingelement(X), N= #min{R: strainread(S,X,R), remainingsubset(S), choice(S)}.

% A chosen subset (strain) S is signed with an element (read) X 
%1) if it is the sole chosen subset covering this element and has rank at most 1
signedsubset1(S,X,R):- choice(S); remainingelement(X); strainread(S,X,R); R<2; {strainread(T,X):choice(T)}1.
%2) if there are exactly two subsets covering the element, 
%   but S has a better support than the other, which is not a signedsubset1 subset
signedsubset2(S,X,R):- choice(S); remainingelement(X); strainread(S,X,R); choice(S2); not signedsubset1(S2,_,_); 
    strainread(S2,X,R2); R<R2; R<2; 2{strainread(T,X):choice(T)}2.

signedsubset(S,X,R):- signedsubset1(S,X,R).
signedsubset(S,X,R):- signedsubset2(S,X,R).
signedsubset(S,X):- signedsubset(S,X,_).

%Statistics on the sum of ranks of elements covered by a chosen signed subset
specranksupport(S,N):- N=#sum{R,X:signedsubset(S,X,R)}; choice(S). 

% Constraints on represented elements (reads) and selected subsets (strains)
%--------------------------------------------------------------------------

% For each non marginal element, some subset must be chosen
:- element(X); not marginal_element(X); not choice(S) : strainread(S,X).

%A subset (strain) is specific of an element (read) if the element is covered by no other chosen strain
specific(S) :- element(X), 1{ choice(Y) : strainread(Y,X) } 1, strainread(S,X).
%A subset (strain) is established if it appears with rank 0 in at least rank0 elements (reads)
established(S):- choice(S), rank0{strainread(S,_,0)}.

%Each selection (strain) must be either established or specific of an element (read)
:- choice(S), not specific(S), not established(S).


% Optimization
%-------------

%Minimize the number of chosen subsets that are not marginal
% and the corresponding total rank support
#minimize{parsimony@5,S:choice(S)}.
#minimize {N@5,X:ranksupport(X,N),N!=#sup}.

%Then Minimize the total specific rank support
#minimize {N@3,S:specranksupport(S,N)}.


%Solution restricted to signed subsets
identified(S):- choice(S); specranksupport(S,N).


%Output
%------
%#show parameter/1.
%#show marginal/1.
%#show nb/1.
#show identified/1.
"""


def get_length(lengthFile: str) -> dict:
    '''Extract the length of each strain from the file'''

    lengths = {}
    try:
        with open(lengthFile, 'r') as f:
            for line in f:
                if len(line.strip().split('\t')) != 3 and len(line.strip().split('\t')) != 2:
                    print(len(line.strip().split('\t')), line)
                    raise ValueError(f'Incorrect \'{lengthFile}\' lengths file !')
                line = line.strip().split('\t')
                genome = line[0]
                length = line[1]
                lengths[genome] = length
    except ValueError:
        raise
    return lengths


def read_species_names(filename: str) -> str:
    '''Get the names of the strains.
    Return each atom species({number},"{name}")'''
    names = []
    try:
        with open(filename) as fd:
            for i, l in enumerate(fd):
                if len(" ".join(l.strip().split('\t')).split(' ')) > 1:
                    raise ValueError(f'Names file {filename} does not contain one strain per line !')
                name = l.strip()
                num = i + 1
                name = name
                names.append(f'strain({num},"{name}").')
    except ValueError:
        raise
    return ('\n'.join(names) + '\n', num)


def read_matrix(filename: str, nbchoices: int, threshold: float, margin: int, nb_strains: int, listname: str) -> str:
    '''Get the strains for each read from the matrix.
    Preprocessing of matrix data, using two parameters, nbchoices and threshold:
    - Deletion of values < threshold;
    - Rounding the weights with three significant digits after the decimal point;
    - Deletion of rows that contain >= nbchoices*(1+ margin/100) values (>= threshold)
    (the read is considered to be ubiquitous and not useful for identification)
    - Flag rows that contain <= nbchoices/(1+ margin/100) values (>= threshold)
    - Sorting the values by descending values and keeping the nbchoices largest ones.
    - Replace values by their rank, rank 0 for the largest value, same rank for ties, rank increasing with decreasing values.

    Return each atom readsstrainsrank(read_number, strain_number, rank).'''

    decimal.getcontext().prec = 3
    l = []  # List of strain ranks for each selected read
    i = 0  # Current number of selected reads
    nbchoices = int(nbchoices)
    mg = 1 + margin / 100

    with open(filename) as fd:
        for line in fd:
            i = i + 1
            sline = line.strip().split(" ")
            if i == 1:
                try:
                    if nb_strains != len(sline):
                        raise ValueError(
                            f'File {filename} does not contains the same  number of strains than the {listname} names file !')
                except ValueError:
                    raise
            vector = [(j + 1, int(decimal.Decimal(w) * 100)) for j, w in enumerate(sline) if
                      (decimal.Decimal(w) * 100) >= int(threshold)]

            if len(vector) > 0 and len(vector) < int(nbchoices * mg):  # read is neither empty or ubiquitous
                if len(vector) <= int(
                        nbchoices / mg):  # compute a flag best or std, depending on the identification ambiguity of reads
                    l.extend(rank_list(nbchoices, "best", i, vector))  # best reads contain few strain choices
                else:
                    l.extend(rank_list(nbchoices, "std", i,
                                       vector))  # std reads contain more strain choices, at most nbchoices
    return ('\n'.join(l) + '\n')


def rank_list(n: int, flag: str, ident: str, x: int) -> str:
    ''' Provides the rank of at most n greatest elements in a list x of strains weights with possible ties
        Returns tuples (flag  for the selected read, id of selected read, index in x (=strain id),
                        weight of the strain element, rank of the strain element in x)'''
    m = len(x)
    mn = min(n, m)
    t = list(range(mn))
    s = sorted(x, key=itemgetter(1), reverse=True)  # sorting by descending values

    t[0] = 0  # the largest value gets  rank 0
    for k in range(mn - 1):  # compute the  rank of values, giving the same rank for ties
        t[k + 1] = t[k] + (s[k + 1][1] != s[k][1])  # the ranks are increasing with values decreasing

    # w is the proportion of qgrams and t[i] the rank of w in the read
    return [f'readsstrainsrank("{flag}",{ident},{j}, {w}, {t[i]}).' for i, (j, w) in
            enumerate(s[:mn])]  # keeps at most n ranks

def treat_strain(strain,dic):
    '''Take a dictionary and a set of strains and add 1 to the set of strains'''
    if frozenset(strain) in dic.keys() :
        dic[frozenset(strain)] += 1
    else:
        dic[frozenset(strain)] = 1

def quantification(file,strains,marginals):
    dic ={"Unknow":0,"Possibly marginal":0,"Possibly error":0}
    with open(file,'r') as f:
        f.readline()
        strain = set()
        marginal = set()
        for line in f:
            line = line.strip()
            if line[0] == '*':
                if line.split(' ')[1] == '0':
                    dic['Unknow'] += 1
                elif strain == set():
                    if marginal == set():
                        dic['Possibly error'] += 1
                    else:                        
                        dic['Possibly marginal'] += 1
                else:
                    treat_strain(strain,dic)
                strain = set()
                marginal = set()
            else:
                s = line.split(' ')[0].split('.fasta')[0]
                if s in strains:
                    strain.add(s)
                elif s in marginals:
                    marginal.add(s)
        if line.split(' ')[1] == '0':
            dic['Unknow'] += 1
        elif strain == set():
            if marginal == set():
                dic['Possibly error'] += 1
            else:                        
                dic['Possibly marginal'] += 1
        else:
            treat_strain(strain,dic)
    #print(dic)

def compute_abundance(out, strains, nb_strains, abundances, reads, length):
    '''Computation of the abundance, this part of the code still needs to be tested and improved'''

    print(f'Start : {abundances}')

    v = vraissemblance(strains, nb_strains, reads, abundances, length)

    for i in range(100):

        # Expectation part
        njs = expectation(strains, nb_strains, reads, abundances)

        # Maximisation part
        ajs = maximization(strains, nb_strains, reads, njs, length)

        modif = 0
        for i in range(nb_strains):
            modif = modif + (abundances[i] - ajs[i])
        abundances = ajs

        v = vraissemblance(strains, nb_strains, reads, abundances, length)

    # Write the results
    with open(out, 'w') as o:
        for i, s in enumerate(strains):
            o.write(f'{s} : {round(abundances[i], 2)}\n')

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def main(args):
    
    ####Tests
    try:
        if not Path(args.matrix).is_file():
            raise FileNotFoundError(f'File \'{args.matrix}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
    
    try:
        if not Path(args.length).is_file():
            raise FileNotFoundError(f'File \'{args.length}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
    
    try:
        if not Path(args.file).is_file():
            raise FileNotFoundError(f'File \'{args.file}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
    
    try:
        if not Path(args.listname).is_file():
            raise FileNotFoundError(f'File \'{args.listname}\' doesn\'t exist !')
    except FileNotFoundError:
        raise
    
    try:
        if not isfloat(args.threshold):
            raise ValueError(f'Threshold \'{args.threshold}\' is not a number !')
        elif float(args.threshold) < 0 or float(args.threshold) > 100:
            raise ValueError(f'Threshold \'{args.threshold}\' is not a valid number ! (must be between 0 and 100 (%))')
    except ValueError:
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
    
    clyngor.CLINGO_BIN_PATH = args.clingo_path
    
    #Get the strains names 
    names,nb_strains = read_species_names(args.listname)
    
    # Margin
    # It is used to filter rows either that have to many null values to be retained
    # or on the opposite the best rows that have a sufficient number of non null values

    # Get the strains for every reads
    strainsreads = read_matrix(args.matrix, args.nbchoices, args.threshold, args.margin, nb_strains, args.listname)

    #Get genomes length
    lengths = get_length(args.length)

    #Initialisation of the constant for the ASP script
    # Initialisation of the constant for the ASP script
    constante = f'#const threshold={args.threshold}. \n#const nbchoices={args.nbchoices}.  \n#const margin={args.margin}.\n#const highpass={args.highpass}. \n#const parsimony=5. \n#const rank0=5.\n'

    #order the ASP script to run in toexec and solve it
    toexec = constante+names+strainsreads+ASP
    try:
        answers = solve(inline=toexec,options='--opt-mode=optN')
    except FileNotFoundError:
        print(f'\nWrong clingo path \'{args.clingo_path}\' !\n')
        raise
    strains_in = set()
    i=0
    models = []
    
    for answer,optimization,optimality,answerNumber in answers.with_answer_number:    
        #print(f'Answer {i+1}',optimization,optimality,answerNumber)
        if optimality:
            models.append(answer)
        i += 1
    
    #If no answer then stop the execution and ask for more reads
    if i == 0:
        print('No strains recognized ! Maybe not enough reads or the strain(s) from the reads is/are not in the index.')
        exit(1)
    
    soluce = set()

    for model in models:
        for atom in model:
            if atom[0] == 'identified':
                strain = atom[1][0][1:-1]
                if '.bf' in strain:
                    strain = strain.split('.bf')[0]
                if '.fasta' in strain:
                    strain = strain.split('.fasta')[0]
                elif '.fna':
                    strain = strain.split('.fna')[0]
                strains_in.add(strain)
                
        soluce.add(frozenset(strains_in))

    strains_in = list(soluce)[-1]
    
    #In case of multiple optima, retain the set with a minimal number of strains
    for s in soluce:
        if len(s) < len(strains_in):
            strains_in = set(s) #strains_in contains the minimum
    
    reads = []
    with open(args.file,'r') as f: #Results file read from HowDeSBT
        f.readline()
        genomes = set()
        for line in f:
            if line[0] == '*':
                reads.append(genomes)
                genomes = set()
            else:
                line = line.strip().split(' ')
                if '.fna' in line[0] or '.fasta' in line[0]:
                    genome = ".".join(line[0].split('.')[:-1])
                else:
                    genome= line[0]
                if '.fna' in line[0] or '.fasta' in line[0]:
                    if genome[-1] == '_':
                        genome = genome[:-1]
                if genome in strains_in:
                    genomes.add(genome)
    reads.append(genomes)
    
    strains = sorted(list(strains_in))
    nb_strains = int(len(strains))
    
    #Start the abundance computation with a uniform distribution of strains
    abondances = []
    for i in range(nb_strains):
        abondances.append(float(1/nb_strains))
        
    #Get the length of each strain in the sample
    longueurs = []
    for s in strains:
        try :
            longueurs.append(float(lengths[s]))
        except KeyError:
            try:
                s = '-'.join(s.split('_'))
                longueurs.append(float(lengths[s]))
            except KeyError:
                print(f'File {args.length} does not contains the same strains than {args.listname}')
                raise
    
    #Compute the abundance for each strain and write it in the output file
    compute_abundance(args.output,strains,nb_strains,abondances,reads,longueurs)
        

if __name__ == "__main__":
    main()
