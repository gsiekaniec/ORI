# ORI (Oxford nanopore Reads Identification) and the bacterial world

ORI (Oxford nanopore Reads Identification) is a software allowing, from long nanopore reads, to identify the bacterial strains present in a sample. 

There is two sub-parts: (1) the creation of the index containing the bacterial strains and (2) the query of this index with long reads in order to identify the strains. 

The index is based on the structure implemented in [HowDeSBT](https://github.com/medvedevgroup/HowDeSBT) (Robert S Harris and Paul Medvedev, Improved representation of sequence bloom trees, Bioinformatics, btz662) modified in order to use qgrams (word from spaced seeds) instead of kmers.

## Installation

Before starting you must install the strains version of HowDeSBT: instruction are in the [HowDeSBT_strains](https://github.com/gsiekaniec/ORI/tree/master/HowDeSBT_strains) repertory.

### Other installation:

networkx must be installed

## How does it work ?

### (1) First step: create your own index

In repertory containning genomes (fasta format) do:

### Create the bloom filters (.bf) for each genome

<img src="attention.png" alt="warning" width="30"/> Warning: How DeBT and ORI uses the name of the files to facilitate this use it is preferable not to have . or _ in these names.

	path/to/howdesbt makebfQ --k=15 --qgram=../seed/seedfile.txt --bits=0.5G *.fasta

We get the names of the bf (bloom filter) files used to create the tree:

	ls *.bf > leafname

### If you want to cluster close strains (not obligatory): the threshold depending on the proximity of your strains

It is sometimes necessary to launch the command once in order to see in the Hamming distance table which threshold would be the most interesting before relaunching to merging the strains.
    
    path/to/howdesbt distance --list=leafname --threshold=0.0002 --merge 
    path/to/python3 cleanMerge.py -n path/to/leafname -r path/to/repository/with/bf/files

### Create the tree

If you have merged your files :

    ls *.bf > leafname
    
Then :
    
    path/to/howdesbt cluster --list=leafname --tree=union.sbt --nodename=node{number} --cull
    path/to/howdesbt build --HowDe --tree=union.sbt --outtree=howde.sbt

time : 110.48 s  
space max : 497336 kbytes  

### (2) Query the tree with reads (fasta/q files)

	/usr/bin/time -v /home/gsiekani/Documents/Softwares/HowQ/howdesbt queryQ --sort --qgram=/home/gsiekani/Documents/MinION/Strains_identification/sequences/TestIndelSeeds/classicSeed.txt /home/gsiekani/Documents/MinION/Strains_identification/sequences/reads/Mixture_test/3CIRM+JIM/CIRM67_4000_better_than_9.fastq --tree=howde.sbt --threshold=0.5 > ../results/CIRM67_4000_better_than_9.txt

time : 84.04 s  
space max : 71236 kbytes 

	for i in `ls *.fasta`; do echo ${i%.fasta}; done > listname.txt
	/usr/bin/time -v python3 Results_treatment.py -f ../results/CIRM67_4000_better_than_9.txt -l ../results/listname.txt -o ../results/matrice_CIRM67_0.5.tsv

time : 0.53 s  
space max : 30612 kbytes 


### Get the length of the genomes

	python3 getLength.py -g ../genomes -o length.txt

### Identification/Quantification

	python3 ExtractFromMatrice.py -m ../results/matrice_CIRM67_0.5.tsv -f ../results/CIRM67_better_than_9.txt -le length.txt -l ../results/listname.txt -c /home/gsiekani/.local/bin/clingo_compiled_from_repo -t 55 -n 10


