# ORI (Oxford nanopore Reads Identification) and the bacterial world

## Installation

Before starting you must install the strains version of HowDeSBT: instruction are in the HowDeSBT_strains repertory.
networkx must be installed

## First step: creating your own index

In repertory containning genomes (fasta format) do:

### Create the bloom filters (.bf) for each genome
	
	path/to/howdesbt makebfQ --k=15 --qgram=../seed/seedfile.txt --bits=0.5G *.fasta

### If you want to cluster close strains (not obligatory): the threshold depending on the proximity of your strains

It is sometimes necessary to launch the command once in order to see in the Hamming distance table which threshold would be the most interesting before relaunching to merging the strains.
        
    ls *.bf > leafname
    
    path/to/howdesbt distance --list=leafname --threshold=0.0002 --merge 

### Create the tree

    path/to/howdesbt cluster --list=leafname --tree=union.sbt --nodename=node{number} --cull
	
time : 0.07 s  
space max : 9996 kbytes  

	3) /usr/bin/time -v /home/gsiekani/Documents/Softwares/HowQ/howdesbt build --HowDe --tree=union.sbt --outtree=howde.sbt

time : 110.48 s  
space max : 497336 kbytes  

### Query this tree with fasta/q file

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


