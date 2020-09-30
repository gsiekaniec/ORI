# ORI (Oxford nanopore Reads Identification)
# and the bacterial world

## Example commandes

### Create the bloom filter (bf) files
	
	/usr/bin/time -v /home/gsiekani/Documents/Softwares/HowQ/howdesbt makebfQ --k=15 --qgram=/home/gsiekani/Documents/MinION/Strains_identification/sequences/TestIndelSeeds/classicSeed.txt --bits=0.5G *.fasta

time : 52.22 s  
space max : 73124 kbytes  

### Create the tree

	1) ls *.bf > leafname
	
	2) /usr/bin/time -v ~/Documents/Softwares/HowQ/howdesbt cluster --list=leafname --tree=union.sbt --nodename=node{number} --cull
	
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

### Identification/Quantification

	python3 getLength -g ../genomes -o length.txt

	python3 ExtractFromMatrice.py -m ../results/matrice_CIRM67_0.5.tsv -f ../results/CIRM67_better_than_9.txt -le length.txt -l ../results/listname.txt -c /home/gsiekani/.local/bin/clingo_compiled_from_repo -t 55 -n 10

