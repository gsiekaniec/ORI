# ORI

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

