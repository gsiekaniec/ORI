# ORI (Oxford nanopore Reads Identification) and the bacterial world

ORI (Oxford nanopore Reads Identification) is a software allowing, from long nanopore reads, to identify the bacterial strains present in a sample. 

There are two sub-parts: (1) the creation of the index containing the bacterial strains and (2) the query of this index with long reads in order to identify the strains. 

The index is based on the structure implemented in [HowDeSBT](https://github.com/medvedevgroup/HowDeSBT) (Robert S Harris and Paul Medvedev, Improved representation of sequence bloom trees, Bioinformatics, btz662) modified in order to use qgrams (word from spaced seeds) instead of kmers.

As previously said, in ORI we replaced kmers by using spaced seeds which introduce don’t care positions in kmers that are not disturbed by errors. To select the best seed pattern for classification of long reads we used the [iedera](https://github.com/laurentnoe/iedera) software (Noe L., Best hits of 11110110111: model-free selection and parameter-free sensitivity calculation of spaced seeds, Algorithms for Molecular Biology, 12(1). 2017 <http://doi.org/10.1186/s13015-017-0092-1>). The best seed for our classification tools seems to be the following spaced seed of size 15 and weight 13: 111111001111111. 

## Installation

Before starting you must install the strains version of HowDeSBT: instruction are in the [HowDeSBT_strains](https://github.com/gsiekaniec/ORI/tree/master/HowDeSBT_strains) repertory.

### Other installation:

python3 is required to run ORI

networkx must be installed 

clingo must be installed: some help [here](https://lucas.bourneuf.net/blog/asp-tuto.html#installation-manuelle)

clyngor is also required

## How does it work ?

### First step: create your own index

<img src="img/attention.png" alt="warning" width="30"/> Warning: refrence genomes must be in .fasta or .fna

In repertory containing reference genomes (fasta format) do:

#### Create the bloom filters (.bf) for each genome

<img src="img/attention.png" alt="warning" width="30"/> Warning: HowDeSBT and ORI use the name of the files to facilitate this use, please avoid "." or "_" in these names.

    path/to/howdesbt makebfQ --k=15 --qgram=../seed/seedfile.txt --bits=0.5G *.fasta
	
Now that the bloom filters are created it is no longer necessary to keep the fastas files, **if it is not necessary to keep them**, they can be deleted to save space.
In addition, if the fastas cannot be completely downloaded on the machine due to lack of space, it is possible to download them little by little and create the filters as you go by deleting the fasta files once in the form of a filter (.bf).

Then we get the names of the bf (bloom filter) files used to create the tree:

    ls *.bf > leafname

As the last quantification step requires the size of the genomes, it is preferable to calculate it when we have our genomes:

    path/to/python3 ORI.py length -g path/to/the/genomes -o path/to/the/output/length.txt

#### If you want to cluster close strains (not obligatory): the threshold depending on the proximity of your strains

It is sometimes necessary to launch the command once in order to see in the Hamming distance table which threshold would be the most interesting before relaunching to merging the strains.
    
    path/to/howdesbt distance --list=leafname --threshold=0.0002 --merge 
    path/to/python3 ORI.py clean_merge -n path/to/leafname -r path/to/repository/with/bf/files -o path/to/the/output/list_number_file.txt

    ls *.bf > leafname_merge

Since the genomes of some strains have been merged, the size of these clusters must also be recalculated:

    path/to/python3 ORI.py merge_length -b path/to/leafname_merge -l path/to/length.txt -c path/to/list_number_file.txt -o path/to/the/output/merge_length.txt

#### Create the tree

To run these commands you must be in the directory containing the .bf files.
    
    path/to/howdesbt cluster --list=leafname --tree=union.sbt --nodename=node{number} --cull
    
    #or if you have merged your files
    path/to/howdesbt cluster --list=leafname_merge --tree=union.sbt --nodename=node{number} --cull
    
    path/to/howdesbt build --HowDe --determined,brief --rrr --tree=union.sbt --outtree=howde.sbt
   
Once the compressed bloom filters have been created, we can delete those that are not compressed:

    ls | grep -Pv 'detbrief.rrr.' | grep '.bf' | xargs rm --

### Query the tree with reads (fasta/q files)

In order to facilitate identification it may be wise to remove reads of too poor quality. For this it is possible to use:

	path/to/python3 python_scripts/suppr_bad_quality_read.py -fq path/to/fastq -q min_quality -l min_length

Then we can start the identification:

	path/to/howdesbt queryQ --sort --qgram=path/to/seedfile.txt --tree=path/to/howde.sbt --threshold=0.5  path/to/fastq_file > path/to/results.txt

	path/to/python3 ORI.py matrice -f path/to/results/from/HowDeSBT -l path/to/leafname/or/leafname_merge -o path/to/results/matrix.tsv

#### Identification/Quantification

	python3 ORI.py identification -m path/to/matrix.tsv -f path/to/results/from/HowDeSBT -le path/to/length.txt -l path/to/leafname/or/leafname_merge -c path/to/clingo


