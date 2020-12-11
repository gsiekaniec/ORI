![linux](https://github.com/gsiekaniec/ORI/workflows/linux/badge.svg)
![macOS](https://github.com/gsiekaniec/ORI/workflows/macOS/badge.svg)
# <img src="img/ORI.png" alt="ORI" width="3000"/>

ORI (Oxford nanopore Reads Identification) is a software using long nanopore reads to identify bacteria present in a sample at the strain level. 

There are two sub-parts in the ORI program: (1) the creation of the index containing the reference genomes of the interest species and (2) the query of this index with long reads from Nanopore sequencing in order to identify the strain(s). 

The index is based on the structure implemented in [HowDeSBT](https://github.com/medvedevgroup/HowDeSBT) [1] modified in order to use qgrams (word from spaced seeds) instead of kmers.

As previously said, we replaced kmers by qgrams. To get qgrams we use spaced seeds which introduce don’t care positions in kmers that won’t be disturbed by sequencing errors. To select the best seed pattern for classification of long reads we used the [iedera](https://github.com/laurentnoe/iedera) software [2][3]. The best seed for our classification tools seems to be the following: size 15 and weight 13, 111111001111111. This seed is then applied to words of size 15 resulting in qgrams used instead of kmers. 

The preconstructed indexes for *Streptococcus thermophilus* strains are available in the ORI github in the directory: [preconstructed_indexes](https://github.com/gsiekaniec/ORI/tree/master/preconstructed_indexes) directory of ORI.

<sub>1. Robert S Harris and Paul Medvedev, Improved representation of sequence bloom trees, Bioinformatics, btz662 <10.1093/bioinformatics/btz662>

<sub>2. Noe L., Best hits of 11110110111: model-free selection and parameter-free sensitivity calculation of spaced seeds, Algorithms for Molecular Biology, 12(1). 2017 <http://doi.org/10.1186/s13015-017-0092-1> 

<sub>3. Kucherov G., Noe L., Roytberg, M., A unifying framework for seed sensitivity and its application to subset seeds, Journal of Bioinformatics and Computational Biology, 4(2):553-569, 2006 <http://doi.org/10.1142/S0219720006001977> 

----

## Installation

The easiest way to install ORI is through [conda](https://github.com/gsiekaniec/ORI/tree/master/conda).

	conda create -p ori_env
	conda activate ori_env
	conda install -c gsiekaniec -c conda-forge ori

----

## How does it work ?

The seedfile.txt can be found in the [seed](https://github.com/gsiekaniec/ORI/tree/master/seed) directory of ORI.

You can also create the seedfile.txt in your chosen repertory:
	
	cd path/to/repertory
	touch seedfile.txt
	echo "111111001111111" > seedfile.txt

### I) First step: create your own index

<img src="img/attention.png" alt="warning" width="30"/> Caution: reference genomes must be in **.fasta** or **.fna** and in the **same directory**.

If you have the fasta files distributed in several subfolders, you have to redirect them to a single directory.

Repertory containing fasta files (**must be the current directory**) to run the following scripts:

	cd path/to/the/fasta/files/repertory

#### 1) Create the bloom filters (.bf) for each genome

	howdesbt makebfQ --k=15 --qgram=path/to/seedfile.txt --bits=0.25G *.fasta

| Parameters | Description |
|----------|:-------------:|
| --k | length of the seed given in --qgram. You may modify your seedfile.txt and this parameter. |
| --qgram | file containing the used spaced seed. It's a text file containing the seed used on the first line (here 111111001111111). |
| --bits | size of bloom filters. A size that is too small will give too many false positives to be usable and a size that is too large will take up too much space and greatly increase the computation times. |

Then we get the names of the bf (bloom filter) files used to create the tree:

	ls *.bf > leafname

As the last quantification step requires the size of the genomes, it is preferable to calculate it when we have our genomes:

	ORI.py length -g path/to/the/genomes -o path/to/the/output/length.txt
	
| Parameters | Description | Required |
|----------|:-------------:|------:|
| -g/--genomes | path to the repertory containing genome (.fna or .fasta). | Yes |
| -o/--outfile | output file containing length of each genome. | No. Default: length.txt |
	
Now that the bloom filters are created it is no longer necessary to keep the fastas files. **If it is not necessary to keep them**, they can be deleted to save space.
In addition, if the fasta files cannot be completely downloaded on the machine due to lack of space, it is possible to download them little by little and create the filters as you go by deleting the fasta files once in the form of a filter (.bf).

#### 1.5) Option: if you want to cluster closely related strains to identify a fine cluster of strains rather than a mixed list of single strains:
	
	howdesbt distance --list=leafname

| Parameters | Description |
|----------|:-------------:|
| --list | list of the bloom filters names (one per line). |
	

<img src="img/attention.png" alt="warning" width="30"/> Caution: the threshold parameter t depends on the genetic proximity of your strains (based on Hamming distances matrix). To adapt this parameter to your bacterial species, just launch the following command once (`ORI.py threshold_determination`). It gives a figure **threshold.png** as output containing the distribution of the distances between the strains of the index. More generally, if your cluster of strains is too large and gives you to many possibilities of identification, try a lower t value (e.g. i know that the strains number 205, 51 and 55 are really closed on a phylogenic tree, but a bit farther away to strains 54 and 78; if a threshold of 0.0002 (default value) gives you a cluster containing the 5 strains, you can lower to -t 0.0001 to obtained two separated clusters). 
	
	ORI.py threshold -m path/to/hamming_matrix.tsv -t 0.0002
	
| Parameters | Description | Required |
|----------|:-------------:|------:|
| -m/--matrix | path to the hamming distance matrix. It's the output of the first *howdesbt distance* | Yes | 
| -t/--threshold | threshold that we want to set to merge close genomes. Be careful not to set this threshold too high or too low. Floating number between 0 and 1. | No. Default: 0.0002 |
	
The default 0.0002 value is the value used to merge *Streptococcus thermophilus* strains using filters of size 0.5G. This value must be modified in the case of using another species and/or another filter size.
	
Once you have defined your own t value, merge your strains in adapted clusters:
	
	howdesbt distance --list=leafname --threshold=0.0002 --matrix=hamming_matrix.bin --merge
	
| Parameters | Description |
|----------|:-------------:|
| --list | list of the bloom filters names (one per line). |
| --matrix | path to the hamming distance matrix. It's the output of the first *howdesbt distance* |
| --threshold | hamming distance threshold between bloom filter for merging them. Floating number between 0 and 1. |
| --merge | merge maximal cliques ? |

	ORI.py clean_merge -n path/to/leafname -r path/to/repository/with/bf/files -o path/to/the/output/list_number_file.txt
	
| Parameters | Description | Required |
|----------|:-------------:|------:|
| -n/--names | list of the bloom filters names (one per line). | Yes |
| -r/--repository | repertory containing the bloom filter file (.bf). | Yes |
| -o/--outfile | output file. Out file containing one id number, one genome name and the corresponding number of sequence per line. | Yes |

	ls *.bf > leafname_merge

Since the genomes of some strains have been merged, the size of these clusters must also be recalculated:

    ORI.py merge_length -b path/to/leafname_merge -l path/to/length.txt -c path/to/list_number_file.txt -o path/to/the/output/merge_length.txt

| Parameters | Description | Required |
|----------|:-------------:|------:|
| -b/--bflist | list of the bloom filters names (one per line). | Yes |
| -l/--lengthfile | file containing length of each genome. It's the output of *ORI.py length* (default : length.txt). | Yes |
| -c/--correspondance | file containing correspondance between numbers and genomes. It's the output of *ORI.py clean_merge* (list_number_file.txt). | Yes |
| -o/--outfile | output file containing length of each genome or genomes cluster. | No. Default: length_merge.txt |

#### 2) Create the tree

To run these commands you must be in the directory containing the .bf files.
    
    howdesbt cluster --list=leafname/or/leafname_merge --tree=union.sbt --nodename=node{number} --cull
    
| Parameters | Description | 
|----------|:-------------:|
| --list | list of the bloom filters names (one per line). |
| --tree | name for tree toplogy file. |
| --nodename | filename template for internal tree nodes this must contain the substring {number}. |
| --cull | remove nodes from the binary tree; remove those for which saturation of determined is more than 2 standard deviations. |
    
    howdesbt build --howde --tree=union.sbt --outtree=howde.sbt

| Parameters | Description |
|----------|:-------------:|
| --howde | create tree nodes as determined/how, but only store active bits. Create the nodes as rrr-compressed bit vector(s). |
| --tree | name for tree toplogy file. |
| --outtree | name of topology file to write tree consisting of the filters built. |
  
Once the compressed bloom filters have been created, we can delete those that are not compressed:

    ls | grep -Pv 'detbrief.rrr.' | grep '.bf' | xargs rm --

### II) Second step: query the tree with reads (fastq files)

#### 0.5) Deletion of poor quality reads (not mandatory)

In order to facilitate identification it may be wise to remove reads of too poor quality. For this it is possible to use:

	ORI.py suppr_bad_reads -fq path/to/fastq -q min_quality_value -l min_length_value

| Parameters | Description | Required |
|----------|:-------------:|------:|
| -fq/--fastq | fastq file. | Yes |
| -q/--qualityMin | Minimum quality (Phred score) threshold to save a read. | No. Default: 9 |
| -l/--lengthMin | minimum length threshold to save a read. | No. Default: 2000 |
| --gzip | use of compressed fastq: fastq.gz. | No |

#### 1) Query part and construction of the {strains x reads} matrix 

As we show in the publication (coming soon), ORI's identification  is better with 4000 reads than with 16000 due to noise related to sequencing errors. It is therefore advisable to reduce the number of reads with for example:

	head -n 16000 fastq_file > fastq_file_4000_reads.fq

Then we can start the identification:

	howdesbt queryQ --sort --qgram=path/to/seedfile.txt --tree=path/to/howde.sbt --threshold=0.5  fastq_file_4000_reads > path/to/results_howde.txt

| Parameters | Description | 
|----------|:-------------:|
| --sort | sort matched leaves by the number of query qgrams present, and report the number of qgrams present. |
| --qgram | file containing the used spaced seed. It's a text file containing the seed used on the first line (here 111111001111111). |
| --tree | name of the tree toplogy file (howde.sbt). |
| --threshold | fraction of query qgrams that must be present in a leaf to be considered as a match; this must be between 0 and 1. |

	ORI.py matrix -f path/to/results/from/HowDeSBT -l path/to/leafname/or/leafname_merge -o path/to/results/matrix.tsv

| Parameters | Description | Required |
|----------|:-------------:|------:|
| -f/--file | Results file from HowDe output. | Yes |
| -l/--list_name |  list of the bloom filters names (one per line) in the same order than in the matrix. | Yes |
| -out/--output | output {strains x reads} matrice file.. | No. Default: matrice.tsv |

#### 2) Identification/Quantification

	ORI.py identification -m path/to/matrix.tsv -f path/to/results/from/HowDeSBT -le path/to/length.txt/or/merge_length.txt -l path/to/leafname/or/leafname_merge -c path/to/clingo/or/$(which clingo)(with the conda installation)

| Parameters | Description | Required |
|----------|:-------------:|------:|
| -m/--matrix | {strains X reads} matrice file. | Yes |
| -f/--file | results file from Howde output. | Yes |
| -le/--length | file with one genome and is length per line. It's the output of *ORI.py length* or *ORI.py merge_length*. | Yes |
| -l/--listname | list of the bloom filters names (one per line) in the same order than in the matrix. | Yes |
| -c/--clingo_path | clingo path. With a conda installation this path is in $(which clingo). | Yes |
| -o/--output | output results file. | No. Default: out.txt |
| -t/--threshold | Minimum percent value in the matrix for association between reads and species (between 0 and 100). | No. Default: 50 |
| -n/--nbchoices | Only the nbchoices maximum values of a row are considered. Warning, must be less or equal to the number of species. | No. Default: 12 |

The results are not very readable (especially in case of merge of close strains), it is possible to have cleaner results:

	ORI.py beautiful_results -f path/to/results/from/ORI -n path/to/list_number_file.txt --pie_chart	

| Parameters | Description | Required |
|----------|:-------------:|------:|
| -f/--file | results file from ORI. | Yes |
| -n/--number_name_list | file containing correspondance between numbers and genomes. It's the output of *ORI.py clean_merge* (merge_length.txt). | Yes |
| -o/--output | output file. | No. Default: clean_results.txt |
| --pie_chart | create a pie chart of the results in png format. | No |


