#!/bin/bash

../../HowDeSBT_strains/howdesbt makebfQ --k=15 --qgram=../../seed/seedfile.txt --bits=0.04G *.fasta

diff ./EPS_test.fasta.bf ./EPS.fasta.bf
diff ./JIM8232_test.fasta.bf ./JIM8232.fasta.bf
diff ./CIRM65_test.fasta.bf ./CIRM65.fasta.bf
diff ./CIRM67_test.fasta.bf ./CIRM67.fasta.bf
