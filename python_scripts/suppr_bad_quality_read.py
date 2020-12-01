#!/usr/bin/env python3
# coding: utf-8

import argparse,os,statistics,math,gzip

def cli_options():
    "Return parsed cli options"
    parser = argparse.ArgumentParser(description=__doc__)
    required = parser.add_argument_group('Required argument')
    optional = parser.add_argument_group('Optional argument')
    required.add_argument("--fastq", "-fq", required=True, help="Fastq file")
    optional.add_argument("--qualityMin", "-q", default=9, help="Quality threshold to save a read")
    optional.add_argument("--lengthMin", "-l", default=2000, help="Length threshold to save a read")
    optional.add_argument("--gzip", required=False, action="store_true", help="fastq.gz")
    return parser.parse_args()

def iterFastq (file: str, threshold: int, length:int, gz: bool, nb_read_tot: int, general_length: int) -> (tuple):
    '''Lis un fichier fastq '''
    if gz:
        f = gzip.open(file, 'rb')
    else :
        f = open(file,'rb')
        
    for line in f:
        if str(line,'utf-8')[0] == '@':
            header = str(line,'utf-8')
            nb_read_tot +=1
            seq = str(f.readline(),'utf-8')
            f.readline()
            qual = []
            seqQual = str(f.readline(),'utf-8')
            general_length += len(seq)
            if len(seq) > length:
                for char in seqQual:
                    q = ord(char)-33
                    error_rate = 10**(-(q/10))
                    qual.append(error_rate)
                read_quality = -10*(math.log10(statistics.mean(qual)))
                if read_quality >= threshold :
                    yield (header, seq, seqQual)
            yield (nb_read_tot,general_length)
            
    f.close()  

def run (file: str, threshold: float, length:int, gz: bool):
    out = file_to_save (file, threshold)
    print (f"Write reads in : {out}")
    o = open(out,'w')
    nb_read_tot = 0
    nb_read_end = 0
    general_length = 0
    length_sup = 0
    for i in iterFastq(file,threshold,length,gz,nb_read_tot,general_length):
        if len(i) == 2:
            nb_read_tot = i[0]
            general_length = i[1]
        elif len(i) == 3 :
            nb_read_end+=1
            length_sup += len(i[1])
            o.write(f'{i[0]}\n{i[1]}\n+\n{i[2]}\n')
        
    print (f"Reads number : {nb_read_tot}\nReads number after filter : {nb_read_end}")
    print (f"Mean reads length : {int(general_length/nb_read_tot)}\nMean reads length after filter : {int(length_sup/nb_read_end)}")
    o.close()
    
def file_to_save (file: str, threshold: int) -> str:
    out = file.split('.fastq')[0]+"_better_than_"+str(int(threshold))+".fastq"
    return out
    
if __name__ == "__main__":
    os.system("date '+%F -> %T'")
    print ("Start")
    options = cli_options()
    run(options.fastq,float(options.qualityMin),int(options.lengthMin),options.gzip)
    print ("End")
    os.system("date '+%F -> %T'")
