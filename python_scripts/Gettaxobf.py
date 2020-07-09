import argparse, os, time
from ete3 import NCBITaxa

ncbi = NCBITaxa()

def cli_options():
    "Return parsed cli options"
    parser = argparse.ArgumentParser(description=__doc__)
    required = parser.add_argument_group('Required argument')
    optional = parser.add_argument_group('Optional argument')
    required.add_argument('--high-taxonomy', '-htax', type=int, dest='htax' , required=True, help='Taxonomy max level. ex : Lactobacillales -> 186826')
    optional.add_argument('--low-taxonomy', '-ltax', type=str, dest='ltax', choices=['phylum', 'class', 'order', 'family', 'genus'], default='genus', required=False, help='Taxonomy min level')
    return parser.parse_args()

def getTaxoList (hightaxo,lowtaxo):
    ''' Among the hightaxo gives all the lowtaxo '''
    descendants = ncbi.get_descendant_taxa(hightaxo, rank_limit=lowtaxo)
    taxa = set()
    for i in descendants:
        if ncbi.get_rank([i])[i] == lowtaxo:
            taxa.add(i)
    rest = set(ncbi.get_descendant_taxa(hightaxo))
    return taxa,rest
    
def downloadNcbiFile():
    ''' Download the assembly_summary_refseq.txt containing all the refseq assembly summary '''
    os.system('wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt')

def createDir(path : str):
    ''' Create the directory if it don't exist'''
    try:
        os.mkdir(path)
    except FileExistsError:
        print ('Warning repertory exist !!!\n')
        pass

def indexSummary():
    '''Create an index of the refseq summary'''
    lineDict = {}
    with open('./assembly_summary_refseq.txt','r') as f:
        l = f.tell()
        line = f.readline()
        while line != '':
            if line[0] != '#':
                line = line.strip().split('\t')
                if int(line[5]) not in lineDict.keys():
                    lineDict[int(line[5])]=[l]
                else:
                    lineDict[int(line[5])].append(l)
            l = f.tell()
            line = f.readline()
    return lineDict

def download (f,lineSummary,tax,genomeid,getComplete):
    '''Download and rename the genomeid file(s)'''
    for pos in lineSummary[genomeid]:
        f.seek(pos,0)
        line = f.readline()
        
        # Get ftp for this strain
        
        ftp = 'ftp'+str(line.split('ftp')[-1].split('\t')[0])
        end = ftp.split('/')[-1]
        
        #if line.split('\t')[11] == 'Complete Genome' :  # If you want only  complete assemblies (not contig or scaffold)
        if line.split('\t')[11] != 'Contig' :
            # Creation of the directory for the low taxonomy
            
            if not os.path.exists('./'+str(tax)):
                createDir('./'+str(tax))
                
            # Download the file if it don't exist
            genomePath = './'+str(tax)+'/'+end+'_genomic.fna.gz'
            
            if line.split('\t')[11] == 'Complete Genome':
                getComplete.add(str(genomeid)+'-'+str(end)+'.fna.gz')
            
            if os.path.exists('./'+str(tax)+'/'+str(genomeid)+'-'+str(end)+'.fna.gz'):
                pass
            else :
                time.sleep(0.5)
                os.system('wget -P ./'+str(tax)+' '+ftp+'/'+end+'_genomic.fna.gz')
                try :
                    os.rename(str(genomePath),'./'+str(tax)+'/'+str(genomeid)+'-'+str(end)+'.fna.gz')
                except FileNotFoundError :
                    compt = 0
                    while not os.path.exists(genomePath):
                        compt += 1
                        time.sleep(1)
                        os.system('wget -P ./'+str(tax)+' '+ftp+'/'+end+'_genomic.fna.gz')
                        if compt == 10 :
                            raise RuntimeError ('Impossible to download this genome : '+str(tax)+' !')
                    os.rename(str(genomePath),'./'+str(tax)+'/'+str(genomeid)+'-'+str(end)+'.fna.gz')
        return getComplete

def treeCreation(HOWQ):
    if len(os.listdir('.')) > 4:
        os.system('ls *.bf > leafnames')
        os.system(str(HOWQ)+'howdesbt cluster --list=leafnames --tree=union.sbt --nodename=node{number} --cull')
        os.system(str(HOWQ)+'howdesbt build --HowDe --tree=union.sbt --outtree=howde.sbt')

def suppFiles():
    for file in os.listdir('.'):
        if (not '.rrr.' in file and '.bf' in file) : #or file.endswith(".fna.gz") :
            '''os.remove(file)'''
    if os.path.exists("leafnames"):
        os.remove("leafnames")
    if os.path.exists("list"):
        os.remove("list")
    if os.path.exists("list_complete"):
        os.remove("list_complete")

def suppEmptyRepertory():
    subfolders = [ f.path for f in os.scandir('.') if f.is_dir()] 
    for directory in subfolders:
        if len(os.listdir(directory)) == 0:
            print (directory)
            os.rmdir(directory)

def manageRest(f,lineSummary,genomeid,highTaxo,diff,getComplete):
    getComplete=set()
    lineage = ncbi.get_lineage(genomeid)
    lowtaxid = lineage[int(lineage.index(highTaxo))+diff]
    try:
        getComplete = download (f,lineSummary,lowtaxid,genomeid,getComplete)
    except KeyError:
        print('No genome available for '+str(ncbi.get_taxid_translator([genomeid])[genomeid]))
    return lowtaxid,getComplete

def createBloomFilterTreeRest(t,HOWQ,seed):
    os.chdir ('./'+str(t[0]))
    print ("Creation of the bloom filters for "+str(str(ncbi.get_taxid_translator([t[0]])[t[0]]))+" sequences")
    os.system('echo ../'+str(t[0])+' > list_complete')
    if t[1] == set():
        os.system('ls *.gz >> list_complete')
    else:
        for complete in t[1] :
            with open('./list_complete','a') as li:
                li.write(str(complete)+'\n')
    if not os.path.exists('../'+str(t[0])+'.bf'):
        print('Execute: '+str(HOWQ)+'howdesbt makebfQ --k=15 --qgram='+str(seed)+' --ubits=1G --bits=0.5G --list=list_complete')                    
        os.system(str(HOWQ)+'howdesbt makebfQ --k=15 --qgram='+str(seed)+' --ubits=1G --bits=0.5G --list=list_complete')
        print('rm ./*.bf')
        os.system('rm ./*.bf')
    os.system('echo all > list')
    os.system('ls *.gz >> list')
    print('Execute: '+str(HOWQ)+'howdesbt makebfQ --k=10 --qgram='+str(seed)+' --bits=0.5G --list=list')
    os.system(str(HOWQ)+'howdesbt makebfQ --k=10 --qgram='+str(seed)+' --bits=0.5G --list=list')
    print('rm all.bf')
    os.system('rm all.bf')
    treeCreation(HOWQ)
    suppFiles()
    os.chdir ('..')

def run(highTaxo: int, lowTaxo: str, seed : str, HOWQ: str, diff : int):
    
    # Download the refseq summary if path do not exist and make an index
    
    if not os.path.exists('./assembly_summary_refseq.txt'):
        downloadNcbiFile()
    createDir('./Bacteria')
    lineSummary = indexSummary()    
    
    # Give the list of all low taxonomy in the high taxonomy
    
    listTax,rest = list(getTaxoList(highTaxo,lowTaxo))
    
    # For low taxonomy how to treat the strains
    with open('./assembly_summary_refseq.txt','r') as f:
        os.chdir ('./Bacteria')
        
        # Test if tree (.rrr. files) already exist
        
        rrr = False
        for file in os.listdir('.'):
            if '.rrr.' in file :
                rrr = True
     
        for t in listTax:
            # Set of Complete genome
            getComplete=set()
            print('##\n##',t,'##\n##')
            print (ncbi.get_taxid_translator([t]),'#####################################################"')
            for genomeid in ncbi.get_descendant_taxa(t):
                rest.remove(genomeid)
                print (ncbi.get_taxid_translator([genomeid]),'----------')
                try:
                    if not os.path.exists('./Bacteria/'+str(t)):
                        getComplete = download(f,lineSummary,t,genomeid,getComplete)
                except KeyError:
                    pass
                    #print('No genome available for '+str(ncbi.get_taxid_translator([genomeid])[genomeid]))
            
            # Create bloom filter for strains
            
            if os.path.exists('./'+str(t)):
                if not rrr:
                    os.chdir ('./'+str(t))
                    print ("Creation of the bloom filters for "+str(str(ncbi.get_taxid_translator([t])[t]))+" sequences")
                    os.system('echo ../'+str(t)+' > list_complete')
                    if getComplete == set():
                        os.system('ls *.gz >> list_complete')
                    else:
                        for complete in getComplete :
                            with open('./list_complete','a') as li:
                                li.write(str(complete)+'\n')
                    if not os.path.exists('../'+str(t)+'.bf'):
                        print('Execute: '+str(HOWQ)+'howdesbt makebfQ --k=15 --qgram='+str(seed)+' --ubits=1G --bits=0.5G --list=list_complete')                    
                        os.system(str(HOWQ)+'howdesbt makebfQ --k=15 --qgram='+str(seed)+' --ubits=1G --bits=0.5G --list=list_complete')
                        print('rm ./*.bf')
                        os.system('rm ./*.bf')
                    
                    os.system('echo all > list')
                    os.system('ls *.gz >> list')
                    print('Execute: '+str(HOWQ)+'howdesbt makebfQ --k=15 --qgram='+str(seed)+' --bits=0.5G --list=list')                    
                    os.system(str(HOWQ)+'howdesbt makebfQ --k=15 --qgram='+str(seed)+' --bits=0.5G --list=list')
                    print('rm all.bf')
                    os.system('rm all.bf')
                    treeCreation(HOWQ)
                    suppFiles()
                    os.chdir ('..')
        
        print('\n\n######################################## End of Genome with good taxonomy ########################################\n\n')
        
        # How we manageme high taxonomy genome without low taxonomy
        lowtaxlist = set()
        for genomeid in rest :
            lowtaxid, getComplete = manageRest(f,lineSummary,genomeid,highTaxo,diff,getComplete)
            lowtaxlist.add((lowtaxid,frozenset(getComplete)))
        for t in lowtaxlist:
            if os.path.exists('./'+str(t[0])) and not rrr:
                createBloomFilterTreeRest(t,HOWQ,seed)
            
        # Create bloom filter for low taxonomy
        
        if not rrr:
            treeCreation(HOWQ)
            suppFiles()
        suppEmptyRepertory()
    os.chdir ('..')

def TestTaxonomy (high, low) -> int:
    taxonomy_rank_oder = {'phylum' : 0, 'class' : 1, 'order' : 2, 'family' : 3, 'genus' : 4}
    if taxonomy_rank_oder[ncbi.get_rank([high])[high]] >= taxonomy_rank_oder[low] :
        raise ValueError('The given high taxonomy is below or equal to the given low taxonomy')
    return (taxonomy_rank_oder[low] - taxonomy_rank_oder[ncbi.get_rank([high])[high]])

def main():
    os.system("date '+%F -> %T'")
    print ("Start\n")
    options = cli_options()
    diff = TestTaxonomy(options.htax, options.ltax)
    HOWdeSBTPath = '/home/gsiekani/Documents/Softwares/HowQ/'
    seed = '/home/gsiekani/Documents/MinION/Strains_identification/sequences/TestIndelSeeds/classicSeed.txt'
    run(int(options.htax), options.ltax, seed, HOWdeSBTPath, diff)
    print ("\nEnd")
    os.system("date '+%F -> %T'")

if __name__ == "__main__":
   main()



