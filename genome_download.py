import os

def kraken2taxid(krakenfile,cutoff,kraken_result):
    cutoff=float(cutoff)
    krakendir=kraken_result
    taxo=['D','P','C','O','F','G','S']
    output_tax_level,output_taxid='',''
    pathi=os.path.join(krakendir,krakenfile)
    filei=open(pathi,'r')
    linei=filei.readlines()
    OTU_dic={}
    for i in linei[-1::-1]:
        percent=float(i.split('\t')[0])
        tax_level=i.split('\t')[3]
        taxname=i.split('\t')[5].strip()
        for t in taxo[-1::-1]:
            #print(t)
            if tax_level == t and percent >=cutoff:
                OTU_dic[t]=taxname
    return OTU_dic
    #print(output_tax_level,output_taxid)

def genome_download(krakenfile,cutoff,input_file,thread,kraken_result,complete,download_need):
    HGT_check=kraken_result.replace('kraken_result','ref_genome')
    pathdownload=os.path.join(HGT_check,'',krakenfile)
    if os.path.exists(HGT_check) == False:
        os.system('mkdir %s' % (HGT_check))
    if os.path.exists(pathdownload) == False:
        os.system('mkdir %s' % (pathdownload))
    kraken2taxid_result=kraken2taxid(krakenfile,cutoff,kraken_result)
    print(kraken2taxid_result)
    if 'S' in kraken2taxid_result:
        if download_need=='need':
            print('------Download Species genome...------')
            if kraken2taxid_result['S'] == 'Escherichia coli':
                os.system('ncbi-genome-download bacteria --output-folder %s -M type -F fasta -l %s --flat-output -g "%s" -P -p %s' % (pathdownload,complete,kraken2taxid_result['S'],thread))
            else:
                os.system('ncbi-genome-download bacteria --output-folder %s -F fasta -l %s --flat-output -g "%s" -P -p %s' % (pathdownload,complete,kraken2taxid_result['S'],thread))
    if 'G' in kraken2taxid_result:
        if download_need=='need':
            print('------Download Genus genome...------')
            os.system('ncbi-genome-download bacteria --output-folder %s -F fasta -M type -l %s --flat-output -g "%s" -P -p %s' % (pathdownload,complete,kraken2taxid_result['G'],thread))
        pathtarget1=os.path.join(input_file,krakenfile.replace('.txt',''))
        pathtarget2=os.path.join(HGT_check,krakenfile,krakenfile.replace('.txt',''))
        os.system('cd %s && gunzip  -f *.gz' % (pathdownload))
        os.system('cp %s %s' % (pathtarget1,pathtarget2))
        return kraken2taxid_result
    else:
        return 'Not classified to genus'

    
