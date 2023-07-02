import os
from genome_download import kraken2taxid

def final_check(tsvfile,target_bgc,output,download_result):
    fout=open(os.path.join(output,'final.txt'),'a')
    network=os.path.join(output,'GCF','network_files')
    network_annotation=os.path.join(network,os.listdir(network)[0],'Network_Annotations_Full.tsv')
    annotationfile=open(network_annotation,"r")
    linei=annotationfile.readlines()
    bgctype1=''
    bgctype2=''
    for b in linei:
        if target_bgc in b:
            bgctype1=b.split('\t')[3]
            bgctype2=b.split('\t')[4]
    annotationfile.close()
    filetsv = open(tsvfile,'r')
    linetsv = filetsv.readlines()
    clan_count={}
    target_clan=''
    for i in linetsv[1:]:
        BGC=i.split('\t')[0]
        clan=i.split('\t')[1].strip()
        if clan in clan_count:
            clan_count[clan]+=1
        else:
            clan_count[clan]=1
        if BGC == target_bgc:
            target_clan = clan
    print(target_bgc,target_clan,tsvfile)
    #print(target_bgc,target_clan,clan_count[target_clan])
    filetsv.close()
    if clan_count[target_clan] > 1:
        print(target_bgc,'Not HGT BGC.')
        fout.write(target_bgc+'\tNot HGT BGC\t'+bgctype1+'/'+bgctype2+'\t'+'Homologous BGCs are found in closely related species.\t'+tsvfile+'\tFamily Number: '+target_clan+'\t'+str(clan_count[target_clan]-1)+' homologous BGCs are found.\n')
    else:
        f_g_s=''
        if 'F' in download_result:
            f_g_s+=download_result['F']
            f_g_s+=' '
        if 'G' in download_result:
            f_g_s+=download_result['G']
            f_g_s+=' '
        if 'S' in download_result:
            f_g_s+=download_result['S']
        ref_BGC_genome=os.path.join(output,'ref_genome',target_bgc.split('_BGC_')[0]+'.txt')
        print (target_bgc,'HGT BGC!')

        if len(os.listdir(ref_BGC_genome)) > 1:
            fout.write(target_bgc+'\tIs HGT BGC\t'+bgctype1+'/'+bgctype2+'\t'+f_g_s+' in '+str(len(os.listdir(ref_BGC_genome))-1)+' genome in same genus, similar BGC were not found.\n')
    fout.close()
        
