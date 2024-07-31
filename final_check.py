
# -- coding:UTF-8 --
import os
from genome_download import kraken2taxid
from phylotest import phylotest


def dic2otu(taxdic):
    taxo=['D','P','C','O','F','G','S']
    otu=''
    for i in taxo:
        if i in taxdic:
            otu+=i
            otu+='__'
            otu+=taxdic[i]
            otu+=';'
    return otu
def BGC_len(target_bgc,pathgbk):
    length='NA'
    pathbgc=os.path.join(pathgbk,target_bgc+'.gbk')
    filebgc=open(pathbgc,'r')
    linebgc=filebgc.readlines()
    length=linebgc[0].split(' bp')[0].split(' ')[-1]
    return length


def final_check(tsvfile,target_bgc,output,download_result,input_file,thread,download_result_donor):
    fout=open(os.path.join(output,'result.txt'),'a')
    phylotest_result=phylotest(target_bgc,input_file,output,thread)
    network=os.path.join(output,'GCF','network_files')
    pathgbk=os.path.join(output,'BGC_gbk')
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
    if os.path.exists(tsvfile) == False:
        fout.write('-\tBigscape Error\n')
    else:
        filetsv = open(tsvfile,'r')
        linetsv = filetsv.readlines()
        clan_count={}
        clan_list={}
        target_clan=''
        for i in linetsv[1:]:
            BGC=i.split('\t')[0]
            clan=i.split('\t')[1].strip()
            if clan in clan_count:
                clan_count[clan]+=1
                clan_list[clan].append(BGC)
            else:
                clan_count[clan]=1
                clan_list[clan]=[]
                clan_list[clan].append(BGC)
            if BGC == target_bgc:
                target_clan = clan
        print(target_bgc,target_clan,tsvfile)
        #print(target_bgc,target_clan,clan_count[target_clan])
        filetsv.close()
        f_g_s=''
        f_g_s_donor=''
        print(download_result)
        f_g_s=dic2otu(download_result)
        f_g_s_donor=dic2otu(download_result_donor)
    
        ref_BGC_genome=os.path.join(output,'ref_genome',target_bgc.split('_BGC_')[0]+'.txt')
        if clan_count[target_clan] > 1:
            print(target_bgc,'Not HGT BGC.')
            ishtbgc=False
            length=BGC_len(target_bgc,pathgbk)
            fout.write(str(ishtbgc)+'\t'+length+'\t'+bgctype1+'/'+bgctype2+'\tRecipient:'+f_g_s+'\tPotential Donor:'+f_g_s_donor+'\t'+str(clan_count[target_clan]-1)+' similar BGCs were found in '+str(len(os.listdir(ref_BGC_genome))-1)+' genome in same genus.\t'+tsvfile+'  Family Number: '+target_clan+'\n')
        else:
    
            
            print(target_bgc,'phylo',phylotest_result)
            if phylotest_result == 'True':
                ishtbgc=True
                print (target_bgc,'HGT BGC!')
            else:
                print(target_bgc,'Not HGT BGC.')
                ishtbgc=False
            #if len(os.listdir(ref_BGC_genome)) > 1:
            length=BGC_len(target_bgc,pathgbk)
            fout.write(str(ishtbgc)+'\t'+length+'\t'+bgctype1+'/'+bgctype2+'\tRecipient:'+f_g_s+'\tPotential Donor:'+f_g_s_donor+'\tSimilar BGC was not found in '+str(len(os.listdir(ref_BGC_genome))-1)+' genome in same genus.\n')
    fout.close()
