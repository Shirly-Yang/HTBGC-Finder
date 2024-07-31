# -- coding:UTF-8 --
import os
import argparse
from fastani_run import fastani_run
from antismash_run import antismash_multirun
from kraken2_run import kraken2_multirun
from bigscape_run import bigscape_run
from bigscape_run import GCF_outliers
from genome_download import genome_download
from final_check import final_check


#htbgcfinder -t 1 -i ./testdata/fna -o ./testdata/result
def main():
    parser = argparse.ArgumentParser(description='This is a tool for mining potential horizontally transferred BGC.')
    parser.add_argument('-i','--input', help='Input dir with genomes. ')
    parser.add_argument('-o','--output', help='Output file. default = [Result]')
    parser.add_argument('-v','--version', action='version', version='BGC_HGT v1.1.4',
                            help="Show BGC_HGT's version")
    parser.add_argument('-a','--anicutoff', help='Remove redundant genomes with ANI cutoff as a threshold. default = [99.9]')
    parser.add_argument('-k','--krakencutoff', help='Kraken cutoff. default = [50]')
    parser.add_argument('-t','--thread', help='Thread. default = [8]')
    parser.add_argument('-p','--krakenparallels',help = 'Number of kraken2 parallels.Excessive quantities may lead to memory exhaustion.default = [1]')
    parser.add_argument('-l','--level', help='Assembly level of reference genomes：[all, complete, chromosome, scaffold, contig]. default = [complete]')
    args = parser.parse_args()

    output_file= 'Result'
    ANI_cutoff=99.9
    kraken_cutoff = 50
    kraken_parallels = 1
    thread=8
    level= 'complete'
    if args.input != None:
        input_file = args.input
    if args.output != None:
        output_file= args.output
    if args.anicutoff != None:
        ANI_cutoff= float(args.anicutoff)
    if args.krakencutoff != None:
        kraken_cutoff= float(args.krakencutoff)
    if args.thread != None:
        thread= int(args.thread)
    if args.level != None:
        level= args.level
    if args.krakenparallels != None:
        kraken_parallels= int(args.krakenparallels)
    
    print('input:',input_file,'output:',output_file,'cutoff',ANI_cutoff,'thread',thread)

    if os.path.exists(output_file) == False:
        os.system('mkdir %s' % (output_file))    
    
    redundant_genome=os.path.join(output_file,'redundant_genome')
    antismash_result=os.path.join(output_file,'antismash_result')
    BGC_gbk=os.path.join(output_file,'BGC_gbk')
    GCF=os.path.join(output_file,'GCF')
    kraken_result=os.path.join(output_file,'kraken_result')
    ref_antismash_result=os.path.join(output_file,'ref_antismash_result')
    ref_BGC_gbk=os.path.join(output_file,'ref_BGC_gbk')
    ref_genome=os.path.join(output_file,'ref_genome')
    ref_bigscape=os.path.join(output_file,'ref_bigscape')  
    pathphylotest=os.path.join(output_file,'phylotest')  

    if os.path.exists(ref_antismash_result) == False:
        os.system('mkdir %s' % (ref_antismash_result))
    if os.path.exists(ref_BGC_gbk) == False:
        os.system('mkdir %s' % (ref_BGC_gbk))
    if os.path.exists(ref_genome) == False:
        os.system('mkdir %s' % (ref_genome))
    if os.path.exists(ref_bigscape) == False:
        os.system('mkdir %s' % (ref_bigscape))
    if os.path.exists(pathphylotest) == False:
        os.system('mkdir %s' % (pathphylotest))
        
    ##Step 1. Remove redundant genomes with ANI
    print('fastANI running... ...')
    fastani_run(input_file,ANI_cutoff,thread,redundant_genome)

    ##Step 2. Predict BGC in non-redundant genomes
    print('antismash running... ...')
    antismash_multirun(input_file,antismash_result,BGC_gbk,thread)

    ##Step 3. Identify taxonomy
    print('kraken2 running... ...')
    kraken2_multirun(input_file,kraken_result,kraken_parallels)
    ##Step 4. Groups BGCs into GCFs according to sequence similarity networks and find outlier genome in each GCF
    print('bigscape running... ...')
    Outlier=GCF_outliers(BGC_gbk,GCF,thread)

    print('Outlier:',Outlier)

    final_result=open(os.path.join(output_file,'result.txt'),'w')
    final_result.write('BGC\tDistance\tP-value\tHTBGC\tLength\tType\tRecipient\tPotential Donor\tSimilar BGC\tDetail\n')
    final_result.close()


    for potential_BGC in Outlier:
        print(potential_BGC)

        result_bigscape_i=os.path.join(ref_bigscape,potential_BGC[0],'network_files')
        ref_genome_i=os.path.join(ref_genome,potential_BGC[0])
        ref_gbkout=os.path.join(ref_BGC_gbk,potential_BGC[0])
        download_need='need'
        if os.path.exists(ref_genome_i) == True:
            if len((os.listdir(ref_genome_i)))>1:
                download_need='not_need'
        print(download_need)

        download_result=genome_download(potential_BGC[0],kraken_cutoff,input_file,thread,kraken_result,level,download_need,potential_BGC[4])
        print(download_result)

        if os.path.exists(result_bigscape_i) == False:

            ##Step 6. Predict BGC in Reference genome
            if download_result != 'Not classified to genus':
                ref_bgcout=os.path.join(ref_antismash_result,potential_BGC[0])
                ref_genome_i=os.path.join(ref_genome,potential_BGC[0])
                ref_bigscape_i=os.path.join(ref_bigscape,potential_BGC[0])

                antismash_multirun(ref_genome_i,ref_bgcout,ref_gbkout,thread)
                os.system('cp %s %s' % (os.path.join(BGC_gbk,potential_BGC[1]+'.gbk'),os.path.join(ref_gbkout,potential_BGC[1]+'.gbk')))
                print(os.path.join(BGC_gbk,potential_BGC[1]+'.gbk'),os.path.join(ref_gbkout,potential_BGC[1]+'.gbk'))
                if len(os.listdir(ref_gbkout))>2:
                    GCF_result=bigscape_run(ref_gbkout,ref_bigscape_i,min(len(os.listdir(ref_gbkout)),thread))
                    if GCF_result != 'no network_files':
                        print(GCF_result,potential_BGC,potential_BGC[1])
                        final_check(GCF_result,potential_BGC[1],output_file,download_result[0],input_file,thread,download_result[1])
                    
        else:
            if len(os.listdir(result_bigscape_i)) >0:
                GCF_result=os.path.join(result_bigscape_i,os.listdir(result_bigscape_i)[0],'mix','mix_clustering_c0.30.tsv')
            else:
                GCF_result = 'no network_files'
            if download_result != 'Not classified to genus' and GCF_result != 'no network_files':
                if len(os.listdir(ref_gbkout))>2:
                    final_check(GCF_result,potential_BGC[1],output_file,download_result[0],input_file,thread,download_result[1])
        print(potential_BGC[1],input_file,output_file,thread)

if __name__ == '__main__':
	main()