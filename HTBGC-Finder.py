import os
import argparse
from fastani_run import fastani_run
from antismash_run import antismash_multirun
from kraken2_run import kraken2_multirun
from bigscape_run import bigscape_run
from bigscape_run import GCF_outliers
from genome_download import genome_download
from final_check import final_check


#python main.py  -t 30 -i /mnt/nfs/5110v5/yangxiao/BGC_HGT/HMP_data/result/finalmags -o /mnt/nfs/5110v5/yangxiao/BGC_HGT/HMP_data/result/main_result
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This is a tool for mining potential horizontally transferred BGCs.')
    parser.add_argument('-i','--input', help='Input dir with genomes. ')
    parser.add_argument('-o','--output', help='Output file. default = [Result]')
    parser.add_argument('-v','--version', action='version', version='BGC_HGT v1.1.0',
                            help="Show BGC_HGT's version")
    parser.add_argument('-c','--cutoff', help='Remove redundant genomes with ANI cutoff as a threshold. default = [99]')
    parser.add_argument('-t','--thread', help='Thread. default = [8]')
    parser.add_argument('-l','--level', help='Assembly level of reference genomes：[all, complete, chromosome, scaffold, contig]. default = [complete]')
    args = parser.parse_args()

    output_file= 'Result'
    ANI_cutoff=99.9
    thread=8
    level= 'complete'
    
    if args.input != None:
        input_file = args.input
    if args.output != None:
        output_file= args.output
    if args.cutoff != None:
        ANI_cutoff= float(args.cutoff)
    if args.thread != None:
        thread= int(args.thread)
    if args.level != None:
        level= args.level
    
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

    if os.path.exists(ref_antismash_result) == False:
        os.system('mkdir %s' % (ref_antismash_result))
    if os.path.exists(ref_BGC_gbk) == False:
        os.system('mkdir %s' % (ref_BGC_gbk))
    if os.path.exists(ref_genome) == False:
        os.system('mkdir %s' % (ref_genome))
    if os.path.exists(ref_bigscape) == False:
        os.system('mkdir %s' % (ref_bigscape))
        
    ##Step 1. Remove redundant genomes with ANI
    fastani_run(input_file,ANI_cutoff,thread,redundant_genome)

    ##Step 2. Predict BGC in non-redundant genomes
    antismash_multirun(input_file,antismash_result,BGC_gbk,thread)

    ##Step 3. Identify taxonomy
    kraken2_multirun(input_file,kraken_result,thread)

    ##Step 4. Groups BGCs into GCFs according to sequence similarity networks and find outlier genome in each GCF
    Outlier=GCF_outliers(BGC_gbk,GCF,thread)

    print(Outlier)
    final_result=open(os.path.join(output_file,'final.txt'),'w')
    
    #final_result.write(str(len(Outlier))+'\n')
    final_result.close()
    #Outlier=[('3300011580_2.fna.txt', '3300011580_2.fna_BGC_c00008_3300011...region003', '21', 1.7320508075688772), ('MARD_SAMN05375008_REFG_MMP05375008.fna.txt', 'MARD_SAMN05375008_REFG_MMP05375008.fna_BGC_c00002_MARD_MM...region001', '46', 1.7320508075688772), ('MARD_SAMN05375008_REFG_MMP05375008.fna.txt', 'MARD_SAMN05375008_REFG_MMP05375008.fna_BGC_c00003_MARD_MM...region001', '86', 1.7320508075688772), ('MARD_SAMN05375008_REFG_MMP05375008.fna.txt', 'MARD_SAMN05375008_REFG_MMP05375008.fna_BGC_c00006_MARD_MM...region001', '87', 1.7320508075688772)]
    for potential_BGC in Outlier:
        print(potential_BGC)
        result_bigscape_i=os.path.join(ref_bigscape,potential_BGC[0],'network_files')
        ref_gbkout=os.path.join(ref_BGC_gbk,potential_BGC[0])
        if os.path.exists(result_bigscape_i) == False:
            download_need='need'
            ##Step 5. Download outliers genome
            download_result=genome_download(potential_BGC[0],10,input_file,thread,kraken_result,level,download_need)

            ##Step 6. Predict BGC in Reference genome
            if download_result != 'Not classified to genus':
                ref_bgcout=os.path.join(ref_antismash_result,potential_BGC[0])
                ref_genome_i=os.path.join(ref_genome,potential_BGC[0])
                ref_bigscape_i=os.path.join(ref_bigscape,potential_BGC[0])

                antismash_multirun(ref_genome_i,ref_bgcout,ref_gbkout,thread)
                os.system('cp %s %s' % (os.path.join(BGC_gbk,potential_BGC[1]+'.gbk'),os.path.join(ref_gbkout,potential_BGC[1]+'.gbk')))
                print(os.path.join(BGC_gbk,potential_BGC[1]+'.gbk'),os.path.join(ref_gbkout,potential_BGC[1]+'.gbk'))
                if len(os.listdir(ref_gbkout))>2:
                    GCF_result=bigscape_run(ref_gbkout,ref_bigscape_i,10)
                    if GCF_result != 'no network_files':
                        print(GCF_result,potential_BGC,potential_BGC[1])
                        final_check(GCF_result,potential_BGC[1],output_file,download_result)
        else:
            download_need='not_need'
            download_result=genome_download(potential_BGC[0],10,input_file,thread,kraken_result,level,download_need)
            if len(os.listdir(result_bigscape_i)) >0:
                GCF_result=os.path.join(result_bigscape_i,os.listdir(result_bigscape_i)[0],'mix','mix_clustering_c0.90.tsv')
            else:
                GCF_result = 'no network_files'
            if download_result != 'Not classified to genus':
                if len(os.listdir(ref_gbkout))>2:
                    final_check(GCF_result,potential_BGC[1],output_file,download_result)

