
# -- coding:UTF-8 --
import os
from multiprocessing import Pool

def kraken2_singlerun(parameter):
    genomei=parameter[0]
    resulti=parameter[1]
    os.system('kraken2 %s --db $KRAKEN2_DB_PATH --report %s.txt --quick --thread 1  > kraken2.log 2>&1' % (genomei,resulti))

def kraken2_multirun(inputdir,kraken_result,thread):
    outputdir=kraken_result
    parameter=[]
    if os.path.exists(outputdir) == False:
        os.system('mkdir %s' % (outputdir))
    for i in os.listdir(inputdir):
        genomei=os.path.join(inputdir,i)
        resulti=os.path.join(outputdir,i)
        para=(genomei,resulti)
        if os.path.exists(resulti+'.txt')== False:
            parameter.append(para)        
    pool=Pool(processes=thread)
    pool.map(kraken2_singlerun,parameter)
    
    
def kraken2_run(inputdir,kraken_result,thread):
    outputdir=kraken_result
    if os.path.exists(outputdir) == False:
        os.system('mkdir %s' % (outputdir))
    for i in os.listdir(inputdir):
        genomei=os.path.join(inputdir,i)
        os.system('kraken2 %s --db $KRAKEN2_DB_PATH --report %s.txt --quick --thread %s' % (genomei,os.path.join(outputdir,i),thread))

