# -- coding:UTF-8 --
import os
from multiprocessing import Pool

def antismash_run(parameter):
    genomei=parameter[0]
    resulti=parameter[1]
    os.system('antismash %s --genefinding-tool prodigal  --minimal  --output-dir %s --cpus 1 > antismash.log 2>&1' % (genomei,resulti))
    return


def antismash_multirun(inputdir,outputdir,gbkdir,thread):
    if os.path.exists(outputdir) == False:
        os.system('mkdir %s' % (outputdir))
    parameter=[]
    for i in os.listdir(inputdir):
        genomei=os.path.join(inputdir,i)
        resulti=os.path.join(outputdir,i)
        para=(genomei,resulti)
        parameter.append(para)
        
    pool=Pool(processes=thread)
    pool.map(antismash_run,parameter)
    bgcdir=gbkdir
    if os.path.exists(bgcdir) == False:
        os.system('mkdir %s' % (bgcdir))
    for i in os.listdir(outputdir):
        pathi=os.path.join(outputdir,i)
        for j in os.listdir(pathi):
            if 'region' in j and '.gbk' in j:
                pathj=os.path.join(pathi,j)
                outj=os.path.join(bgcdir,i+'_BGC_'+j)
                os.system('cp %s %s' % (pathj,outj))
