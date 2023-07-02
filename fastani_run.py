import os

def fastani_run(inputdir,ani_cutoff,thread,redundant_genome):
    output='ANI_result.txt'
    os.system('cd %s && ls > MAGs_list.txt && fastANI --ql MAGs_list.txt --rl MAGs_list.txt -o %s -t %s && rm MAGs_list.txt' % (inputdir,output,thread))
    aniresult=open(os.path.join(inputdir,output),'r')
    lineani=aniresult.readlines()
    MAGs_list=[]
    redundant=[]
    for i in lineani:
        MAG1=i.split('\t')[0]
        if MAG1 not in MAGs_list:
            MAGs_list.append(MAG1)
    for i in lineani:
        MAG1=i.split('\t')[0]
        MAG2=i.split('\t')[1]
        ani=float(i.split('\t')[2])
        if MAG1 != MAG2 and MAG1 in MAGs_list and MAG2 in MAGs_list and ani > ani_cutoff:
            MAGs_list.remove(MAG2)
            redundant.append(MAG2)
    pathredu=redundant_genome
    if os.path.exists(pathredu) == False:
        os.system('mkdir %s' % (pathredu))
    for i in redundant:
        os.system('mv %s %s' % (os.path.join(inputdir,i),os.path.join(pathredu,i)))
    print(MAGs_list)
    #os.system('rm %s' % (os.path.join(inputdir,output)))


#fastani_run('MAGs','ani_result.txt',99)
