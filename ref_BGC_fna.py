import os



def ref_BGC_fna(pathdir,patho):
    for i in os.listdir(pathdir):
        path2=os.path.join(pathdir,i)
        for m in os.listdir(path2):
            if m.starswith('GCF'):
               path3=os.path.join(pathdir,i,m)
               pathresult=os.path.join(patho,i,m)
               fileresult=open(pathresult,'w')
               for q in path3:
                   if q.startswith('NC') and endswith('gbk'):
                       pathgbk=os.path.join(path3,q)
                       filebgk=open(pathgbk,'r')
                       switch='off'   
                       seq='' 
                       for linei in filegbk.readlines():
                           if i[0:6] == 'ORIGIN': 
                               switch='on'
                           if switch == 'on':
                               i=i.split() 
                               for j in i[1:]: 
                                   seq+=j.upper() 
                                   nam=m+'_'+q
                                   fileresult.write(('>'+nam+'\n'+seq+'\n')
               fileresult.close()
