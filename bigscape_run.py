import os
import numpy as np
from scipy.spatial.distance import pdist
from sklearn.neighbors import NearestNeighbors


def calculate_outlier_scores(data, k):
    nbrs = NearestNeighbors(n_neighbors=k).fit(data)
    distances, _ = nbrs.kneighbors(data)
    outlier_scores = np.sum(distances[:, -1:], axis=1)
    return outlier_scores

def outliers_caculation_back(dist_matrix):
    mean_dist = np.mean(dist_matrix, axis=1)
    std_dist = np.std(dist_matrix, axis=1)
    if np.std(mean_dist) == 0:
        return 'No outliers',0
    else:
        z_scores = (mean_dist - np.mean(mean_dist)) / np.std(mean_dist)
        threshold = 2
        outliers = np.where(z_scores == np.max(z_scores))[0]
        #print('outliers: ',np.where(z_scores > threshold)[0],outliers,z_scores[outliers][0])
        if z_scores[outliers][0]>threshold:
            #print(dist_matrix,type(dist_matrix),outliers[0],z_scores[outliers][0])
            return outliers[0],z_scores[outliers][0]
        else:
            return 'No outliers',z_scores[outliers][0]

def outliers_caculation(dist_matrix):
    data = np.array(dist_matrix)
    k = len(data)-1
    print(len(data))
    outlier_scores = calculate_outlier_scores(data, k)
    dic_score={}
    sum_score=0
    max_score=(0,0)
    for i, score in enumerate(outlier_scores):
        print(f"样本 {i}: {score}")
        dic_score[i]=score
        sum_score+=float(score)
        if score > max_score[0]:
            max_score=(score,i)
    average=sum_score/(len(dic_score)+1)
    outscore=('No outliers',(max_score[0],average))
    if max_score[0]>average*1.5 and (max_score[0] - average) > 10:
        outscore=((max_score[1],max_score[0]))
    return outscore






def OTU_distance(kraken2_result,MAG1,MAG2):
    OTU_dis=999
    path1=os.path.join(kraken2_result,MAG1)
    path2=os.path.join(kraken2_result,MAG2)
    taxo=['D','P','C','O','F','G','S']
    dict1={}
    dict2={}
    for t in taxo:
        dict1[t]='None1'
        dict2[t]='None2'
    file1=open(path1,'r')
    line1=file1.readlines()
    for l in line1[-1::-1]:
        dict1[l.split('\t')[3]]=l.split('\t')[-1].strip()
    file2=open(path2,'r')
    line2=file2.readlines()
    for l in line2[-1::-1]:
        dict2[l.split('\t')[3]]=l.split('\t')[-1].strip()
    for t in taxo:
        if dict1[t] == dict2[t]:
            OTU_dis=7-taxo.index(t)
            #OTU_dis=OTU_dis**2
            OTU_dis=2**OTU_dis
    return OTU_dis


def bigscape_run(gbkdir,outputdir,thread):
    if os.path.exists(outputdir) == False:
        os.system('mkdir %s' % (outputdir))
    else:
        os.system('rm -r %s' % (os.path.join(outputdir,'network_files')))
    os.system('bigscape.py -i %s -o %s --pfam_dir $PFAM_PATH --mix --no_classify --include_gbk_str .gbk -c %s --cutoffs 0.9' % (gbkdir,outputdir,thread))
    print('bigscape outputdir',outputdir)
    network=os.path.join(outputdir,'network_files')
    if len(os.listdir(network)) >0:
      bigscape_result=os.path.join(network,os.listdir(network)[0],'mix','mix_clustering_c0.90.tsv')      
      return bigscape_result
    else:
      return 'no network_files'

    
    

def GCF_outliers(gbkdir,outputdir,thread):
    kraken_result=gbkdir.replace('BGC_gbk','kraken_result')
    outlier_result=open(gbkdir.replace('BGC_gbk','outlier.txt'),'w')
    if len(os.listdir(gbkdir))>1:
        bigscape_result=bigscape_run(gbkdir,outputdir,thread)
        filei=open(bigscape_result,'r')
        linei=filei.readlines()
        n=1
        BGC_dict={}
        outlier_list=[]
        BGC_pre=[]
        clan_pre=''
        for l in linei:
            BGC=l.split('\t')[0].split('_BGC_')[0]+'.txt'
            clan=l.split('\t')[1].strip()
            BGC_dict[(BGC,clan)]=l.split('\t')[0]
            if clan == clan_pre and BGC not in BGC_pre:
                BGC_pre.append(BGC)
                n+=1
            else:
                if n>2:
                    #print(BGC_pre)
                    dist_matrix=np.full((len(BGC_pre),len(BGC_pre)), 0)
                    for index1 in range(len(BGC_pre)):
                        for index2 in range(len(BGC_pre)):
                            mag1=BGC_pre[index1]
                            mag2=BGC_pre[index2]
                            dist_matrix[index1][index2]=OTU_distance(kraken_result,mag1,mag2)
                    result=outliers_caculation(dist_matrix)
                    if result[0] != "No outliers":
                        outlier=result[0]
                        print('outlier',outlier)
                        out_MAG=BGC_pre[outlier]
                        output=(out_MAG,clan_pre,result[1])
                    else:
                        output=("No outliers",clan_pre,result[1])
                    outlier_result.write(str(output[0])+'\t'+str(output[1])+'\t'+str(output[2])+'\n')
                    outlier_list.append(output)
                    print(dist_matrix,result,output,BGC_pre)
                n=1
                BGC_pre=[]
                BGC_pre.append(BGC)
            clan_pre=clan
        returnlist=[]
        outlier_result.close()
        for outl in outlier_list:
            if outl[0] != "No outliers":
                genome=outl[0]
                out_bgc=BGC_dict[(genome,outl[1])]
                returnlist.append((genome,out_bgc,outl[1],outl[2]))
        return(returnlist)
    else:
        return('Only 1 region')



