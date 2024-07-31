# -- coding:UTF-8 --
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
    #outlier_scores = calculate_outlier_scores(data, k)
    dic_score={}
    sum_score=0
    max_score=(0,0)
    print('data,:',data)
    for i in range(len(data)):
        score=max(data[i])-min(data[i])
        dic_score[i]=score
        if score > max_score[0]:
            max_score=(score,i)
    average=sum_score/(len(dic_score)+1)
    outscore=[]
    
    print(dist_matrix,max_score)
    #if max_score[0]>average*1.5 and (max_score[0] - average) > 10:
    for i in range(len(data)):
        distlist=[]
        maxindex=i
        for j in range(len(data[i])):
            if i != j:
                distlist.append(data[i][j] - data[i][i])
                if data[i][j] > data[i][maxindex]:
                    maxindex = j
        if min(distlist) > 4:
            rt_tmp=(i,min(distlist),maxindex)
            outscore.append(rt_tmp)
    if outscore==[]:
        outscore=(('No outliers',(max_score[0],average)))
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
        os.system('bigscape -i %s -o %s --pfam_dir $PFAM_PATH --mix --clans-off --no_classify --include_gbk_str .gbk -c %s --cutoffs 0.3  > bigscape.log 2>&1' % (gbkdir,outputdir,thread))
    else:
        if len(os.listdir(outputdir)) == 0:
            os.system('bigscape -i %s -o %s --pfam_dir $PFAM_PATH --mix --clans-off --no_classify --include_gbk_str .gbk -c %s --cutoffs 0.3  > bigscape.log 2>&1' % (gbkdir,outputdir,thread))
    print('bigscape outputdir',outputdir)
    network=os.path.join(outputdir,'network_files')
    if len(os.listdir(network)) >0:
      bigscape_result=os.path.join(network,os.listdir(network)[0],'mix','mix_clustering_c0.30.tsv')      
      return bigscape_result
    else:
      return 'no network_files'

    
    

def GCF_outliers(gbkdir,outputdir,thread):
    kraken_result=gbkdir.replace('BGC_gbk','kraken_result')
    outlier_result=open(gbkdir.replace('BGC_gbk','outlier.txt'),'w')
    GCF_result=open(gbkdir.replace('BGC_gbk','MAGs_in_GCF.txt'),'w')
    if len(os.listdir(gbkdir))>1:
        bigscape_result=bigscape_run(gbkdir,outputdir,thread)
        filei=open(bigscape_result,'r')
        linei=filei.readlines()
        n=1
        BGC_dict={}
        outlier_list=[]
        BGC_pre=[]
        BGC_pre2=[]
        clan_pre=''
        for l in linei:
            BGC=l.split('\t')[0].split('_BGC_')[0]+'.txt'
            BGC2=l.split('\t')[0]
            clan=l.split('\t')[1].strip()
            BGC_dict[(BGC,clan)]=l.split('\t')[0]
            if clan == clan_pre and BGC not in BGC_pre:
                BGC_pre.append(BGC)
                BGC_pre2.append(BGC2)
                n+=1
            else:
                if n>1:
                    #print(BGC_pre)
                    for mags in BGC_pre2:
                        GCF_result.write(mags+'\t')
                    GCF_result.write('\n')
                    dist_matrix=np.full((len(BGC_pre),len(BGC_pre)), 0)
                    for index1 in range(len(BGC_pre)):
                        for index2 in range(len(BGC_pre)):
                            mag1=BGC_pre[index1]
                            mag2=BGC_pre[index2]
                            dist_matrix[index1][index2]=OTU_distance(kraken_result,mag1,mag2)
                    outliers_result=outliers_caculation(dist_matrix)
                    print('result:',outliers_result)
                    if outliers_result[0] != "No outliers":
                        for result in outliers_result:
                            outlier=result[0]
                            donor = result[2]
                            print('outlier',outlier,'donor',donor)
                            out_MAG=BGC_pre[outlier]
                            dor_MAG=BGC_pre[donor]
                            print('out_MAG',out_MAG,'dor_MAG',dor_MAG)
                            output=(out_MAG,clan_pre,result[1],dor_MAG)

                            outlier_result.write(str(output[0])+'\t'+str(output[1])+'\t'+str(output[2])+'\t'+str(n)+'\t'+str(output[3])+'\n')
                            outlier_list.append(output)
                    else:
                        output=("No outliers",clan_pre,outliers_result[1])
                        outlier_result.write(str(output[0])+'\t'+str(output[1])+'\t'+str(output[2])+'\t'+str(n)+'\n')
                        outlier_list.append(output)                                            
                    print(dist_matrix,outliers_result,'output:',output,BGC_pre)
                n=1
                BGC_pre=[]
                BGC_pre.append(BGC)
                BGC_pre2=[]
                BGC_pre2.append(BGC2)
            clan_pre=clan
        returnlist=[]
        outlier_result.close()
        GCF_result.close()
        for outl in outlier_list:
            if outl[0] != "No outliers":
                genome=outl[0]
                out_bgc=BGC_dict[(genome,outl[1])]
                returnlist.append((genome,out_bgc,outl[1],outl[2],outl[3]))
        print('returnlist:',returnlist)
        return(returnlist)
    else:
        return('Only 1 region')



