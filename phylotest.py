import os
from multiprocessing import Pool 
import numpy as np
from scipy import stats
from scipy.stats import shapiro, wilcoxon, ttest_rel
import statistics
#def gbk2faa(f):

def perform_statistical_test(list1, list2):
    differences = [y - x for x, y in zip(list1, list2)]
    if len(differences) < 3:
        stat_value = sum(differences) / len(differences)
        _, p_value = stats.wilcoxon(list1, list2)
        
    else:
        _, p_value_shapiro = stats.shapiro(differences)
        if p_value_shapiro >= 0.05:
            _, p_value = stats.ttest_rel(list1, list2)
            stat_value = sum(differences) / len(differences)
        else:
             _, p_value = stats.wilcoxon(list1, list2)
             stat_value = statistics.median(differences)
    ishtbgc='False'
    if stat_value>0 and p_value<0.05:
         ishtbgc='True'
    return(stat_value, p_value,ishtbgc)

def paired_difference_test(list1, list2):
    differences = [b - a for a, b in zip(list1, list2)]
    _, p = shapiro(differences)
    if p > 0.05:
        mean = np.mean(differences)
        #std = np.std(differences)
        stat = f"Mean = {mean:.2f}"
        test_method = "Paired T-Test"
        _, p_value = ttest_rel(list1, list2)
    else:
        median = np.median(differences)
        stat = f"Median = {median:.2f}"
        test_method = "Wilcoxon Signed-Rank Test"
        _, p_value = wilcoxon(list1, list2)

    return stat, test_method, p_value

def gbk2faa(fgbk,ffaa):
    ctg_list=[]
    filei=open(fgbk,'r')
    fout=open(ffaa,'w')
    linei=filei.readlines()
    filei.close()
    ctg=''
    seq=''
    count=0
    switch=False
    for i in linei:
        if switch:
            if count == 2:
                if len(seq) > 20:
                    print('ctg',ctg)
                    ctg_list.append(ctg)
                    fout.write('>candidate_'+ctg+'\n')
                    #print('seq',seq)
                    fout.write(seq+'\n')
                ctg=''
                seq=''
                count=0
                switch=False
            if count <2:
                if i.count('"')==0:
                    seq+=i.strip()
                if i.count('"')==1:
                    seq+=i.split('"')[0].strip()
                    count+=1

        if '/locus_tag=' in i:
            ctg=i.split('"')[1]
        if '/translation=' in i:
            if i.count('"')==2:
                seq=i.split('"')[1]
                count=2
            if i.count('"')==1:
                seq=i.split('"')[1].strip()
                switch=True
                count=1
    fout.close()
    return(ctg_list)

def tblastn_run(parameter):
    ffaa=parameter[0]
    fgenome=parameter[1]
    fresult=parameter[2]
    os.system('tblastn -query %s -subject %s -out %s -max_target_seqs 1 -evalue 1e-05 -outfmt "6 qacc sacc qstart qend length pident qcovs evalue sseq"  > tblastn.log 2>&1' % (ffaa,fgenome,fresult))
    #print('tblastn',parameter)

def mafft_run(parameter):
    fastain=parameter[0]
    fastaout=parameter[1]
    thread = parameter[2]
    #print(fastain,fastaout)
    os.system('mafft --quiet --thread %s %s > %s' % (thread,fastain,fastaout))

def iqtree_run(parameter):
    filealn=parameter[0]
    filetree=parameter[1]
    thread = parameter[2]
    seqcount=0
    aln=open(filealn,'r')
    linealn=aln.readlines()
    for i in linealn:
    	if '>' in i:
    		seqcount+=1
    aln.close()
    if seqcount >4:
    	os.system('iqtree -s %s -m MFP --quiet --prefix %s -T AUTO' % (filealn,filetree))
    else:
    	os.system('iqtree -s %s -m MFP --quiet --prefix %s -T AUTO' % (filealn,filetree))

def fasttree_run(parameter):
    filealn=parameter[0]
    filetree=parameter[1]+'.mldist'
    thread = parameter[2]
    seqcount=0
    aln=open(filealn,'r')
    linealn=aln.readlines()
    aln.close()
    for i in linealn:
    	if '>' in i:
    		seqcount+=1
    os.system('fasttree -makematrix %s > %s' % (filealn,filetree))



def mldist(fmldist):
    #print(fmldist)
    filedist=open(fmldist,'r')
    linedist=filedist.readlines()
    filedist.close()
    maxdist=0
    rindex=[]
    mindex=[]
    candist=[]
    for i in range(len(linedist)):        
        if linedist[i].startswith('candidate_'):
            candist=linedist[i].strip().split()
        if linedist[i].startswith('ref_'):
            rindex.append(i)
        if linedist[i].startswith('mag_'):
            mindex.append(i)
        for j in linedist[i].strip().split()[1:]:
            maxdist=max(maxdist,float(j))
    rmin=max(maxdist,1)
    mmin=max(maxdist,1)
    #rmin=maxdist
    #mmin=maxdist
    if len(candist)>0:
        for r in rindex:
            rmin=min(rmin,float(candist[r]))
        for m in mindex:
            mmin=min(mmin,float(candist[m]))
    #print(fmldist,'r',rmin,'m',mmin)
    result=(mmin,rmin)
    return(result)

def phylotest(BGC,inputfile,pathresult,thread):
	mag=BGC.split('_BGC_')[0]
	BGC=BGC+'.gbk'
	pathref=os.path.join(pathresult,'ref_genome',mag+'.txt')
	fgbk=os.path.join(pathref.replace('ref_genome','ref_BGC_gbk'),BGC)
	
	pathphylo=os.path.join(pathresult,'phylotest',BGC.replace('.gbk',''))
	pathblastref=os.path.join(pathphylo,'blastref')
	pathblastmag=os.path.join(pathphylo,'blastmag')
	pathblastfaa=os.path.join(pathphylo,'faa')
	pathblastaln=os.path.join(pathphylo,'aln')
	pathblasttree=os.path.join(pathphylo,'tree')
	ffaa=os.path.join(pathphylo,'candidateBGC.faa')
	if os.path.exists(pathphylo) == False:
	    os.system('mkdir %s' % (pathphylo))
	if os.path.exists(pathblastref) == False:
	    os.system('mkdir %s' % (pathblastref))
	if os.path.exists(pathblastmag) == False:
	    os.system('mkdir %s' % (pathblastmag))
	if os.path.exists(pathblastfaa) == False:
	    os.system('mkdir %s' % (pathblastfaa))
	if os.path.exists(pathblastaln) == False:
	    os.system('mkdir %s' % (pathblastaln))
	mldist_exists=False
	if os.path.exists(pathblasttree) == False:
	    os.system('mkdir %s' % (pathblasttree))
	else:
	    if len(os.listdir(pathblasttree))>0:
	        if os.listdir(pathblasttree)[0].endswith('.mldist'):
	            mldist_exists=True
	if mldist_exists==False:
		ctg_list=gbk2faa(fgbk,ffaa)
		blastlist=[]
		for i in os.listdir(pathref):
		    if '.txt' not in i:
		        fgenome=os.path.join(pathref,i)
		        fresult=os.path.join(pathblastref,i+'_blastresult.txt')
		        blastlist.append((ffaa,fgenome,fresult))
		
		fgcf=open(os.path.join(pathresult,'MAGs_in_GCF.txt'),'r')
		linegcf=fgcf.readlines()
		fgcf.close()
		for i in linegcf:
		    if BGC.replace('.gbk','') in i:
		        isp=i.strip().split('\t')
		        for j in isp:
		            if j != BGC.replace('.gbk',''):
		                gcfmags=j.split('_BGC_')[0].strip()
		                fgenome=os.path.join(inputfile,gcfmags)
		                fresult=os.path.join(pathblastmag,gcfmags+'_blastresult.txt')
		                blastlist.append((ffaa,fgenome,fresult))
		
		with Pool(processes=thread) as pool:
		    pool.map(tblastn_run,blastlist)
		
		print('Generate protein sequences...')
		for i in ctg_list:
		    seq_dict={}
		    for j in os.listdir(pathblastref):
		        fileblastresult=open(os.path.join(pathblastref,j),'r')
		        linej=fileblastresult.readlines()
		        fileblastresult.close()
		        for l in linej:
		            if i+'\t' in l:
		                seq=l.split('\t')[-1].strip()
		                if mag in j:
		                    genenam='candidate_'+j+'_'+i
		                else:
		                    genenam='ref_'+j+'_'+i
		                seq_dict[genenam]=seq
		                break
		    for j in os.listdir(pathblastmag):
		        fileblastresult=open(os.path.join(pathblastmag,j),'r')
		        linej=fileblastresult.readlines()
		        fileblastresult.close()
		        for l in linej:
		            if i+'\t' in l:
		                seq=l.split('\t')[-1].strip()
		                genenam='mag_'+j+'_'+i
		                seq_dict[genenam]=seq
		    #print(seq_dict)
		    foutfaa=open(os.path.join(pathblastfaa,i+'.faa'),'w')
		    for m in seq_dict:
		        foutfaa.write('>'+m+'\n'+seq_dict[m]+'\n')
		    foutfaa.close()
		
		mafftlist=[]
		print('Aligning...')
		for i in os.listdir(pathblastfaa):
		    ctgfaa=os.path.join(pathblastfaa,i)
		    ctgaln=os.path.join(pathblastaln,i.replace('.faa','_aln.faa'))
		    mafftlist.append((ctgfaa,ctgaln,thread))
		with Pool(processes=thread) as pool:
		    pool.map(mafft_run,mafftlist)
		    
		iqtreelist=[]
		print('Constructing phylogenetic tree...')
		for i in os.listdir(pathblastaln):
		    ctgaln=os.path.join(pathblastaln,i)
		    ctgtree=os.path.join(pathblasttree,i.replace('_aln.faa',''))
		    ctgdistfile=ctgtree+'.mldist'
		    if os.path.exists(ctgdistfile) == False:
		        iqtreelist.append((ctgaln,ctgtree,thread))
		with Pool(processes=thread) as pool:
		    pool.map(fasttree_run,iqtreelist)
	
	mindist_mr=[]
	for i in os.listdir(pathblasttree):
	    if i.endswith('.mldist'):
	        filemldist=os.path.join(pathblasttree,i)
	        mindist_mr.append(mldist(filemldist))
	print(mindist_mr)
	
	listm=[]
	listr=[]
	for i in mindist_mr:
	    listm.append(i[0])
	    listr.append(i[1])
	
	stat, p_value,ishtbgc = perform_statistical_test(listm, listr)
	

	print("Statistic:", stat)
	print("P-Value:", p_value)
	print('HTBGC:',ishtbgc)

	fout=open(os.path.join(pathresult,'result.txt'),'a')
	fout.write(BGC+'\t'+str(stat)+'\t'+str(p_value)+'\t')
	fout.close()
	return ishtbgc

