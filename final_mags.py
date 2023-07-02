import os

pathbinning='/mnt/nfs/5110v5/yangxiao/BGC_HGT/HMP_data/result_metaspades'
pathi=os.path.join(pathbinning,'bin_refinement')
pathfinalmags='/mnt/nfs/5110v5/yangxiao/BGC_HGT/HMP_data/HMP_metaspades/mags'
if os.path.exists(pathfinalmags)== False:
   os.system('mkdir %s' % (pathfinalmags))
#pathsrrna=os.path.join(pathbinning,'srrna.txt')
#srrna=open('/mnt/nfs/5110v5/yangxiao/BGC_HGT/HMP_data/result/srrna.txt',"w")
count=0
for i in os.listdir(pathi):
  pathbins=os.path.join(pathi,i,'metawrap_50_10_bins')
  if os.path.exists(pathbins):
    for m in os.listdir(pathbins):
      #bins=i+m
      copy1=os.path.join(pathbins,m)
      copy2=os.path.join(pathfinalmags,i+'_'+m)
      #os.system('cp %s %s' % (copy1,copy2))
      count+=1
  else:
    srrna=[]
    srrna.append(pathbins)
    print(srrna)
print(count)
    
    
