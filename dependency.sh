##dependencies
conda create -n HTBGCFinder python=3.7.0
conda activate HTBGCFinder
conda install -y -c bioconda fastani
conda install -y biopython=1.78
conda install -y -c bioconda antismash=6.1.1
conda install -y -c bioconda bigscape
conda install -y -c bioconda ncbi-genome-download 
conda install -y -c bioconda kraken2


##Add kraken2's species database
echo 'export KRAKEN2_DB_PATH="your/path/to/kraken2_db"' >> ~/.bashrc
source ~/.bashrc
conda activate HTBGCFinder
mkdir $KRAKEN2_DB_PATH
cd $KRAKEN2_DB_PATH
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz
tar zxvf minikraken2_v1_8GB_201904.tgz -C $KRAKEN2_DB_PATH
mv minikraken2_v1_8GB/* $KRAKEN2_DB_PATH

##Add PFAM database
which antismash
#your/path/to/HTBGCFinder/bin/antismash
export PFAM_PATH="your/path/to/HTBGCFinder/lib/python3.7/site-packages/antismash/databases/pfam/34.0" >> ~/.bashrc



channels:
  #- https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
  #- https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
  #- https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
  #- https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
  #- https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
  #- https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
  #- bioconda
  - r
  - defaults
  - conda-forge
  #- https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
  #- https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
  - bioconda


