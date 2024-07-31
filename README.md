# HTBGC-Finder
A tool for identifying potential horizontally transferred BGCs.
- Version：HTBGC-Finder V1.1.4
- Update：2024/07/31
- Project homepage: https://github.com/Shirly-Yang/HTBGC-Finder

## Installation

### conda (recommended)：

- ` conda create -n htbgcfinder python=3.7.0 `
  
- ` conda activate htbgcfinder `

- ` conda install MicrobiomeX::htbgcfinder `

- `echo 'export PFAM_PATH="/your/path/to/miniconda3/envs/htbgcfinder/lib/python3.7/site-packages/antismash/databases/pfam/34.0"' >> ~/.bashrc `

The database of kraken2 is required but not supported by conda. If you have previously built the database, please skip this step and proceed to the next one.
- ` kraken2-build --standard --threads 8 --db /your/path/to/kraken2database `

 
- `echo 'export KRAKEN2_DB_PATH="/your/path/to/kraken2database/standard"'>> ~/.bashrc `
  
- `source ~/.bashrc `
  
Tips：
If you encounter the following error while building the database: `rsync_from_ncbi.pl: unexpected FTP path (new server?)`, 
you can try the following solution: 
change `    ` in the file /your/path/to/miniconda3/envs/htbgcfinder/libexec/rsync_from_ncbi.pl from `   ` to `   `.


### github：
- ` git clone https://github.com/Shirly-Yang/HTBGC-Finder `
- dependencies
    - python=3.7.0
    - biopython=1.78
    - antismash=6.1.1
    - fastANI 
    - bigscape
    - ncbi-genome-download
    - kraken2
  - An exemplary shell script [dependencies.sh](https://github.com/Shirly-Yang/HTBGC-Finder/blob/master/dependency.sh) can be found in the root directory.
  This script file should be adjusted if you wish to place some files to different locations or use absolute path names.

## Usage

### conda：
- `htbgcfinder -i ./testdata/mags -o ./testdata/result`

### github：
- ` python HTBGC-Finder.py -i ./testdata/mags -o ./testdata/result `
you can use `python HTBGC-Finder.py -h` to view available parameters.

## Script
-  HTBGCFinder.py: identify potential horizontally transferred BGCs
-  fastani_run.py: Remove redundant MAGs
-  kraken2_run.py: Perform species annotation
-  antismash_run.py: Identify BGCs
-  bigscape_run.py：Clustering BGC as gene cluster family (GCF)
-  genome_download.py: Download corresponding reference genomes
-  ref_BGC_fna.py: Identify BGCs
-  phylotest.py: Perform phylogenetic analysis
-  final_check.py: identify horizontally transferred BGCs


If you If used this script, please cite the following article：


