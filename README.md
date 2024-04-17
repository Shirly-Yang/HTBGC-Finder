# HTBGC-Finder
A tool for detecting potential horizontally transferred BGCs.
- Version：HTBGC-Finder V1.1.0
- Update：2023/7/20
- Project homepage: https://github.com/SharonYXiao/HTBGC-Finder

## Install
- dependencies
    - python=3.7.0
    - biopython=1.78
    - antismash=6.1.1
    - fastANI 
    - bigscape
    - ncbi-genome-download
    - kraken2
  - If you encounter problems during the installation process, you can refer to the file [dependencies.txt](https://github.com/SharonYXiao/HTBGC-Finder/blob/master/dependencies.txt). 
  - An exemplary shell script [dependencies.sh](https://github.com/SharonYXiao/HTBGC-Finder/blob/master/dependency.sh) can be found in the root directory.
  This script file should be adjusted if you wish to place some files to different locations or use absolute path names.

- HTBGC-Finder
git clone https://github.com/SharonYXiao/HTBGC-Finder

## Use
- python HTBGC-Finder.py -i ./testdata/mags -o ./testdata/result
- you can use "python HTBGC-Finder.py -h" to view available parameters.

## Script
-  fastani_run.py:Remove redundant MAGs
-  kraken2_run.py:Perform species annotation
-  antismash_run.py: Identify BGC
-  bigscape_run.py：Clustering BGC as gene cluster family (GCF)
-  genome_download.py:Download the corresponding reference genomes
-  final_check.py: Determination of horizontally transferred BGCs by comparison



