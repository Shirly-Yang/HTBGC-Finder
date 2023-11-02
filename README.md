# HTBGC-Finder
a tool for mining potential horizontally transferred BGCs.

-Version：HTBGC-Finder V1.1.0

-Update：2023/7/20

-Project homepage: https://github.com/SharonYXiao/HTBGC-Finder

## Install

git clone https://github.com/SharonYXiao/HTBGC-Finder


## Use

-github安装包下载使用：

python main.py -i ./testdata/mags -o ./testdata/result
python main.py -h 查看所有可用参数及命令

## Script

- 使用说明：分析常用脚本类型
    - .sh为Shell脚本，使用/bin/bash命令执行；
    - .py为Python脚本，使用python执行

- 脚本功能：宏基因组数据中识别潜在的水平转移而来的BGCs
    - fastani_run.py:去除冗余的MAGs
    - kraken2_run.py:进行物种注释
    - antismash_run.py: 识别BGC
    - bigscape_run.py：聚类BGC为gene cluster family（GCF)
    - genome_download.py:下载相应的参考基因组
    - final_check.py: 通过比对确定潜在的水平转移的BGC



