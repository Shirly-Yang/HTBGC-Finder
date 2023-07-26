# HTBGC-Finder
a tool for mining potential horizontally transferred BGCs.

版本(Version)：HTBGC-Finder V1.1.0

更新时间(Update)：2023/7/20

项目主页(Project homepage): https://github.com/SharonYXiao/HTBGC-Finder

## 软件安装(Install)

软件已打包上传至anaconda，下载即可使用。
The software has been packaged and uploaded to anaconda., downloaded and added to the environment variable to make it available.

也可使用git下载，可使用wget或主页中直接下载压缩包，解压后使用。
git clone https://github.com/SharonYXiao/HTBGC-Finder


## 使用方法

anaconda下载使用：

conda install -c HTBGC-Finder
HTBGC-Finder -i ./testdata/mags -o ./testdata/result

github安装包下载使用：

python main.py -i ./testdata/mags -o ./testdata/result
python main.py -h 查看所有可用参数及命令

## 脚本 

- 使用说明：分析常用脚本类型
    - .sh为Shell脚本，使用/bin/bash命令执行；
    - .py为Python脚本，使用python执行

- 脚本功能：微生物组数据分析
    - fastani_run.py:去除冗余的MAGs
    - kraken2_run.py:进行物种注释
    - antismash_run.py: 识别BGC
    - bigscape_run.py：聚类BGC为gene cluster family（GCF)
    - genome_download.py:下载相应的参考基因组
    - final_check.py: 通过比对确定潜在的水平转移的BGC

使用此脚本，请引用下文：

If used this script, please cited:



