# 20.440_Project
# Transcriptional landscapes of eRNA in cancer and their comparative analysis with Model Organisms

<h2>Overview</h2>

At a high level, what does your repo contain / do? 
Include citations for any scripts your team did not write that are 
included directly in the repo and any non-standard methods that your 
analysis is based on. 

This repo contains data files and scripts to study transcriptional landscapes of enhancer RNAs (eRNAs) in cancer and in mouse models. Specifically, we investigated eRNA expressions in tumors and normal tissues and identified potential diagnostic and therapeutic eRNA targets. Furthermore, we identified hits with conserved expression levels in mouse models which can lead to meaningful cancer studies on homologous human eRNA targets.

<h2>Data</h2>

At a high level, how was the data generated? 
If it’s too large to upload to your GitHub, where can it be accessed?
Include citations, if any.

Please access the data at (https://www.dropbox.com/sh/rj1o071ufja6c04/AADhtQBcz6WH2l00BzyUl7RQa?dl=0).

It was sourced from the following locations:
Human eRNA sequences and expression levels were exported from the Human Enhancer RNA Atlas (HeRA) (Zhang et al, 2021)
> https://hanlab.tamhsc.edu/HeRA
Human eRNA expression data in cancers were exported from the Enhancer RNA in Cancer (eRic) (Zhang et al., 2019)
> https://hanlab.uth.edu/eRic/
Mouse eRNA sequences and expression levels were exported from the Animal eRNA Database (Animal-eRNAdb) (Jin et al., 2022)
> http://gong_lab.hzau.edu.cn/Animal-eRNAdb/


<h2>Folder structure </h2>

At a high level, what is in each folder and subfolder?

Folder `script` contains a python script for the project and a R script used to convert .rds files from original sources to .csv files in the data folder on Dropbox. \
Folder `table` contains .csv files generated from the python script.\
Folder `figure` contains figures generated from the python script.\

<h2>Version </h2>
The Python script was written using Python 3.9.7 and the R script was written using R version 4.1.2.


<h2>Installation </h2>

How do I run your code?

To run the script, download the data folder and the script folder. Create empty directories ‘table’ and ‘figure’, and run the python script in terminal
```
python <path to script>  
```

What software and package versions do I need to install?

Package dependencies include
> bioinfokit==2.0.8 \
> matplotlib==3.4.3 \
> numpy==1.20.3 \
> pandas==1.3.4 \
> scikit_learn==1.0.2 \
> scipy==1.7.1 \
> statsmodels==0.12.2 \
> venn==0.1.3

All packages can be installed using pip for Python 3.
Example to install a package
```
pip install matplotlib==3.4.3
```


