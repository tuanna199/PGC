## PGC (Pan-genome analysis of Chinese)



### Description

Characterization of structural variants and *de novo* assembly of Chinese population by long-read sequencing.


In order to make sure the results are reproduceable, the pipeline is performed using framework [**Snakemake**](https://snakemake.readthedocs.io/en/stable/) coupled with the environment conducted by [**Anoconda**](https://www.anaconda.com/). And the pipeline can be used in other population with long-read sequencing.

The pipeline mainly contains two parts:

* Characterization of structural variants and population genetic analysis
* *De novo* genome assembly and non-reference sequence analysis



### Requirements

##### Softwares:
  * assembly-stats=1.0.1
  * bedtools=2.27.1
  * biopython=1.73
  * blast=2.5.0
  * minimap2=2.15
  * mosdepth=0.2.5
  * nanofilt=2.2.0
  * nanoplot=1.20.0
  * nanoqc=0.8.1
  * nanovar=1.3.6
  * nanosv==1.2.4
  * samtools=1.9
  * seqkit=0.10.1
  * snakemake=5.4.2
  * sniffles=1.0.10
  * ucsc-liftover=377
  * kaiju=1.7.3
  * wtdbg2=2.4
  * diamond=0.9.21
  * augustus=3.3


##### Python packages:
  * biopython=1.73
  * gseapy=0.9.15
  * pysam=0.15.2
  * tinyfasta=0.1.0
  * networkx=2.3
  * numpy=1.16.4

##### R packages:
  * r-argparser=0.4
  * r-base=3.5.1
  * r-dplyr=0.7.8
  * r-futile.logger=1.4.3
  * r-ggplot2=2.2.1
  * r-ggpubr=0.2.4
  * r-gridextra=2.3
  * r-qqman=0.1.4
  * r-readr=1.3.1
  * r-reshape2=1.4.3
  * r-scales=1.0.0
  * r-scrime=1.3.5
  * r-venndiagram=1.6.20

##### Other packages should be manually installed:
  * [CD-HIT=v4.8.1](http://weizhongli-lab.org/cd-hit/)
  * [MAKER2=v2.31.1](http://www.yandell-lab.org/software/maker.html)
  * [RepeatMasker=4.0.9](http://www.repeatmasker.org/)
  * [dna-brnn=v0.1-r65](https://github.com/lh3/dna-nn)
  * [InterProScan=v5.40-77.0](https://www.ebi.ac.uk/interpro/)
  * [EIGENSOFT=v7.2.1](https://github.com/DReichLab/EIG)
  * [IGV=v2.8.6](http://software.broadinstitute.org/software/igv/)


The softweres and packages with corresponding versions were tested and used in our pipeline and manuascript.
You can install the packages with different versions if you make sure these packages can run.


## Configure the environment

### manually install packages

* First [install **Annoconda**](https://docs.anaconda.com/anaconda/install/linux/). 

For instance:
```
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
bash Anaconda3-2020.07-Linux-x86_64.sh
```


* You can install all softwares and packages in one anaconda environment, or in multiple environments.
If you change to one new environment, you can type this common:
```
source activate new-environment
```


* Creat the environment and install the packages
```
conda create -n NanoSV
```

* Most packages can be found in [anaconda](https://anaconda.org/).
For instance, **bedtools**, you can install like this:
```
conda install -c bioconda bedtools
```

### create environment based on config file

Create the environment and install the packages using the config file:
```
conda env create -f environment.yml
```







## Quick start for the pipeline


**Just need clone and copy this package to your destination path.**

In the directory **pipeline** of the package, there are pipeline file **Population.pipeline.py** and the config file **Population.pipeline.yaml**.

For example, we install all packages in the environment with name **NanoSV**,
then we can run the pipeline:
```   
source activate NanoSV
snakemake -p  -s ~/github/NanoHub/pipeline/Population.pipeline.py --configfile ~/github/NanoHub/pipeline/Population.pipeline.yaml -j 24
```




   
## Introduction of pipeline

#### Config file
The config file of parameters [**Population.pipeline.yaml**](https://github.com/ZhikunWu/PGC/blob/master/pipeline/Population.pipeline.yaml)

For instance:
```
### Parameters
CondaENV: /home/wuzhikun/anaconda3/envs/NanoSV
PIPENV: /home/wuzhikun/github/PGC
ProjectPath: /home/wuzhikun/Project/Population

### target samples
SAMPLES:
  - CN196
  - CN068

```

First three parameters are most import. 

* **CondaENV** indicates the anoconda environment which contain installed softwares and packages.
* **PIPENV** indicates the path of this pipeline.
* **ProjectPath** indicates the path of project which contain the data for analysis.


#### Pipeline file
The pipeline file [**Population.pipeline.py**](https://github.com/ZhikunWu/PGC/blob/master/pipeline/Population.pipeline.py)

For example:
```
### config parameter
CondaENV = config["CondaENV"]

### rule of pipeline
include: RULE_DIR + '/BaseFunction.rule.py'

### outcomes of pipeline
rule all:
    input:
        IN_PATH + "/QualityControl/Samples_quality_summary.xls",
```


The pipeline file mainly contain three parts:

* Parameter from **.yaml** config file
* Rule used in pipeline
* Control the pipeline and get the results



#### Rule files

The rule file is based on the function module:

* databaseFormat.rule.py: Change format of compared datasets
* BaseFunction.rule.py: Base function of python
* PopQuality.rule.py: Quality control for nanopore reads
* PopMultipleSV.rule.py: SV calling from multiple callers for each sample
* PopSVFilt.rule.py: High-confidence SVs after three filtering steps
* PopFrequency.rule.py: Allele frequency analysis for non-redundant SVs of the population
* populationGeno.rule.py: Genotypes of SVs in population
* PopStatistics.rule.py: Statistics for non-redundant SVs
* PopStructure.rule.py: Population structure analysis based on SVs
* PopAnnotation.rule.py: Annotation for SVs relative to genomic location
* PopOverlap.rule.py: Comparison to other published SV dataasets
* PopMergeSeq.rule.py: Analysis of SVs associated sequences
* NanoAssembly.rule.py: *De novo* assembly of genome for each individual




## Contributions:

* ZhikunWu: wuzhikun88@gmail.com
* TongLi: tli.aioniya@gmail.com
* ZehangJiang: zehang.kong@gmail.com