## PGC (Pan-genome analysis of Chinese)



### Description

Characterization of structural variants and de novo assembly of Chinese population by long-read sequencing.


In order to make sure the results are reproduceable, the pipeline is performed using framework **Snakemake** coupled with the environment conducted by **Anoconda**. And the pipeline can be used in other population of the long-read sequencing



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

##### Other packages should be installed with manual:
  * [CD-HIT=v4.8.1](http://weizhongli-lab.org/cd-hit/)
  * [MAKER2=v2.31.1](http://www.yandell-lab.org/software/maker.html)
  * [RepeatMasker=4.0.9](http://www.repeatmasker.org/)
  * [dna-brnn=v0.1-r65](https://github.com/lh3/dna-nn)
  * [InterProScan=v5.40-77.0](https://www.ebi.ac.uk/interpro/)
  * [EIGENSOFT=v7.2.1](https://github.com/DReichLab/EIG)
  * [IGV=v2.8.6](http://software.broadinstitute.org/software/igv/)


The softweres and packages with corresponding versions were used in our pipeline and manualscript.
You can install the packages with different versions if these packages can 


## Installation

### manually install packages

Create the anoconda environment. 

You can install all softwares and packages in one anaconda environment, or in multiple environments.


Creat the environment and install the packages
```
conda create -n NanoSV
```

Most packages can be found in [anaconda](https://anaconda.org/).
For instance, **bedtools**, you can install like this:
```
conda install -c bioconda bedtools
```

### create environment based on config file

Create the environment using the config file:
```
conda env create -f environment.yml
```



## Quick start for the pipeline

In the directory **pipeline** of the package, there are pipeline file **Population.pipeline.py** and the config file **Population.pipeline.yaml**.

For example, we install all packages in the environment with name **NanoSV**,
then we can run the pipeline:
```   
source activate NanoSV
snakemake -p  -s ~/github/NanoHub/pipeline/Population.pipeline.py --configfile ~/github/NanoHub/pipeline/Population.pipeline.yaml -j 24
```




   
   
Functions of the pipeline:

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
* NanoAssembly.rule.py: De novo assembly of genome for each individual




## Contributions:

* [ZhikunWu](https://github.com/ZhikunWu): 598466208@qq.com
* TongLi: 
* ZehangJiang: