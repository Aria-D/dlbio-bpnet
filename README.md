

# DLBio - 'MethylBPNet' for Autoimmune Disorders [Final Project]
### Authors: Shenni Liang, Aria Dandawate, Tona Gonzalez

_This is the repository for the DLBio final project -- training BPNet models for use with methylation data to study autoimmune disorders._

## Conda Environment:
We use the BPNet conda environment for the training and post-processing of training results: https://github.com/kundajelab/bpnet.

## Methylation Data:

The methylation data used in our project was retrieved from ADEx (https://adex.genyo.es/). The RA and SLE samples were downloaded and used for the project. The Illumina reference files (see section below), were retrieved from the Illumina website for the BeadChip assays. 

## Data Preprocessing:
First, our methylation data must be preprocessed to be compatible with BPNet. The code for preprocessing is located in the `scripts/preprocessing` directory. First, one should change the paths to the directory where the data is stored. This directory should also have paths to the Illumina 450k reference files that have informtation about the loci of each CpG site. Once these paths are corrected, you can just run the script:

<pre> ```python3 preprocess_data.py ``` </pre>

This script will output .bedgraph and .bw files that will be in the output directory specfied inside the script. 

