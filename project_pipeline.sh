#!/bin/bash

# Script Name: project_pipeline.sh
# Description: Important example scripts for reproducing the projects (excluding data preprocessing).
# Author: Shenni Liang
# Date: 2025-04-14
# Version: 1.0
# Usage: ./project_pipeline.sh
# Notes: Install bpnet env before use. Change directory to your local before running.

# Here we train a RA disease kernel size - 50 model and do follow up modisco analysis as an example

# bpnet train <dataspec.yml> <output dir> [optional flags]`
bpnet train /workdir/shl4034/class_project/dlbio-bpnet/dataspec_test.yml /workdir/shl4034/class_project/RA_disease_model_50_kernel  --override='train.epochs=5;DilatedConv1D.conv1_kernel_size=50'

# compute contribution scores
bpnet contrib /workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696 --method=deeplift --memfrac-gpu=1 --contrib-wildcard='*/profile/wn' /workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696/contrib.deeplift.h5
bpnet contrib /workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696 --method=deeplift --shuffle-seq --max-regions 5000 --contrib-wildcard='*/profile/wn' /workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696/contrib.deeplift.null.h5

# exporting files and motif finding
bpnet export-bw /workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696 /workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696/bigwigs/ --contrib-method=deeplift --scale-contribution
bpnet modisco-run /workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696/contrib.deeplift.h5 \
    --null-contrib-file=/workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696/contrib.deeplift.null.h5 \
    --contrib-wildcard=RA_disease/profile/wn --premade=modisco-50k --override='TfModiscoWorkflow.min_metacluster_size=1000' \
    --only-task-regions /workdir/shl4034/class_project/RA_disease_model_50_kernel/19749f53-9581-4c79-868f-becc7270a696/modisco/RA_disease/ \
    --overwrite
# Note that we don't actually have modisco results for the RA k50 model - the last command is just an example illustration of running modisco