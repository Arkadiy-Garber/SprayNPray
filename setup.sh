#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

chmod +x spray-and-pray.py

gzip -d taxmap_slv_ssu_ref_nr_138.1.txt.gz

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

## adding conda channels
conda config --add channels defaults 2> /dev/null
conda config --add channels bioconda 2> /dev/null
conda config --add channels conda-forge 2> /dev/null
conda config --add channels au-eoed 2> /dev/null

## creating GToTree environment and installing dependencies
conda create -n sprayandpray -c r r-tm r-snowballc r-rcolorbrewer r-xml r-argparse r-wordcloud diamond prodigal metabat2 --yes

## activating environment
source activate sprayandpray

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding path to executable script
export PATH=\"$(pwd):"'$PATH'\"" \

echo '#!/bin/sh'" \

export rscripts=\"$(pwd)\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# re-activating environment so variable and PATH changes take effect
source activate sprayandpray


printf "\n        ${GREEN}DONE!${NC}\n\n"
