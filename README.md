# SprayNPray
## Rapid and simple taxonomic profiling of genome and metagenome contigs 
#### Bin validation
#### Contaminant identification and filtering
#### Endosymbiont identification
#### Identification of bacteria-to-eukaryote horizontal gene transfers

### Citing SprayNPray
If you found this software useful to your research please cite as follows:

Garber, AI. 2020: SprayNPray: Rapid and simple taxonomic profiling of genome and metagenome contigs, GitHub repository: https://github.com/Arkadiy-Garber/SprayNPray.

Please also cite the various dependencies used by spray-and-pray: [DIAMOND](https://pubmed.ncbi.nlm.nih.gov/25402007/), [Prodigal](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/), and [Metabat2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/).


Provide either raw contigs or ORFs (amino acid) in FASTA format.
A reference protein dataset also needs to be provided. Ideally, this should be NCBI's RefSeq or nr database from NCBI, 
which can be downloaded via:
    
    wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz

This program will then DIAMOND BLAST each ORF against the reference database
(This takes about 30 minutes for a typical bacteria genome of ~4Mbp).

The program will then output a file in CSV format summarizing the dominant taxanomic groups matched to each contig.
This allows the user to visually inspect the data, seeing what the closest taxonomic group to each contig is.
Additionally, the user can specify a set of criteria (e.g. GC-content, read coverage, coding density, closest taxonomic hits) to re-write the provided contigs into a new FASTA file.


## easy-installation:
  
    git clone https://github.com/Arkadiy-Garber/SprayNPray.git
    cd SprayNPray
    bash setup.sh
    source activate sprayandpray

Do not worry about the dependencies after conda installation. Just enter `source deactivate` when finished using the program.


## Installation without conda (not recommended):

    git clone https://github.com/Arkadiy-Garber/SprayNPray.git
    export PATH=$PATH:$(pwd)

(be sure to put spray-and-pray.py into your $PATH in your bash profile)

### Install Dependencies:

* [DIAMOND](https://github.com/bbuchfink/diamond)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [Metabat](https://bitbucket.org/berkeleylab/metabat)
* [Python3](https://www.python.org/download/releases/3.0/)


## Usage

### quick-start

    spray-and-pray.py -g genomeContigs.fa -out genomeContigs -ref /path/to/nr/nr.faa

### Decontaminating a Pseudomonas assembly

    ParaHunter.sh -a pseudomonas_crude.fa -out pseudomonas_clean.fa -Genus Pseudomonas -species aeruginosa -perc 50 --fa -ref /path/to/nr/nr.faa

In the above command, we require that at least 50% of contig's genes have a top DIAMOND hit to Pseudomonas aruginosa. Contigs matching these parameters will be written to a new FASTA file: pseudomonas_clean.fa


### Pulling out endosymbiont genomes from an assembly of the mealybug Maconellicoccus hirsutus

    ParaHunter.sh -a M_hirsutus_assembly.fa -out endosymbionts.fa -cd 0.5 -L 1000000 -perc 50 --fa -ref /path/to/nr/nr.faa -Domain Bacteria

In the above comomand, we require a gene density of 0.5 genes per kb, maximum length of 1 Mb, top DIAMOND hit to be to a bacterial gene, and that at least 50% of the contig's genes to be of bacterial origin to be written to endosymbionts.fa


### Identifying putative HGTs in an assembly of the mealybug Planococcus citri

    ParaHunter.sh -a P_citri.fa -out putative_hgts.csv --hgt -ref /path/to/nr/nr.faa


