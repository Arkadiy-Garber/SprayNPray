# SprayNPray
## Rapid and simple taxonomic profiling of genome and metagenome contigs 
#### - Bin validation
#### - Metagenome profiling
#### - Contaminant identification and filtering
#### - Endosymbiont identification
#### - Identification of bacteria-to-eukaryote horizontal gene transfers

## Citing SprayNPray
If you found this software useful to your research please cite as follows:

Garber, A. I., Armbruster, C. R., Lee, S. E., Cooper, V. S., Bomberger, J. M., & McAllister, S. M. (2021). SprayNPray: user-friendly taxonomic profiling of genome and metagenome contigs. In bioRxiv (p. 2021.07.17.452725). https://doi.org/10.1101/2021.07.17.452725

Please also cite the various dependencies used by spray-and-pray: [DIAMOND](https://pubmed.ncbi.nlm.nih.gov/25402007/), [Prodigal](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/), [Metabat2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/), and [Biopython](https://biopython.org/)


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
    conda activate sprayandpray

Please also place spray-and-pray.py into the PATH variable: e.g. ~/bin/SprayNPray. You can update the PATH variable in a file called .bash_profile or .profile in your home directory.

Do not worry about the dependencies after conda installation. Just enter `conda deactivate` when finished using the program.


## Installation without conda:

    git clone https://github.com/Arkadiy-Garber/SprayNPray.git
    export PATH=$PATH:$(pwd)

(be sure to put spray-and-pray.py into your $PATH in your bash profile)

### Install Dependencies:

* [DIAMOND](https://github.com/bbuchfink/diamond)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [Metabat](https://bitbucket.org/berkeleylab/metabat)
* [Python3](https://www.python.org/download/releases/3.0/)
* [Biopython3](https://biopython.org/)


## Usage

### quick-start

    spray-and-pray.py -g genomeContigs.fa -out genomeContigs -ref /path/to/nr/nr.faa

### Decontaminating a Pseudomonas assembly

    spray-and-pray.py -a pseudomonas_crude.fa -out pseudomonas_clean.fa -Genus Pseudomonas -species aeruginosa -perc 50 --fa -ref /path/to/nr/nr.faa

In the above command, we require that at least 50% of contig's genes have a top DIAMOND hit to Pseudomonas aruginosa. Contigs matching these parameters will be written to a new FASTA file: pseudomonas_clean.fa


### Pulling out endosymbiont genomes from an assembly of the mealybug Maconellicoccus hirsutus

    spray-and-pray.py -a M_hirsutus_assembly.fa -out endosymbionts.fa -cd 0.5 -L 1000000 -perc 50 --fa -ref /path/to/nr/nr.faa -Domain Bacteria

In the above comomand, we require a gene density of 0.5 genes per kb, maximum length of 1 Mb, top DIAMOND hit to be to a bacterial gene, and that at least 50% of the contig's genes to be of bacterial origin to be written to endosymbionts.fa


### Identifying putative HGTs in an assembly of the mealybug Planococcus citri

    spray-and-pray.py -a P_citri.fa -out putative_hgts.csv --hgt -ref /path/to/nr/nr.faa


