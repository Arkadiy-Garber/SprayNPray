# SprayNPray

See the Wiki Tab for installation instructions and a brief tutorial on how to use this software.


Provide either raw contigs or ORFs (amino acid) in FASTA format.
A reference protein dataset also needs to be provided. Ideally, this should be NCBI's RefSeq or nr database from NCBI, 
which can be downloaded via: wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz

This program will then DIAMOND BLAST each ORF against the reference database
(This takes about 30 minutes for a typical bacteria genome of ~4Mbp).

The program will then output a file in CSV format summarizing the dominant taxanomic groups matched to each contig.
This allows the user to visually inspect the data, seeing what the closest taxonomic group to each contig is.
Additionally, the user can specify a set of criteria (e.g. GC-content, read coverage, coding density, closest taxonomic hits) to re-write the provided contigs into a new FASTA file.

SprayNPray is developed by Arkadiy Garber, Arizona State University, Tempe, Arizona, USA.

There is no official publication for spray-and-pray. If it was useful for your work, you can cite it as: A. Garber 2020: SprayNPray, GitHub repository: https://github.com/Arkadiy-Garber/SprayNPray.

Please also cite dependencies used by spray-and-pray.
