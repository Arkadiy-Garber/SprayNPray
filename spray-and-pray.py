#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import numpy as np
import sys
import statistics


def firstNonspace(ls):
    for i in ls:
        if i != "":
            break
    return i


def gc(seq):
    gc = 0
    for bp in seq:
        if bp == "C" or bp == "G":
            gc += 1
    return gc/len(seq)


def Dictparser(Dictionary):
    lowest = float(1000)
    for i in Dictionary:
        if float(Dictionary[i]) < float(lowest):
            lowest = Dictionary[i]
            key = i
    return [i, lowest]


def reverseComplement(seq):
    out = []
    for i in range(len(seq)-1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def ribosome(seq):
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 4:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def ave(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count/len(ls)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count/len(seq)


def reject_outliers(data):
    m = 2
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
    return filtered


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


parser = argparse.ArgumentParser(
    prog="spray-and-pray.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
    

    
    Developed by Arkadiy Garber; University of Montana, Biological Sciences
    Please send comments and inquiries to arkadiy.garber@mso.umt.edu
    ************************************************************************
    '''))

parser.add_argument('-g', type=str, help="Input bin/assembly in FASTA format", default="NA")

parser.add_argument('-o', type=str, help="Input ORFs in FASTA amino acid format", default="NA")

parser.add_argument('-ref', type=str, help="Input reference protein database (recommended: nr). Could be FASTA file or "
                                           "DIAMOND database file (with extension .dmnd)", default="NA")

parser.add_argument('-bam', type=str, help="Input sorted BAM file with coverage info (optional)", default="NA")

parser.add_argument('-out', type=str, help="Basename for output files", default="NA")

parser.add_argument('-lvl', type=str, help="Level of the taxonomic hierarchy to include in the summary file (Domain, Phylum, Class, Genus, species)", default="NA")

parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST", default=1)

parser.add_argument('--makedb', type=str, help="if the DIAMOND database does not already exist "
                                                    "(i.e. file with extension .dmnd), and you would like the program t"
                                               "o run  diamond makedb, provide this flag", const=True, nargs="?")

parser.add_argument('--spades', type=str, help="is this a SPAdes assembly, with the original SPAdes headers? If so, "
                                               "then you can provide this flag, and BinBlaster will summarize using the coverage "
                                               "information provided in the SPAdes headers", const=True, nargs="?")

parser.add_argument('--meta', type=str, help="contigs are from a mixed community of organisms", const=True, nargs="?")

parser.add_argument('--hgt', type=str, help="provide this flag if you'd like the program to output potential HGTs into a separate file. "
                                            "This feature is designed for eukaryotic contigs expected to have HGTs of bacterial origin.", const=True, nargs="?")

parser.add_argument('--fa', type=str, help="write subset of contigs that match user-specified parameters to a separate FASTA file", const=True, nargs="?")

parser.add_argument('-blast', type=str, help="DIAMOND BLAST output file from previous run", default="NA")

parser.add_argument('-Domain', type=str, help="domain expected among hits to provided contigs, to be written to FASTA file (e.g. Bacteria, Archaea, Eukaryota)", default="NA")

parser.add_argument('-Phylum', type=str, help="phylum expected among hits to provided contigs, to be written to FASTA file (e.g. Proteobacteria). "
                                              "If you provide this name, please be sure to also provide the domain name via -domain", default="NA")

parser.add_argument('-Class', type=str, help="class name expected among hits to provided contigs, to be written to FASTA file (e.g. Gammaproteobacteria). "
                                             "If you provide this name, please be sure to also provide the domain and phylum names", default="NA")

parser.add_argument('-Genus', type=str, help="genus name expected among hits to provided contigs, to be written to FASTA file (e.g. Shewanella). "
                                             "If you provide this name, please be sure to also provide the domain, phylum, and class names", default="NA")

parser.add_argument('-species', type=str, help="species name expected among hits to provided contigs, to be written to FASTA file (e.g. oneidensis, coli, etc.). "
                                               "If you provide this name, please be sure to also provide the domain, phylum, class, and genus names", default="NA")

parser.add_argument('-perc', type=float, help="percentage of total hits to the contig that must be to the specified genus/species for writing to FASTA", default=0)

parser.add_argument('-gc', type=float, help="minimum GC-content of contigs to write to FASTA (default = 0)", default=0)

parser.add_argument('-GC', type=float, help="maximum GC-content of contigs to write to FASTA (default = 100)", default=100)

parser.add_argument('-cov', type=float, help="minimum coverage of contigs to write to FASTA (default = 0)", default=0)

parser.add_argument('-COV', type=float, help="maximum coverage of contigs to write to FASTA (default = 100000000)", default=100000000)

parser.add_argument('-cd', type=float, help="minimum coding density (in hits/kb) to write to FASTA (default = 0.25)", default=0.25)

parser.add_argument('-CD', type=float, help="maximum coding density (in hits/kb) to write to FASTA (default = 5)", default=5)

parser.add_argument('-l', type=float, help="minimum length of contig to write to FASTA (default = 1000)", default=1000)

parser.add_argument('-L', type=float, help="maximum length of contig to write to FASTA (default = 100000000)", default=100000000)

parser.add_argument('-aai', type=float, help="minimum average amino acid identity (percent) to reference proteins (default 35)", default=35)

# parser.add_argument('-key', type=str, help="Path to the taxmap_slv_ssu_ref_nr_138.1.txt file, which should be in the repository containing this program", default="NA")

args = parser.parse_args()

print(".")
# checking paramters:
if args.ref == "NA":
    print("Please provide a reference file via -ref")
    raise SystemExit
else:
    print("Reference file: " + args.ref)

if args.lvl != "NA":
    if args.lvl in ["Domain", "Phylum", "Class", "Genus", "species"]:
        print("Taxonomic level included in the summary file: " + args.lvl)
    else:
        print("Invalud taxonomic level provided via -lvl")
        raise SystemExit

if args.blast != "NA":
    print("Provided BLAST output file: " + args.blast)

os.system("which spray-and-pray.py > mainDir.txt")

file = open("mainDir.txt")
location = os.getcwd()
for i in file:
    location = i.rstrip()
location = allButTheLast(location, "/")

silvaFile = location + "/taxmap_slv_ssu_ref_nr_138.1.txt"


if args.fa:
    print("SprayNPray will write a FASTA file with contigs matching user-specified metrics: " + args.out + "-contigs.fa")
    print("SprayNPray will write a FASTA file with contigs not matching user-specified metrics: " + args.out + "-unmatched.contigs.fa\n")
    if args.species != "NA":
        if "NA" in [args.Genus, args.Class, args.Phylum, args.Domain]:
            print("If species name is provided, please provide also the Genus, Class, Phylum, and Domain names")
            raise SystemExit
        else:
            print("species restriction: " + args.species)

    if args.Genus != "NA":
        if "NA" in [args.Class, args.Phylum, args.Domain]:
            print("If Genus name is provided, please provide also the Class, Phylum, and Domain names")
            raise SystemExit
        else:
            print("Genus restriction: " + args.Genus)

    if args.Class != "NA":
        if "NA" in [args.Phylum, args.Domain]:
            print("If Class name is provided, please provide also the Phylum and Domain names")
            raise SystemExit
        else:
            print("Class restriction: " + args.Class)

    if args.Phylum != "NA":
        if "NA" in [args.Domain]:
            print("If Phylum name is provided, please provide also the Domain name")
            raise SystemExit
        else:
            print("Phylum restriction: " + args.Phylum)

    if args.Domain:
        print("Domain restriction: " + args.Domain)


if args.o != "NA":
    file = open(args.o)
    file = fasta2(file)
    if args.makedb:
        print("Running DIAMOND: making DIAMOND database")
        os.system("diamond makedb --in %s --db %s.dmnd" % (args.ref, args.ref))

    print("Running DIAMOND BLAST")
    os.system(
        "diamond blastp --db %s.dmnd --query %s-proteins.faa "
            "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle "
            "--out %s.blast --max-target-seqs 50 --evalue 1E-15 --threads %d --query-cover 50 --subject-cover 50"
        % (args.ref, args.o, args.o, args.t))

    print("Preparing summary: %s" % args.out)

    aaiDict = defaultdict(list)
    blastDict = defaultdict(list)
    blast = open("%s.blast" % args.o)
    for i in blast:
        ls = i.rstrip().split("\t")
        orf = ls[0]
        name = (ls[12])
        name = name.split("]")[0]
        name = name.split("[")[1]
        blastDict[orf].append(name)
        aai = ls[2]
        aaiDict[orf].append(float(aai))

    out = open(args.out, "w")
    out.write(
        "ORF" + "," + "Average_AAI" + "," + "closest_blast_hits" + "\n")
    for i in file.keys():

        hitsList = blastDict[i]
        try:
            AAI = statistics.mean(aaiDict[i])
        except statistics.StatisticsError:
            AAI = "NA"
        out.write(
            i + "," + str(AAI) + ",")

        for j in hitsList:
            try:
                out.write(j + "; ")
            except TypeError:
                pass
        out.write("\n")

    print("Finished!")


else:

    file = open(args.g)
    file = fasta2(file)
    total = 0
    for i in file.keys():
        total += len(file[i])


    if total < 20000:
        if args.meta:
            pass
        else:
            print("looks like there are less than 20000 characters in your provided sequences file. Please re-run the script with the --meta flag")
            raise SystemExit

    if args.blast == "NA":

        print("Running Prodigal: calling ORFs from provided contigs")
        if args.meta:
            os.system("prodigal -i %s -a %s-proteins.faa -p meta" % (args.g, args.g))
        else:
            os.system("prodigal -i %s -a %s-proteins.faa" % (args.g, args.g))

        if args.makedb:
            print("Running Diamond: making DIAMOND BLAST database")
            os.system("diamond makedb --in %s --db %s.dmnd" % (args.ref, args.ref))


        print("Running Diamond BLAST")
        os.system(
            "diamond blastp --db %s.dmnd --query %s-proteins.faa "
            "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle "
            "--out %s.blast --max-target-seqs 50 --evalue 1E-15 --threads %d --query-cover 50 --subject-cover 50"
            % (args.ref, args.g, args.g, args.t))

        blastFile = "%s.blast" % args.g
    else:
        blastFile = args.blast

    if args.bam != "NA":
        print("Extracting coverage information from the provided BAM files")
        os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s" % (args.g, args.bam))

    print("Calculating GC-content")
    gcDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    GC = 0
    total = 0
    for i in file.keys():
        seq = file[i]
        total += len(seq)
        gc = 0
        for bp in seq:
            if bp == "C" or bp == "G":
                GC += 1
                gc += 1
        gcDict[i] = str( float(gc/len(seq)) * 100 )

    print("Preparing summary: %s" % args.out)

    aaiDict = defaultdict(list)
    blastDict = defaultdict(list)
    redunDict = defaultdict(list)
    blast = open(blastFile)
    for i in blast:
        ls = i.rstrip().split("\t")
        contig = allButTheLast(ls[0], "_")
        name = ls[12]
        try:
            name = name.split("]")[0]
            name = name.split("[")[1]
        except IndexError:
            name = "NA"
        aai = ls[2]
        if ls[0] not in redunDict.keys():
            redunDict[ls[0]].append(name)
            blastDict[contig].append(name)
            aaiDict[contig].append(float(aai))

    if args.bam != "NA":
        depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        depth = open("%s.depth" % args.g)
        for i in depth:
            ls = i.rstrip().split("\t")
            if ls[1] != "contigLen":
                depthDict[ls[0]]["length"] = int(ls[1])
                depthDict[ls[0]]["depth"] = ls[2]

    # reading silva headers
    silvaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    silva = open(silvaFile)
    for i in silva:
        ls = i.rstrip().split("\t")
        if ls[0] != "primaryAccession":

            Domain = (ls[3].split(";")[0])
            Phylum = (ls[3].split(";")[1])
            if Phylum == "":
                Phylum = Domain

            try:
                Class = (ls[3].split(";")[2])
            except IndexError:
                Class = (ls[3].split(";")[0])

            if Domain in ["Bacteria", "Archaea"]:
                Genus = lastItem(ls[3].split(";"))

                if Genus == "Oikopleura":
                    Domain = "Eukaryota"
                if Genus == "Diplosphaera":
                    Domain = "Eukaryota"
                if Genus == "Planococcus":
                    Domain = "Eukaryota"

            elif Domain in ["Eukaryota"]:
                Genus = ls[4].split(" ")[0]

                if lastItem(ls[3].split(";")) == "Chloroplast":
                    Genus = "Chloroplast" + "_" + ls[4].split(" ")[0]

                if lastItem(ls[3].split(";")) == "Mitochondria":
                    Genus = "Mitochondria" + "_" + ls[4].split(" ")[0]

                if Genus == "uncultured":
                    Genus = "uncultured_" + lastItem(ls[3].split(";"))

                if Genus == "Labrys":
                    Domain = "Bacteria"
                if Genus == "Halofilum":
                    Domain = "Bacteria"
                if Genus == "Bacillus":
                    Domain = "Bacteria"
                if Genus == "Lactobacillus":
                    Domain = "Bacteria"
                if Genus == "Pseudomonas":
                    Domain = "Bacteria"
                if Genus == "Arthrobacter":
                    Domain = "Bacteria"
                if Genus == "Paracoccus":
                    Domain = "Bacteria"
                if Genus == "Ensifer":
                    Domain = "Bacteria"
                if Genus == "Arthrobacter":
                    Domain = "Bacteria"
                if Genus == "Aeromonas":
                    Domain = "Bacteria"
                if Genus == "Acinetobacter":
                    Domain = "Bacteria"
                if Genus == "Edwardsiella":
                    Domain = "Bacteria"
                if Genus == "Mesorhizobium":
                    Domain = "Bacteria"
                if Genus == "Kitasatospora":
                    Domain = "Bacteria"
                if Genus == "Clostridium":
                    Domain = "Bacteria"
                if Genus == "Rhodocista":
                    Domain = "Bacteria"
                if Genus == "Actinomyces":
                    Domain = "Bacteria"

            silvaDict[Genus]["Domain"] = Domain
            silvaDict[Genus]["Phylum"] = Phylum
            silvaDict[Genus]["Class"] = Class
            silvaDict[Genus]["Genus"] = Genus

            silvaDict[Class]["Phylum"] = Phylum
            silvaDict[Class]["Domain"] = Domain
            silvaDict[Class]["Class"] = Class

            silvaDict[Phylum]["Phylum"] = Phylum
            silvaDict[Phylum]["Domain"] = Domain

            silvaDict[Domain]["Domain"] = Domain

    out = open(args.out + ".csv", "w")
    out.write("contig" + "," + "contig_length" + "," + "hits_per_contig" + "," + "cov" + "," + "GC-content" + "," + "Average_AAI" + "," + "closest_blast_hits" + "\n")
    for i in file.keys():
        if args.bam != "NA":
            depth = depthDict[i]["depth"]
            length = depthDict[i]["length"]
        elif args.spades:
            depth = lastItem(i.split("_"))
            length = len(file[i])
        else:
            depth = "Unknown"
            length = len(file[i])
        gc = gcDict[i]
        hitsList = blastDict[i]
        try:
            AAI = statistics.mean(aaiDict[i])
        except statistics.StatisticsError:
            AAI = "NA"
        out.write(i + "," + str(length) + "," + str(len(hitsList) / (length / 1000)) + "," + str(depth) + "," + str(gc) + "," + str(AAI) + ",")

        if args.lvl == "NA":
            for j in hitsList:
                try:
                    out.write(j + "; ")
                except TypeError:
                    pass
        else:
            for j in hitsList:
                try:
                    Genus = j.split(" ")[0]
                    if re.findall(r'symbiont', Genus):
                        Genus = "Bacteria"

                    try:
                        species = j.split(" ")[1]
                        if species == "sp.":
                            species = j.split(" ")[2]
                    except IndexError:
                        species = "unclassified"

                    Domain = silvaDict[Genus]["Domain"]
                    Phylum = silvaDict[Genus]["Phylum"]
                    Class = silvaDict[Genus]["Class"]

                    if len(silvaDict[Genus]) == 0:
                        Domain = "unclassifed"
                        Phylum = "unclassifed"
                        Class = "unclassifed"

                    if args.lvl == "Domain":
                        out.write(Domain + "; ")

                    elif args.lvl == "Phylum":
                        out.write(Phylum + "; ")

                    elif args.lvl == "Class":
                        out.write(Class + "; ")

                    elif args.lvl == "Genus":
                        out.write(Genus + "; ")

                    elif args.lvl == "species":
                        out.write(species + "; ")

                    else:
                        break

                except TypeError:
                    pass

        out.write("\n")
    out.close()

    if args.fa:
        summary = open(args.out)
        out = open(args.out + '-contigs.fa', "w")
        out2 = open(args.out + '-unmatched.contigs.fa', "w")

        for i in summary:
            ls = i.rstrip().split(",")

            if ls[1] != "contig_length":
                length = float(ls[1])
                hitsperkb = float(ls[2])
                gc = float(ls[4])
                try:
                    aai = float(ls[5])
                except ValueError:
                    aai = 100

                if ls[3] != "Unknown":
                    cov = float(ls[3])
                else:
                    cov = 0

                if args.Domain != "NA":
                    # doing the math
                    hits = ls[6].split("; ")
                    totalHits = len(hits)
                    matches = 0
                    for j in hits:
                        Genus = j.split(" ")[0]

                        try:
                            species = j.split(" ")[1]
                            if species == "sp.":
                                species = j.split(" ")[2]
                        except IndexError:
                            species = "unclassified"

                        Domain = silvaDict[Genus]["Domain"]
                        Phylum = silvaDict[Genus]["Phylum"]
                        Class = silvaDict[Genus]["Class"]
                        if len(silvaDict[Genus]) > 0:
                            domain = "unclassifed"
                            phylum = "unclassifed"
                            Class = "unclassifed"

                        if args.Domain != "NA":
                            if args.Phylum != "NA":
                                if args.Class != "NA":
                                    if args.Genus != "NA":
                                        if args.species != "NA":
                                            if species == args.species:
                                                matches += 1
                                        else:
                                            if Genus == args.Genus:
                                                matches += 1
                                    else:
                                        if Class == args.Class:
                                            matches += 1
                                else:
                                    if Phylum == args.Phylum:
                                        matches += 1
                            else:
                                if Domain == args.Domain:
                                    matches += 1
                        else:
                            matches += 1

                    perc = (matches / totalHits) * 100

                else:
                    perc = 100

                if perc >= args.perc and gc >= args.gc and gc <= args.GC and length >= args.l and length <= args.L and cov >= args.cov and cov <= args.COV and aai >= args.aai:
                    out.write(">" + ls[0] + "\n")
                    out.write(file[ls[0]] + "\n")
                else:
                    out2.write(">" + ls[0] + "\n")
                    out2.write(file[ls[0]] + "\n")

        out.close()

    if args.hgt:
        summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        summary = open(args.out + ".csv")
        for i in summary:
            ls = i.rstrip().split(",")
            if ls[2] != "hits_per_contig":
                counter = 0
                for j in ls[6].split("; "):
                    if re.findall(r'Bacteria', j):
                        counter += 1

                if counter > -1:
                    if float(ls[2]) < 0.25:
                        summaryDict[ls[0]] = ls

        prots = open("%s-proteins.faa" % args.g)
        prots = fasta2(prots)

        outSeq = open("%s.hgt.fasta" % args.out, "w")
        out = open("%s.hgt.csv" % args.out, "w")
        out.write("contig,orf,blastHit,seq\n")
        blastDict = defaultdict(lambda: defaultdict(list))
        blast = open(blastFile)
        for i in blast:
            ls = i.rstrip().split("\t")
            contig = allButTheLast(ls[0], "_")
            if contig in summaryDict.keys():
                try:
                    name = (ls[12].split("[")[1])
                    name = name[0:len(name) - 1]
                    Genus = (name.split(" ")[0])
                    Domain = (silvaDict[Genus]["Domain"])
                    blastDict[ls[0]][Domain].append(replace(ls[12], [","], ";"))

                    # if Domain == "Bacteria":
                    #     outSeq.write(">" + ls[0] + "\n")
                    #     outSeq.write(prots[ls[0]] + "\n")
                    #     out.write(
                    #         contig + "," + ls[0] + "," + str(replace(ls[12], [","], ";")) + "," + prots[ls[0]] + "\n")

                except IndexError:
                    pass

        for i in blastDict.keys():
            Bacteria = blastDict[i]["Bacteria"]
            Eukaryota = blastDict[i]["Eukaryota"]
            if len(Bacteria) > len(Eukaryota):
                contig = allButTheLast(i, "_")
                outSeq.write(">" + i + "\n")
                outSeq.write(prots[i] + "\n")
                out.write(contig + "," + i + "," + str(Bacteria[0]) + "," + prots[i] + "\n")

        outSeq.close()
        out.close()

    print("Finished!")






