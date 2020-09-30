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
parser.add_argument('-out', type=str, help="Name output file", default="NA")
parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST", default=1)
parser.add_argument('--skip_makedb', type=str, help="if the DIAMOND database already exists "
                                                    "(i.e. file with extension .dmnd), and you provided that file to "
                                                    "the -ref argument, then include this flag to skip the time-consuming "
                                                    "step of building the DIAMOND BLAST database", const=True, nargs="?")
parser.add_argument('--spades', type=str, help="is this a SPAdes assembly, with the SPAdes headers? If so, "
                                               "then you can provide this flag, and BinBlaster will summarize using the coverage "
                                               "information provided in the SPAdes headers", const=True, nargs="?")
parser.add_argument('--meta', type=str, help="meta procedure for prodigal", const=True, nargs="?")

args = parser.parse_args()

if args.o != "NA":
    file = open(args.o)
    file = fasta2(file)
    if not args.skip_makedb:
        print("Running Diamond: making DIAMOND BLAST database")
        os.system("diamond makedb --in %s --db %s.dmnd" % (args.ref, args.ref))

        print("Running Diamond: running DIAMOND BLAST")
        os.system("diamond blastp --db %s.dmnd --query %s --outfmt 6 --out %s.blast --max-target-seqs 15 --evalue 1E-6 --threads %d" % (args.ref, args.o, args.o, args.t))

    else:
        print("Running Diamond: skipping makedb and running DIAMOND BLAST")
        os.system(
            "diamond blastp --db %s.dmnd --query %s --outfmt 6 --out %s.blast --max-target-seqs 15 --evalue 1E-6 --threads %d"
            % (args.ref, args.o, args.o, args.t))

    print("extracting DIAMOND BLAST hit information")
    os.system("blast-to-fasta.sh %s.blast %s %s.blast-fasta" % (args.o, args.ref, args.o))

    print("cutting...")
    os.system("cut -f2 %s.blast > ids.txt" % (args.o))
    print("running seqtk...")
    os.system("seqtk subseq %s ids.txt > %s.blast-fasta" % (args.ref, args.o))
    os.system("rm ids.txt")

    print("Preparing summary: %s" % args.out)
    nameDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blastFasta = open("%s.blast-fasta" % args.o)
    for i in blastFasta:
        if re.match(r'>', i):
            line = (i.rstrip()[1:])
            ls = line.split(" ")
            id = (ls[0])
            try:
                name = (allButTheFirst(line[0:150], " "))
                name = name.split("]")[0]
                name = name.split("[")[1]
                nameDict[id] = name
            except IndexError:
                pass
    aaiDict = defaultdict(list)
    blastDict = defaultdict(list)
    blast = open("%s.blast" % args.o)
    for i in blast:
        ls = i.rstrip().split("\t")
        orf = ls[0]
        name = (nameDict[ls[1]])
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

    print("Running Prodigal: calling ORFs from provided contigs")
    if args.meta:
        os.system("prodigal -i %s -a %s-proteins.faa -p meta" % (args.g, args.g))
    else:
        os.system("prodigal -i %s -a %s-proteins.faa" % (args.g, args.g))

    if not args.skip_makedb:
        print("Running Diamond: making DIAMOND BLAST database")
        os.system("diamond makedb --in %s --db %s.dmnd" % (args.ref, args.ref))

        print("Running Diamond: running DIAMOND BLAST")
        os.system("diamond blastp --db %s.dmnd --query %s-proteins.faa --outfmt 6 --out %s.blast --max-target-seqs 1 --evalue 1E-6 --threads %d" % (args.ref, args.g, args.g, args.t))

    else:
        print("Running Diamond: skipping makedb and running DIAMOND BLAST")
        os.system(
            "diamond blastp --db %s.dmnd --query %s-proteins.faa --outfmt 6 --out %s.blast --max-target-seqs 1 --evalue 1E-6 --threads %d"
            % (args.ref, args.g, args.g, args.t))

    print("extracting DIAMOND BLAST hit information")
    os.system("blast-to-fasta.sh %s.blast %s %s.blast-fasta" % (args.g, args.ref, args.g))


    print("cutting...")
    os.system("cut -f2 %s.blast > ids.txt" % (args.g))
    print("running seqtk...")
    os.system("seqtk subseq %s ids.txt > %s.blast-fasta" % (args.ref, args.g))
    os.system("rm ids.txt")

    if args.bam != "NA":
        print("Extracting coverage information from the provided BAM files")
        os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s" % (args.g, args.bam))

    print("Calculating GC-content")
    gcDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    GC = 0
    total = 0
    assembly = open(args.g)
    assembly = fasta2(assembly)
    for i in assembly.keys():
        seq = assembly[i]
        total += len(seq)
        gc = 0
        for bp in seq:
            if bp == "C" or bp == "G":
                GC += 1
                gc += 1
        gcDict[i] = str(float(gc/len(seq)))

    print("Preparing summary: %s" % args.out)
    nameDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blastFasta = open("%s.blast-fasta" % args.g)
    for i in blastFasta:
        if re.match(r'>', i):
            line = (i.rstrip()[1:])
            ls = line.split(" ")
            id = (ls[0])
            try:
                name = (allButTheFirst(line[0:150], " "))
                name = name.split("]")[0]
                name = name.split("[")[1]
                nameDict[id] = name
            except IndexError:
                pass


    aaiDict = defaultdict(list)
    blastDict = defaultdict(list)
    blast = open("%s.blast" % args.g)
    for i in blast:
        ls = i.rstrip().split("\t")
        contig = allButTheLast(ls[0], "_")
        name = (nameDict[ls[1]])
        aai = ls[2]
        blastDict[contig].append(name)
        aaiDict[contig].append(float(aai))

    if args.bam != "NA":
        depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        depth = open("%s.depth" % args.g)
        for i in depth:
            ls = i.rstrip().split("\t")
            depthDict[ls[0]]["length"] = int(ls[1])
            depthDict[ls[0]]["depth"] = ls[2]

    out = open(args.out, "w")
    out.write("contig" + "," + "contig_length" + "," + "hits_per_contig" + "," + "cov" + "," + "GC-content" + "," + "Average_AAI" + "," + "closest_blast_hits" + "\n")
    for i in assembly.keys():
        if args.bam != "NA":
            depth = depthDict[i]["depth"]
            length = depthDict[i]["length"]
        elif args.spades:
            depth = lastItem(i.split("_"))
            length = len(assembly[i])
        else:
            depth = "Unknown"
            length = len(assembly[i])
        gc = gcDict[i]
        hitsList = blastDict[i]
        try:
            AAI = statistics.mean(aaiDict[i])
        except statistics.StatisticsError:
            AAI = "NA"
        out.write(i + "," + str(length) + "," + str(len(hitsList)/length) + "," + str(depth) + "," + str(gc) + "," + str(AAI) + ",")
        for j in hitsList:
            try:
                out.write(j + "; ")
            except TypeError:
                pass
        out.write("\n")

    print("Finished!")







