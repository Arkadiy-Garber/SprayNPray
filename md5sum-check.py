#!/usr/bin/env python3
from collections import defaultdict
import re
import statistics
import numpy as np
import os
import random
from Bio import pairwise2
import argparse
import sys
import textwrap


stopAdjacent = {
    "TTA": "TGA",
    "TCA": "TGA",
    "TCG": "TAG",
    "TAC": "TAG",
    "TAT": "TAG",
    "TGG": "TAG",
    "AAG": "TAG",
    "GAG": "TGA",
    "AAA": "TAA",
    "CAG": "TAG",
    "CGA": "TGA",
    "CAA": "TAA",
    "AGA": "TGA",
}


startAdjacent = {
    "ATG": ["ATA", "ATT", "ATC", "ATT", "ACG"],
    "TTG": ["TTA", "TGG", "TTC", "TCG"],
    "GTG": ["GTA", "GGG", "GTC", "GCG"]
}


stopAdjacent2 = {
    "TAG": ["TAC", "TCG", "TAT", "AAG", "GAG", "AAG", "TTG"],
    "TAA": ["TCA", "TAC", "CAA", "GAA", "AAA"],
    "TGA": ["AGA", "CGA", "GGA", "TGC", "TGG", "TCA"]
}


def Parrot(seq):
    ls = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        ls.append(codon)
    return ls


def Stopper(seq):
    potentialStops = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon in stopAdjacent.keys():
            potentialStops.append(codon)

    randint = (random.randint(0, len(potentialStops)))

    potentialStops = []
    ls = []
    count = 0
    trunctation = 0
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon in stopAdjacent.keys():
            potentialStops.append(stopAdjacent[codon])
            if len(potentialStops) == randint:
                trunctation = count
                ls.append(stopAdjacent[codon])
            else:
                ls.append(codon)
                count += 3
        else:
            ls.append(codon)
            count += 3
    newseq = "".join(ls)
    return newseq, (trunctation/len(seq))*100


def Silencer(seq):
    start = seq[0:3]
    therest = seq[3:]
    try:
        newseq = random.choice(startAdjacent[start]) + therest
        start = "+"
    except KeyError:
        newseq = start + therest
        start = "-"

    ls = (Parrot(newseq))
    oldlength = (len(ls)*3)
    count = 0
    newlength = 0
    for i in ls:
        if i in ["ATG", "GTG", "TTG"]:
            newlength = ((len(ls) - count)*3)
            break
        else:
            count += 1

    return newseq, start, (newlength/oldlength)*100


def Unstopper(seq):
    stop = seq[len(seq)-3:len(seq)]
    therest = seq[0:len(seq)-3]
    try:
        return therest + random.choice(stopAdjacent2[stop])
    except KeyError:
        return therest + random.choice(stop)


def Shifter(seq):
    ls = []
    inLenList = [1, 2, 3, 4, 5, 6]
    weightList = [1, 1, 1, 1, 1, 1]

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        ls.append(codon)

    randint = (random.randint(0, 15))
    randints = random.choices(inLenList, weights=weightList, k=randint)
    insertList = []
    for j in randints:
        nucList = (random.choices(["A", "G", "T", "C"], weights=[1, 1, 1, 1], k=j))
        insert = "".join(nucList)
        insertList.append(insert)

    try:
        bottom = (len(ls) / len(insertList))
        top = 1
    except ZeroDivisionError:
        top = 0
        bottom = 1

    ls = []
    inList = []
    delList = []
    for i in range(0, len(seq), 6):
        codon = seq[i:i + 6]
        binary = random.choices(["T", "F"], weights=[top, bottom], k=1)
        if binary[0] == "T":
            binary2 = random.choices(["in", "del"], weights=[1, 2], k=1)
            if binary2[0] == "in":
                randomInsert = random.choice(insertList)
                inList.append(len(randomInsert))
                ls.append(randomInsert)
                ls.append(codon)
            else:
                delLen = random.choice([1, 2, 3, 4, 5, 6])
                delList.append(delLen)
                ls.append(codon[0:len(codon)-delLen])
        else:
            ls.append(codon)

    newseq = "".join(ls)
    ls = []
    for i in range(0, len(newseq), 3):
        codon = newseq[i:i + 3]
        ls.append(codon)

    newseq = "".join(ls)
    alignment = pairwise2.align.globalxx(ribosome(seq), ribosome(newseq))
    prot = (alignment[0][0])
    newprot = (alignment[0][1])
    aai = (AAI(prot, newprot))
    perc = (len(newprot.split("*")[0]) / len(prot)) * 100
    return newseq, aai, perc, sum(inList), sum(delList)


def AAI(seq1, seq2):
    counter = 0
    for i in range(len(seq1)):
        if "-" not in [seq1[i], seq2[i]]:
            if seq1[i] == seq2[i]:
                counter += 1
    return (counter/len(seq1))*100


def RemoveLeadingSpaces(line):
    counter = 0
    newLine = ''
    for i in line:
        if i == " ":
            if counter == 0:
                pass
            else:
                newLine += i
        else:
            counter += 1
            newLine += i
    return newLine


def Unique(ls):
    unqList = []
    for i in ls:
        if i not in unqList:
            unqList.append(i)
    return unqList


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
            prot.append("")
    protein = ("".join(prot))
    return protein


def stabilityCounter6(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = str(0) + (str(0) + str(int))
        return (string)
    if len(str(int)) == 5:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 5:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


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
    return x[0:len(x)-1]


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
    prog="abu.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************

    Developed by Arkadiy Garber; agarber4@asu.edu
    ************************************************************************
    '''))


parser.add_argument('-i', type=str, help="input dir", default="NA")

parser.add_argument('-o', type=str, help="output dir", default="NA")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

DIR = os.listdir(args.o)
bins = os.listdir(args.o + "/bins")

os.system("which md5 > md5.path")
path = open("md5.path")
md5 = 0
for i in path:
    if len(i) > 0:
        md5 = "1"


# if len(DIR) != 9:
#     print(DIR)
#     print("One or more files is missing from the output.")
#     raise SystemExit


# if len(bins) != 4:
#     print(bins)
#     print("One or more files is missing from the output.")
#     raise SystemExit

md5Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
md5 = open("%s/md5sum.txt" % args.i)
for i in md5:
    ls = i.rstrip().split("\t")
    md5Dict[ls[0]] = ls[1]

if md5 == "1":
    for i in DIR:
        if i != "bins" and not re.findall(r'tiff', i):
            os.system("md5 %s/%s >> %s/%s.md5sum" % (args.o, i, args.o, i))
            file = open("%s/%s.md5sum" % (args.o, i))
            original = md5Dict[i]
            test = ""
            for j in file:
                test = lastItem(j.rstrip().split(" "))
            if test != original:
                print(i)
                print(test)
                print(original)
                print("One or more files do not have it's expected contents...")
                raise SystemExit
            os.system("rm %s/%s.md5sum" % (args.o, i))

    for i in bins:
        os.system("md5 %s/bins/%s >> %s/bins/%s.md5sum" % (args.o, i, args.o, i))
        file = open("%s/bins/%s.md5sum" % (args.o, i))
        original = md5Dict[i]
        test = ""
        for j in file:
            test = lastItem(j.rstrip().split(" "))
        if test != original:
            print(i)
            print(test)
            print(original)
            print("One or more files do not have it's expected contents...")
            raise SystemExit
        os.system("rm %s/bins/%s.md5sum" % (args.o, i))

else:
    for i in DIR:
        if i != "bins" and not re.findall(r'tiff', i):
            os.system("md5sum %s/%s >> %s/%s.md5sum" % (args.o, i, args.o, i))
            file = open("%s/%s.md5sum" % (args.o, i))
            original = md5Dict[i]
            test = ""
            for j in file:
                test = j.rstrip().split(" ")[0]
            if test != original:
                print(i)
                print(test)
                print(original)
                print("One or more files do not have it's expected contents...")
                raise SystemExit
            os.system("rm %s/%s.md5sum" % (args.o, i))

    for i in bins:
        os.system("md5sum %s/bins/%s >> %s/bins/%s.md5sum" % (args.o, i, args.o, i))
        file = open("%s/bins/%s.md5sum" % (args.o, i))
        original = md5Dict[i]
        test = ""
        for j in file:
            test = j.rstrip().split(" ")[0]
        if test != original:
            print(i)
            print(test)
            print(original)
            print("One or more files do not have it's expected contents...")
            raise SystemExit
        os.system("rm %s/bins/%s.md5sum" % (args.o, i))

print("All output files look accurate. Seems like SprayNPray and all its dependencies are installed ad functional!")














