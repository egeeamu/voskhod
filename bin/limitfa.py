#! /usr/bin/env python
# -*- coding: UTF8 -*-

from Bio import SeqIO
from optparse import OptionParser
from random import shuffle

parser = OptionParser()
parser.add_option("-i", "--input", dest="inp", help ="")
parser.add_option("-o", "--output", dest="out", help ="")
parser.add_option("-n", "--number", dest="num", help ="")
(options, args) = parser.parse_args()
inp = options.inp
out = options.out
num = int(options.num)

fasta = str(out)
fastaw = open(fasta, "a+b")

nbrtotal = 0
count = -1

for record in SeqIO.parse(inp, "fastq"):
    nbrtotal += 1

# print nbrtotal
# print num

if num > nbrtotal:
    num = nbrtotal

# print num

listseq = list(xrange(nbrtotal))
shuffle(listseq)
listprelev = listseq[0:num]
listprelev = sorted(listprelev)

for record in SeqIO.parse(inp, "fastq"):
    count += 1
    if count == listprelev[0]:
        listprelev.pop(0)
        sqname = record.description
        sq = str(record.seq)
        writeoupt = ">" + sqname + "\n" + sq + "\n"
        fastaw.write(writeoupt)
    if len(listprelev) == 0:
        break

fastaw.close()
