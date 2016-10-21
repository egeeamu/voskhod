#! /usr/bin/env python
# -*- coding: UTF8 -*-

from Bio import SeqIO
import os
from sys import stdout
import sys
import math
import datetime
from collections import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--forward", dest="forward", help="")
parser.add_option("-r", "--reverse", dest="reverse", help="")
(options, args) = parser.parse_args()
forwardfile = str(options.forward)
reversefile = str(options.reverse)

print """
       .
      /|\               .
     / | \            ./|\,
  ,-' \|/ `-.        <-=O=->
<'--==<O>==--`>       '\|/`
  `-. /|\ ,-'           '
     \ | /
      \|/
       '

Version 20160911
Sync R1 & R2  keep only passed pw
Part of the Voskhod project

(C) Arnaud Ungaro
contact@arnaud-ungaro.fr"
"""


LOG_EVERY_N = 100000
count1 = 0
count2 = 0
countbad1 = 0
countbad2 = 0
dico_seq = OrderedDict()
print "Star parsing R1"
for record in SeqIO.parse(forwardfile, "fastq"):
    count1 += 1
    if (count1 % LOG_EVERY_N) == 0:
        print count1
    seq1 = str(record.seq)
    if seq1 == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
        countbad1 += 1
        dico_seq[str(count1)] = {"state": "bad"}
    else :
        dico_seq[str(count1)] = {"state": "ok"}

print "Done"
print "Star parsing R2"
for record in SeqIO.parse(reversefile, "fastq"):
    count2 += 1
    if (count2 % LOG_EVERY_N) == 0:
        print count2
    seq2 = str(record.seq)
    if seq2 == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
        countbad2 += 1
        dico_seq[str(count2)] = {"state": "bad"}
    # else :
    #     dico_seq[str(count2)] = {"state": "ok"}
# print listbat1
# print listbat2
print "Done"
counttoto = 0
for i in dico_seq.keys():
    toto = dico_seq[i]["state"]
    if toto == "bad":
        counttoto += 1
print "Bad in R1: " + str(countbad1)
print "Bad in R2: " + str(countbad2)
print "Common BAD :" + str(counttoto)
count1 = 0
count2 = 0
dict_qual = SeqIO.QualityIO._phred_to_sanger_quality_str
print "Start writing R1"
for record in SeqIO.parse(forwardfile, "fastq"):
    count1 += 1
    if (count1 % LOG_EVERY_N) == 0:
        print count1
    triger = dico_seq[str(count1)]["state"]
    if triger == "ok":
        seq1 = str(record.seq)
        seq1name = str(record.description)
        liste_qual1 = record.letter_annotations["phred_quality"]
        qual_fastaq1 = ""
        for i in liste_qual1:
            qual_fastaq1 += dict_qual[i]
        name1 = "@" + seq1name + "\n" + str(seq1) + "\n" + "+" + "\n" + qual_fastaq1 + "\n"
        #print name
        with open("cleantrinityR1.fq", "a") as myfile:
            myfile.write(name1)

        #print qual_fastaq1
print "Done"
print "Start writing R2"
count1 = 0
for record in SeqIO.parse(reversefile, "fastq"):
    count1 += 1
    if (count1 % LOG_EVERY_N) == 0:
        print count1
    triger = dico_seq[str(count1)]["state"]
    if triger == "ok":
        seq1 = str(record.seq)
        seq1name = str(record.description)
        liste_qual1 = record.letter_annotations["phred_quality"]
        qual_fastaq1 = ""
        for i in liste_qual1:
            qual_fastaq1 += dict_qual[i]

        name1 = "@" + seq1name + "\n" + str(seq1) + "\n" + "+" + "\n" + qual_fastaq1 + "\n"
        #print name
        with open("cleantrinityR2.fq", "a") as myfile:
            myfile.write(name1)

        #print qual_fastaq1
print "Done"
