#! /usr/bin/env python
# -*- coding: UTF8 -*-

from Bio import SeqIO
import os
from sys import stdout
import sys
import math
import datetime
import hashlib
from Bio.Seq import Seq
from collections import *


keep_pairwise_sequences = "N"
minimum_sequence_length = 30
maximum_percent_bad_score = 10.
bad_score = 13.
trim_score = 13.
average_sequence_score = 20.
minimum_score_in_sequence = 5.
deduplicate = "N"
len_head_deduplication = 50

qq = str(sys.argv[1])
#qq = str("./test100k.fq")
indiv = qq.split("/")[-1]

class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
os.system("mkdir -p ./logs")
logfile = str("./logs/Vosklean_" + str(indiv) + '_' + str(datetime.datetime.now().strftime("%d-%m-%y")) + "_log_Vosklean.txt")
sys.stdout = Logger(logfile)

# qq = "./BC1_ACAGTG_L001_R1.fastq"


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
 _    __           __   __
| |  / /___  _____/ /__/ /__  ____ _____
| | / / __ \/ ___/ //_/ / _ \/ __ `/ __ \\
| |/ / /_/ (__  ) ,< / /  __/ /_/ / / / /
|___/\____/____/_/|_/_/\___/\__,_/_/ /_/

¤ Vosklean Fastq file cleaner & duplication remover

Version 20170721
Voskhod Pipeline version V1.2
Part of the Voskhod project
https://github.com/egeeamu/voskhod

(GPL-3.0) 
Arnaud Ungaro contact@arnaud-ungaro.fr

"""

serv = os.popen("""hostname""").read().split("\n")
serv.pop(-1)
dat = os.popen("""date""").read().split("\n")
dat.pop(-1)
distrib = os.popen("""cat /etc/issue""").read().split("\n")
distrib.pop(-1)
kernel = os.popen("""uname -a""").read().split("\n")
kernel.pop(-1)
user = os.popen("""whoami""").read().split("\n")
user.pop(-1)
pwd = os.popen("""pwd""").read().split("\n")
pwd.pop(-1)
cpu = os.popen('''cat /proc/cpuinfo | grep "model name" | cut -d ":" -f2''').read().split("\n")
cpu.pop(-1)


print "\n"
print "File name:       " + str(qq)
print "Dir:             " + str(pwd[0])
print "Server:          " + str(serv[0])
print "Distro:          " + str(distrib[0])
print "Kernel:          " + str(kernel[0])
print "CPU:             " + str(cpu[0])
print "Date:            " + str(dat[0])
print "User:            " + str(user[0])
print "\n"


print "Keep_pairwise_sequences:   " + str(keep_pairwise_sequences)
print "Minimum_sequence_length:   " + str(minimum_sequence_length)
print "Maximum_percent_bad_score: " + str(maximum_percent_bad_score)
print "Bad_score:                 " + str(bad_score)
print "Trim_score:                " + str(trim_score)
print "Average_sequence_score:    " + str(average_sequence_score)
print "Minimum_score_in_sequence: " + str(minimum_score_in_sequence)
print "Deduplicate:               " + str(deduplicate)
print "Len_head_deduplication:    " + str(len_head_deduplication)
print "\n"
print "Checking input file.."

def stat_moyenne(ech):
    if len(ech) == 0:
        return 0
    return sum(ech) / float(len(ech))


def stat_variance(ech):
    if len(ech) == 0:
        return 0
    t1 = float(1. / (len(ech) - 1 ))
    moy = stat_moyenne(ech)
    t2 = 0.
    for i in ech:
        add = (i - moy) * (i - moy)
        t2 += add
    return t1 * t2


def stat_ecart_type(ech):
    if len(ech) == 0:
        return 0
    variance = stat_variance(ech)
    return math.sqrt(variance)


def stat_coeffvar(ech):
    if len(ech) == 0:
        return 0
    return stat_ecart_type(ech) / stat_moyenne(ech)


def stat_mediane(ech):
    if len(ech) == 0:
        return 0
    ech.sort()
    if len(ech) / 2 == len(ech) / 2.:
        n1 = float(ech[(len(ech) / 2) - 1])
        n2 = float(ech[(len(ech) / 2)])
        return (n1 + n2) / 2
    else:
        return ech[(len(ech)) / 2]


nbrseq = 0
lenhash = len_head_deduplication
dico_cdnainfo = OrderedDict()
for record in SeqIO.parse(qq, "fastq"):
    if nbrseq % 50000 == 0:
        print "."
    nbrseq += 1
    sequence = str(record.seq)

    headseq = sequence[0:lenhash]
    tailseq = sequence[-lenhash:]

    # headseqrevcomp = str(Seq(headseq).reverse_complement())
    # tailseqrevcomp = str(Seq(tailseq).reverse_complement())
    #
    # headinvert = headseq[::-1]
    # tailinvert = tailseq[::-1]
    # headinvertrevcomp = str(Seq(headinvert).reverse_complement())
    # tailinvertrevcomp = str(Seq(tailinvert).reverse_complement())
    #
    hash_head = hashlib.md5(headseq).hexdigest()
    # hash_tail = hashlib.md5(tailseq).hexdigest()
    #
    # hash_head_revcomp = hashlib.md5(headseqrevcomp).hexdigest()
    # hash_tail_revcomp = hashlib.md5(tailseqrevcomp).hexdigest()
    #
    # hash_headinvert = hashlib.md5(headinvert).hexdigest()
    # hash_tailinvert = hashlib.md5(tailinvert).hexdigest()
    # hash_headseqrevcomp = hashlib.md5(headinvertrevcomp).hexdigest()
    # hash_tailseqrevcomp = hashlib.md5(tailinvertrevcomp).hexdigest()



    try:
        dico_cdnainfo[hash_head]["hits"] += 1
    except:
        dico_cdnainfo[hash_head] = {"hastag": hash_head}
        dico_cdnainfo[hash_head]["hits"] = 1
        dico_cdnainfo[hash_head]["status"] = "raw"

    # try:
    #     dico_cdnainfo[hash_tail]["hits"] += 1
    # except:
    #     dico_cdnainfo[hash_tail] = {"hastag": hash_tail}
    #     dico_cdnainfo[hash_tail]["hits"] = 1
    #     dico_cdnainfo[hash_tail]["status"] = "raw"
    #
    # try:
    #     dico_cdnainfo[hash_head_revcomp]["hits"] += 1
    # except:
    #     dico_cdnainfo[hash_head_revcomp] = {"hastag": hash_head_revcomp}
    #     dico_cdnainfo[hash_head_revcomp]["hits"] = 1
    #     dico_cdnainfo[hash_head_revcomp]["status"] = "raw"
    #
    # try:
    #     dico_cdnainfo[hash_tail_revcomp]["hits"] += 1
    # except:
    #     dico_cdnainfo[hash_tail_revcomp] = {"hastag": hash_tail_revcomp}
    #     dico_cdnainfo[hash_tail_revcomp]["hits"] = 1
    #     dico_cdnainfo[hash_tail_revcomp]["status"] = "raw"
    #
    # try:
    #     dico_cdnainfo[hash_headinvert]["hits"] += 1
    # except:
    #     dico_cdnainfo[hash_headinvert] = {"hastag": hash_headinvert}
    #     dico_cdnainfo[hash_headinvert]["hits"] = 1
    #     dico_cdnainfo[hash_headinvert]["status"] = "raw"
    #
    # try:
    #     dico_cdnainfo[hash_tailinvert]["hits"] += 1
    # except:
    #     dico_cdnainfo[hash_tailinvert] = {"hastag": hash_tailinvert}
    #     dico_cdnainfo[hash_tailinvert]["hits"] = 1
    #     dico_cdnainfo[hash_tailinvert]["status"] = "raw"
    #
    # try:
    #     dico_cdnainfo[hash_headseqrevcomp]["hits"] += 1
    # except:
    #     dico_cdnainfo[hash_headseqrevcomp] = {"hastag": hash_headseqrevcomp}
    #     dico_cdnainfo[hash_headseqrevcomp]["hits"] = 1
    #     dico_cdnainfo[hash_headseqrevcomp]["status"] = "raw"
    #
    # try:
    #     dico_cdnainfo[hash_tailseqrevcomp]["hits"] += 1
    # except:
    #     dico_cdnainfo[hash_tailseqrevcomp] = {"hastag": hash_tailseqrevcomp}
    #     dico_cdnainfo[hash_tailseqrevcomp]["hits"] = 1
    #     dico_cdnainfo[hash_tailseqrevcomp]["status"] = "raw"

print "\nOk.."
print "Cleaning.."
lqualbrute = []
lqualnet = []
lqualbad = []
lensqinit = []
lensqbad = []
lensqgood = []

#resultq = file("./query.fq", "w")

resultq = file(qq + "_query.fq", "w")
megaraw = 0
megaclean = 0
countgood = 0
countbad = 0
count = 0.
countdup = 0
l = range(1, 101)
for record in SeqIO.parse(qq, "fastq"):
    seqq = str(record.seq)[0:lenhash]
    hashsqq = hashlib.md5(seqq).hexdigest()
    valudicohits = dico_cdnainfo[hashsqq]["hits"]
    valudicostatus = dico_cdnainfo[hashsqq]["status"]
    if valudicostatus == "saw":
        countdup += 1
    dico_cdnainfo[hashsqq]["status"] = "saw"

    if deduplicate.upper() == "N":
        valudicostatus = "raw"
        valudicohits = 1

    count += 1
    if int(((count / nbrseq) * 100)) in l:
        l.remove(int(((count / nbrseq) * 100)))
        prt = str(int((count / nbrseq) * 100))
        stdout.write("\r%s" % str(prt) + " %")
        stdout.flush()
    liste_seq = list(record.seq)
    megaraw += len(str(record.seq))
    lensqinit.append(len(record.seq))
    liste_qual = record.letter_annotations["phred_quality"]
    lqualbrute.append(stat_moyenne(liste_qual))

    while int(liste_qual[0]) < trim_score and len(liste_seq) != 1:
        liste_seq.pop(0)
        liste_qual.pop(0)
    while int(liste_qual[-1]) < trim_score and len(liste_seq) != 1:
        liste_seq.pop(-1)
        liste_qual.pop(-1)
    sq = ''.join(liste_seq)
    moyseq = stat_moyenne(liste_qual)
    nbr_elements = len([i for i in liste_qual if i < bad_score])
    nbr_elementsmin = len([i for i in liste_qual if i < minimum_score_in_sequence])
    dict_qual = SeqIO.QualityIO._phred_to_sanger_quality_str
    qual_fastaq = ""
    for i in liste_qual:
        qual_fastaq += dict_qual[i]

    if moyseq > average_sequence_score and (float(nbr_elements) / len(liste_qual)) * 100. < maximum_percent_bad_score and nbr_elementsmin < 1 and len(sq) >= minimum_sequence_length and valudicostatus != "saw":
        megaclean += len(sq)
        lqualnet.append(moyseq)
        lensqgood.append(len(sq))
        resultq.write("@" + str(record.description) + "\n" + str(sq) + "\n" + "+" + "\n" + qual_fastaq + "\n")
        countgood += 1
    else:

        #lqualnet.append(moyseq)
        #lensqgood.append(len(sq))
        if keep_pairwise_sequences.upper() == "Y":
            resultq.write("@" + str(record.description) + "\n" + str("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN") + "\n" + "+" + "\n" + "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" + "\n")
        countbad += 1
stdout.write("\n")

print "\nSummary:"
print "Quality before trim: " + str(round(stat_moyenne(lqualbrute), 3))
print "Quality after trim:  " + str(round(stat_moyenne(lqualnet), 3))
print "Sequences PASSED:    " + str(countgood) + "  " + str(
    round((float(countgood) / (countbad + countgood)) * 100, 2)) + "%"
print "Sequences REJECTED:  " + str(countbad) + "  " + str(
    round((float(countbad) / (countbad + countgood)) * 100, 2)) + "%"

print "\nBefore trim :"
if len(lqualbrute) > 0:
    print "Mean sequence score:      " + str(round(stat_moyenne(lqualbrute), 3))
    print "Median sequence score:    " + str(round(stat_mediane(lqualbrute), 3))
    print "Best sequence score:      " + str(round(max(lqualbrute), 5))
    print "Worst sequence score:     " + str(round(min(lqualbrute), 5))
    print "Coef of var (score):      " + str(round(stat_coeffvar(lqualbrute), 3))
    print "Std dev (score):          " + str(round(stat_ecart_type(lqualbrute), 3))

    print "\nMean sequences lenght:    " + str(round(stat_moyenne(lensqinit), 1))
    print "Median sequences lenght:  " + str(round(stat_mediane(lensqinit), 1))
    print "Max sequence lenght:      " + str(max(lensqinit))
    print "Min sequence lenght:      " + str(min(lensqinit))
    print "Coef of var (lenght):     " + str(round(stat_coeffvar(lensqinit), 3))
    print "Std dev (lenght):         " + str(round(stat_ecart_type(lensqinit), 3))
    print "Megabases                 " + str(round(megaraw / 1000000., 3))

print "\nTrimed :"
if len(lqualnet) > 0:
    print "Mean sequence score:      " + str(round(stat_moyenne(lqualnet), 3))
    print "Median sequence score:    " + str(round(stat_mediane(lqualnet), 3))
    print "Best sequence score:      " + str(round(max(lqualnet), 5))
    print "Worst sequence score:     " + str(round(min(lqualnet), 5))
    print "Coef of var (score):      " + str(round(stat_coeffvar(lqualnet), 3))
    print "Std dev (score):          " + str(round(stat_ecart_type(lqualnet), 3))

    print "\nMean sequences lenght:    " + str(round(stat_moyenne(lensqgood), 1))
    print "Median sequences lenght:  " + str(round(stat_mediane(lensqgood), 1))
    print "Max sequence lenght:      " + str(max(lensqgood))
    print "Min sequence lenght:      " + str(min(lensqgood))
    print "Coef of var (lenght):     " + str(round(stat_coeffvar(lensqgood), 3))
    print "Std dev (lenght):         " + str(round(stat_ecart_type(lensqgood), 3))
    print "Deduplicated sequences:   " + str(countdup)
    print "Megabases:                " + str(round(megaclean / 1000000., 3))
    print "Kept (in bases):          " + str(round((float(megaclean)/megaraw)*100, 2)) + "%"

os.system("rm -rf ./good.fasta")
os.system("rm -rf ./bad.fasta")
os.system("rm -rf ./good.fq")
os.system("rm -rf ./good.fq")
os.system("rm -rf ./bad.fq")
print "\n"

