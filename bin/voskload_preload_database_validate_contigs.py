#! /usr/bin/env python
# -*- coding: UTF8 -*-

import os
import sqlite3
import sys
import datetime
import math


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


class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
os.system("mkdir -p ./logs")
logfile = str("./logs/Voskload_" + str(datetime.datetime.now().strftime("%d-%m-%y")) + "_log_preload.txt")
sys.stdout = Logger(logfile)

os.system('rm -rf ./db')
os.system('rm -rf ./data_input/db.fasta')
os.system('rm -rf ./results')
os.system('mkdir ./db')
os.system('mkdir ./results')




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
 _    __           __   __                __
| |  / /___  _____/ /__/ /___  ____ _____/ /
| | / / __ \/ ___/ //_/ / __ \/ __ `/ __  /
| |/ / /_/ (__  ) ,< / / /_/ / /_/ / /_/ /
|___/\____/____/_/|_/_/\____/\__,_/\__,_/

¤ Voskload, Voskhod Preloader

Version 20170331
Voskhod Pipeline version V1.1
Part of the Voskhod project
https://github.com/egeeamu/voskhod

(CC-BY-NC-ND 4.0 International license) 
Arnaud Ungaro contact@arnaud-ungaro.fr

"""

conn = sqlite3.connect("./data_input/cdna_infos.db")
c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
conn.commit()

#dico_cdnainfo = OrderedDict()
c.execute("""SELECT DISTINCT transcript_id,gene_id,gene_name,species,transcript_sequence FROM RESULT""")
conn.commit()
match = c.fetchone()
fasta = str("./data_input/db.fasta")
fastaw = open(fasta, "a+b")
listsizes = []
print "Reading cdna_info.db.. "
while match is not None:
    Ensdart = str(match[0])
    Ensdarg = str(match[1])
    Gene_name = str(match[2])
    Specie = str(match[3])
    Seq = str(match[4])
    #dico_cdnainfo[Ensdart] = {"Ensdart": Ensdart, "Ensdarg": Ensdarg, "Gene_name": Gene_name, "Specie": Specie}
    #npres = Seq.upper().count("N")
    npres = 0
    if npres < 4:
        fastaw.write(">" + Ensdart + "\n" + Seq + "\n")
    if npres > 3:
        print "Too many N in Transcript (" + str(npres) + ") : " + str(Ensdart) + " " + str(Seq)
    listsizes.append(len(Seq))
    match = c.fetchone()
fastaw.close()
print "Formating blast database.."

os.system('../../bin/ncbi-blast/bin/makeblastdb -dbtype nucl -in ./data_input/db.fasta -out ./db/db -title Voskhod_DB')
print "Formating blast database ok.. \n"
os.system('rm -rf ./data_input/db.fasta')

print "Mean size:                " + str(int(stat_moyenne(listsizes)))
print "Median size:              " + str(int(stat_mediane(listsizes)))
print "Coefficient of variation: " + str(round(stat_coeffvar(listsizes), 3))
print "Max size:                 " + str(int(max(listsizes)))
print "Min size:                 " + str(int(min(listsizes)))
print "Megabases:                " + str(round(sum(listsizes)/1000000., 3))
