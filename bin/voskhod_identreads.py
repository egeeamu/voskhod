#! /usr/bin/env python
# -*- coding: UTF8 -*-

from collections import *
from Bio import SeqIO
import os
import time
from Bio.Blast import NCBIXML
from sys import stdout
import sqlite3
import csv
import sys
import datetime
import hashlib
import multiprocessing
import traceback

qq = str(sys.argv[1])
indiv = qq.split("/")[-1]
rundir = "run_" + indiv

step = 10000
cpu = 8 # dont touch this ! For developmental purpose only
hsp_read_percent = 70
wordsize = 9

################################
##DEV PARAMETERS DONT TOUCH !!!!
lenhash = 50 # desactivated dont touch !!
hashseq = True # desactivated dont touch !!
keep_multihits = True
crosspe = "Y" # Y or N
savesequences = "Y"


if crosspe == "Y":
    blastlunch = str("""../../bin/ncbi-blast/bin/blastn -task blastn -num_threads """ + str(cpu) + """ -query ./""" + rundir + """/tmp.fasta -db ./db/db -out ./""" + rundir + """/out.txt -qcov_hsp_perc """ + str(hsp_read_percent) + """ -evalue 0.05 -window_size 0 -outfmt "6 qseqid stitle pident sstrand evalue score qseq nident" -dust yes -soft_masking yes -lcase_masking -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -word_size """ + str(wordsize))
else:
    blastlunch = str("""../../bin/ncbi-blast/bin/blastn -task blastn -num_threads """ + str(cpu) + """ -query ./""" + rundir + """/tmp.fasta -db ./db/db -out ./""" + rundir + """/out.txt -qcov_hsp_perc """ + str(hsp_read_percent) + """ -evalue 0.05 -window_size 0 -outfmt "6 qseqid stitle pident sstrand evalue score qseq nident" -dust yes -soft_masking yes -lcase_masking -word_size """ + str(wordsize))


class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
os.system("mkdir -p ./logs")
logfile = str("./logs/Voskhod_" + indiv + "_" + str(datetime.datetime.now().strftime("%d-%m-%y")) + "_log.txt")
sys.stdout = Logger(logfile)

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
 _    __              __    __                __
| |  / /____   _____ / /__ / /_   ____   ____/ /
| | / // __ \ / ___// //_// __ \ / __ \ / __  /
| |/ // /_/ /(__  )/ ,<  / / / // /_/ // /_/ /
|___/ \____//____//_/|_|/_/ /_/ \____/ \__,_/

¤ Automatic Blast / Reads tagger

Version 20170331
Voskhod Pipeline version V1.1
Part of the Voskhod project
https://github.com/egeeamu/voskhod

(CC-BY-NC-ND 4.0 International license) 
Arnaud Ungaro contact@arnaud-ungaro.fr

"""


#crosspe = "Y"
# indiv = str(raw_input("Enter project name: "))
# step = int(raw_input("Enter stepping (100-10000): "))
# cpu = int(raw_input("Enter core to use (1-n): "))
# wordsize = int(raw_input("Enter word_size (4-30): "))
# crosspe = str(raw_input("Cross-Species ? (Y/N): ")).upper()

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
cpu2 = os.popen('''cat /proc/cpuinfo | grep "model name" | cut -d ":" -f2''').read().split("\n")
cpu2.pop(-1)
blastversion = os.popen("""../../bin/ncbi-blast/bin/blastn -version""").read().split("\n")
blastversion.pop(-1)
blastversion = blastversion[0].split(" ")[1]
#blastpath = os.popen("""whereis blastn""").read().split("\n")
#blastpath.pop(-1)

print "\nIndiv:         " + indiv
print "Blast:           " + str(blastversion)
#print "Blast:           " + str(blastpath).replace("[", "").replace("]", "").replace("'", "").split(" ")[-1]
print "File:            " + str(qq)
print "Dir:             " + str(pwd[0])
print "Server:          " + str(serv[0])
print "Distro:          " + str(distrib[0])
print "Kernel:          " + str(kernel[0])
print "CPU:             " + str(cpu2[0])
print "Date:            " + str(dat[0])
print "User:            " + str(user[0])

print "\nParameters:"
print "Indiv:           " + str(indiv)
print "Stepping:        " + str(step)
print "Core(s):         " + str(cpu)
print "Word-Size:       " + str(wordsize)
print "Cross-Spe:       " + str(crosspe)
print "Hsp/Read percent:" + str(hsp_read_percent)
print "Hash's length:   " + str(lenhash)
print "Blast commands   :" + str(blastlunch)


countreads_start = 0
print "\nChecking input file.."
for record in SeqIO.parse(qq, "fastq"):
    if countreads_start % 100000 == 0 and countreads_start != 0:
        print "."
    countreads_start += 1
print "\n"
if step > int(countreads_start / 5.):
    step = int(countreads_start / 5.)

os.system('rm -rf ./' + rundir)
os.system('rm -rf ./results/table_data_' + indiv + '.db')

os.system('mkdir ./' + rundir)

conn2 = sqlite3.connect('./results/table_data_' + indiv + '.db')
c2 = conn2.cursor()
c2.execute('''CREATE TABLE result (sample varchar, read_name varchar, read_sequence varchar, gene_id varchar, species varchar, read_size int, read_hsp_identity int, identity float, evalue float, gene_name varchar, strand int, hsp varchar, hsp_size int, quality varchar, transcript_id varchar, score float, multi_hits varchar)''')
c2.execute('''PRAGMA synchronous = OFF''')
c2.execute('''PRAGMA journal_mode = OFF''')
c2.execute('''PRAGMA cache_size = 4096''')
conn2.commit()

conn = sqlite3.connect("./data_input/cdna_infos.db")
c = conn.cursor()
c.execute('''PRAGMA synchronous = OFF''')
c.execute('''PRAGMA journal_mode = OFF''')
c.execute('''PRAGMA cache_size = 4096''')
conn.commit()

dico_cdnainfo = OrderedDict()
c.execute("""SELECT DISTINCT transcript_id,gene_id,gene_name,species,transcript_sequence FROM RESULT""")
conn.commit()
match = c.fetchone()

print "Reading cdna_info.db.. "
while match is not None:
    Ensdart = str(match[0])
    Ensdarg = str(match[1])
    Gene_name = str(match[2])
    Specie = str(match[3])
    Seq = str(match[4])
    dico_cdnainfo[Ensdart] = {"Ensdart": Ensdart, "Ensdarg": Ensdarg, "Gene_name": Gene_name, "Specie": Specie}
    match = c.fetchone()


print "Analysis of " + str(countreads_start) + " reads starting.."
dico_reads = OrderedDict()
pctage2 = -1.
hspadd = 0.
readsizeadd = 0.
counthsp = 0.
countloop = 0.
identityadd = 0
identitycount = 0
countreads = 0
countreads_total = 0
count_miss = 0
timestart = round((time.time()), 0)
timestart2 = round((time.time()), 0)
seqsec = 0.
counterror = 0
for record in SeqIO.parse(qq, "fastq"):
    countreads += 1
    countreads_total += 1
    liste_qual = record.letter_annotations["phred_quality"]
    qual_fastaq = ""
    for i in liste_qual:
        qual_fastaq += str(i) + " "
    qual_fastaq = qual_fastaq[:-1]
    sq = str(record.seq)
    readsize = len(sq)
    if hashseq == True:
        seqqhash = sq[0:lenhash]
        hashsqq = hashlib.md5(seqqhash).hexdigest()
    else:
        hashsqq = ""
    Reads_name = str(record.description).replace(" ", "_").replace("\t", "_")
    #Reads_name = str(countreads_total)
    dico_reads[Reads_name] = {"indiv": indiv, "read_name": Reads_name, "read": sq, "ensdarg": "", "specie": "", "read_size": readsize, "read_hsp_identity": "", "identity": 0, "evalue": "", "gene_name": "", "sens": "", "hsp": "", "hsp_size": 0, "qual": qual_fastaq, "ensdart": "", "score": 0.0, "md5_hash": hashsqq, "multihits": []}

    if countreads == step:
        os.system("rm -rf ./" + rundir + "/")
        os.system("mkdir ./" + rundir + "/")
        fasta = str("./" + rundir + "/tmp.fasta")
        fastaw = open(fasta, "a+b")
        for i in dico_reads:
            readname = dico_reads[i]["read_name"]
            read = dico_reads[i]["read"]
            npres = read.upper().count("N")
            ratioofn = float(npres) / len(read)
            if ratioofn <= 0.25:
                fastaw.write(">" + str(readname) + "\n" + read + "\n")
            else:
                print "Too many N in read (" + str(npres) + ") : " + str(readname) + " " + str(read)
        fastaw.close()
        try:
            os.system(blastlunch)
        except Exception:
            print "ERROR WITH BLAST NEAR READ " + str(countreads_total)
            print blastlunch
            counterror += 1
        outputtsv = list(csv.reader(open('./' + rundir + '/out.txt', 'rb'), delimiter='\t'))

        try:
            for i in outputtsv:
                try:
                    query = str(i[0])
                    precscore = dico_reads[query]["score"]
                    score = int(i[5])

                    if score >= precscore and keep_multihits == True:
                        Ensdartmulti = str(i[1])
                        Ensdargmulti = dico_cdnainfo[Ensdartmulti]["Ensdarg"]
                        Speciemulti = dico_cdnainfo[Ensdartmulti]["Specie"]
                        mutihitsstr = Speciemulti + "+" + Ensdartmulti + "+" + Ensdargmulti
                        dico_reads[query]["multihits"].append(mutihitsstr)
                        #print Ensdartmulti,Ensdargmulti,Speciemulti

                    if score > precscore:
                        hsp_blast = str(i[6]).replace("-", "")
                        hsp_size = len(hsp_blast)
                        Evalue = float(i[4])
                        sqsens = str(i[3]).replace("minus", "-1").replace("plus", "1")
                        Ensdart = str(i[1])
                        identity = float(i[2]) / 100.
                        Ensdarg = dico_cdnainfo[Ensdart]["Ensdarg"]
                        Specie = dico_cdnainfo[Ensdart]["Specie"]
                        Gene_name = dico_cdnainfo[Ensdart]["Gene_name"]
                        read_size = dico_reads[query]["read_size"]
                        dico_reads[query]["ensdarg"] = Ensdarg
                        dico_reads[query]["specie"] = Specie
                        dico_reads[query]["read_hsp_identity"] = int(i[7])
                        dico_reads[query]["identity"] = identity
                        dico_reads[query]["evalue"] = Evalue
                        dico_reads[query]["gene_name"] = Gene_name
                        dico_reads[query]["sens"] = sqsens
                        dico_reads[query]["hsp"] = hsp_blast
                        dico_reads[query]["hsp_size"] = hsp_size
                        dico_reads[query]["ensdart"] = Ensdart
                        dico_reads[query]["score"] = score
                        #dico_reads[query]["md5_hash"] = ""
                        identitycount += 1
                        identityadd += identity
                except:
                    print "ERROR PARSING NEAR READ " + str(countreads_total)
                    counterror += 1
        except Exception:
            print "ERROR PARSING NEAR READ " + str(countreads_total)
            counterror += 1
        for i in dico_reads:
            try:
                multihits = "|".join(dico_reads[i]["multihits"])
                read_name = dico_reads[i]["read_name"]
                ensdarg = dico_reads[i]["ensdarg"]
                specie = dico_reads[i]["specie"]
                read_size = dico_reads[i]["read_size"]
                read_hsp_identity = dico_reads[i]["read_hsp_identity"]
                identity = dico_reads[i]["identity"]
                evalue = dico_reads[i]["evalue"]
                gene_name = dico_reads[i]["gene_name"]
                sens = dico_reads[i]["sens"]
                hsp_size = dico_reads[i]["hsp_size"]
                ensdart = dico_reads[i]["ensdart"]
                score = dico_reads[i]["score"]
                md5_hash = dico_reads[i]["md5_hash"]
                #md5_hash = ""
                if savesequences == "Y":
                    read = dico_reads[i]["read"]
                    hsp = dico_reads[i]["hsp"]
                    qual = dico_reads[i]["qual"]
                else:
                    read = ""
                    hsp = ""
                    qual = ""
                if float(hsp_size) > 0:
                    hspadd += hsp_size
                    counthsp += 1
                countloop += 1
                readsizeadd += read_size
                if ensdarg == "":
                    identity = ""
                    hsp_size = ""
                    score = ""
                    read_size = ""
                if ensdarg == "":
                    count_miss +=1
                c2.execute('''INSERT INTO result VALUES ("''' + str(indiv) + '''","''' + str(read_name) + '''","''' + str(read) + '''","''' + str(ensdarg) + '''","''' + str(specie) + '''","''' + str(read_size) + '''","''' + str(read_hsp_identity) + '''","''' + str(identity) + '''","''' + str(evalue) + '''","''' + str(gene_name) + '''","''' + str(sens) + '''","''' + str(hsp) + '''","''' + str(hsp_size) + '''","''' + str(qual) + '''","''' + str(ensdart) + '''","''' + str(score) + '''","''' + str(multihits) + '''")''')
            except:
                counterror += 1
                print "ERROR DICO " + str(i)
        dico_reads.clear()
        conn2.commit()
        identity = round(identityadd/identitycount, 3)
        lenghtnz = len(str(countreads_total))
        nsr = countreads_start - countreads_total
        tmpsr = round(nsr / (float(countreads_total) / (time.time() - timestart)), 0)
        timeremaining = "Remaining time : " + str(datetime.timedelta(seconds=tmpsr))
        seqsec = str(round(float(countreads_total) / (time.time() - timestart), 1))
        seqsecinstant = str(round(float(countreads) / (time.time() - timestart2), 1))
        timestart2 = round((time.time()), 0)
        pctage = round((float(((countreads_total / float(countreads_start)) * 100))), 1)
        prt = "[" + str(int(countreads_total)).center(lenghtnz) + "/" + str(int(countreads_start)).center(lenghtnz) + " - " + '%.1f' % pctage + " %]" + "  Identity: " + '%.1f' % (identity * 100) + " % Identified: " + str( 100 - round(count_miss / float(countreads_total) * 100., 1)) + " % " + timeremaining + "  Seq/s MEAN: " + seqsec + "  Seq/s  : " + seqsecinstant + "  mean_READ : " + str(round(readsizeadd / countloop, 3)) + "  mean_HSP : " + str(round(hspadd / counthsp, 3))
        countreads = 0
        if (countreads_total + step) > countreads_start:
            countreads = (countreads_total + step) - countreads_start
        os.system("rm -rf ./" + rundir + "/")
        if pctage2 != pctage:
            print prt
            pctage2 = pctage
        #stdout.write("\r%s" % str(prt))
        #stdout.flush()
#stdout.write("\n")
print "\n"
conn2.commit()
dat = os.popen("""date""").read().split("\n")
dat.pop(-1)

print "\n[DONE] -  " + str(dat)
print "Speed:      " + str(seqsec)
print "Mean HSP:   " + str(round(hspadd / counthsp, 3))
print "Mean Read:  " + str(round(readsizeadd / countloop, 3))
print "Total time  :" + str(datetime.timedelta(seconds=(round(time.time() - timestart, 0))))
print "Total erros :" + str(counterror)
