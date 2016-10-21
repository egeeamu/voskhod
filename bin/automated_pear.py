#! /usr/bin/env python
# -*- coding: UTF8 -*-


import os



listfiles = os.popen("""ls ./*_R1_*.f*q""").read().split("\n")[0:-1]
for i in listfiles:
    forwardfile = i
    reversefile = i.replace("_R1_", "_R2_")
    mergefile = i.replace("_R1_", "_PEAR_Merged_")

    #print forwardfile
    #print reversefile
    "-f forward -r revers -o forawrd_merged -j 8 -n 30 -v 8 | tee forward"
    os.system("./pear-0.9.10-bin-64 -j 8 -n 30 -v 8 -f " + forwardfile + " -r " + reversefile + " -o " + mergefile + " | tee " + mergefile + "_logs.txt")
