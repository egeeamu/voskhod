#! /usr/bin/env python
# -*- coding: UTF8 -*-


import os



listfiles = os.popen("""ls ./*_R1_*.f*q""").read().split("\n")[0:-1]
for i in listfiles:
    forwardfile = i
    reversefile = i.replace("_R1_", "_R2_")
    mergefile = i.replace("_R1_", "_concatenated_")

    #print forwardfile
    print "concatenation in progress : " + mergefile
    os.system("rm -rfv " + mergefile)
    os.system("cat " + forwardfile + " > " + mergefile)
    os.system("cat " + reversefile + " >> " + mergefile)
