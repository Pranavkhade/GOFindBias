#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:40:36 2017

@author: pranavk
"""
import Xlibrary as X
import argparse
import os


parser=argparse.ArgumentParser(description="This is X V0.1. Please use X -h for the help.")
parser.add_argument('gafname', metavar='NameOfTheGapFile', type=str, help='Please provide the name of the .gaf file')
parser.add_argument('loginput', metavar='Log10[1 or 0]', type=int, help='Parse 1 if you want the scale to be Log to the base 10', default=0)
parser.add_argument('topstat', metavar='TopStatistics', type=int, help='How many top entires you want to include?', default=10)
args=parser.parse_args()

gafname=args.gafname
log10=args.loginput
top=args.topstat
###############################################################################Defaults
ontologies=["F","C","P"]
terms=['GO term','PMID']
###############################################################################Set-up
fh=open("Shannon's Statistics.txt","w")
fh.close()
fh=open("Shannon's Statistics.txt","a")
if not os.path.exists('graph_output'):
    os.mkdir('graph_output')

###############################################################################Main
def main(ontology,term):
    file=open(gafname,"r")
    D=X.GOTermCounter(file,ontology,term)
    file.close()
    
    x,y,xlabels,ylabels,logtext=X.GetTopNFromDictionary(D,top,log10)
    X.GenerateBarPlot(x,y,xlabels,ylabels,ontology,term,logtext)
    
    fh.write("%r\t%r\t%r\n"%(term,ontology,round(X.GetShannonEquitability(D),4)))
    return
###############################################################################CallMain
for i in terms:
    for j in ontologies:
        main(j,i)

###############################################################################Post Cleaning
fh.close()