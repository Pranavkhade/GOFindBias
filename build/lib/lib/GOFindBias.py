#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:40:36 2017

@author: pranavk
@version: 1.2.3b1
"""
import argparse
import os
import re
import operator
import matplotlib.pyplot as plt
import math
import Bio.UniProt.GOA as goa
import sys
###############################################################################Library Functions
"""This function sorts the dictionary according to the values"""
def ValueSortDictionary(D):
    sorted_D = sorted(D.items(), key =operator.itemgetter(1))
    return sorted_D

"""This function returns all the necessary elements for the plot of top N elements"""
def GetTopNFromDictionary(D,top,log10):
    sorted_D = ValueSortDictionary(D)
    xlabels=[]
    ylabels=[]
    x=range(1,top+1,1) #Changes according to top parameter
    y=[]

    for i in range(-1,-top-1,-1):
        xlabels.append(sorted_D[i][0])
        ylabels.append(sorted_D[i][1])
        if(log10==1):
            y.append(math.log(sorted_D[i][1],10))
            logtext="Log of "
        else:
            y.append(sorted_D[i][1])
            logtext=""
    return x,y,xlabels,ylabels,logtext

    

"""This function reads the .gaf file and converts it into a dictionary. Please note that root_terms are excluded in this"""
def GOTermCounter(file,ontology,term):
    D={}
    gaf=goa.gafiterator(file)
    for entry in gaf:
        if(entry['GO_ID'] in root_terms or entry['Evidence'] not in EC):
            """print entry['GO_ID']"""
        else:
            if(entry['Aspect']==ontology and term=='GO term'):
                if(entry['GO_ID'] not in D.keys()):
                    D[entry['GO_ID']]=1
                else:
                    D[entry['GO_ID']]+=1
            if(entry['Aspect']==ontology and term=='PMID'):
                for refs in entry['DB:Reference']:
                    if(re.match("PMID",refs)):
                        if(refs not in D.keys()):
                            D[refs]=1
                        else:
                            D[refs]+=1
    return D


"""This function generates a bar plot of top n entries in decreasing order of their occurance"""
def GenerateBarPlot(x,y,xlabels,ylabels,ontology,term,logtext):
    plt.figure(figsize=(20,20))
    plt.bar(x,y)
    if(term=='GO term'):
        plt.ylabel(logtext+'GO term count')
        plt.xlabel('GO IDs')
    elif(term=='PMID'):
        plt.ylabel(logtext+'PMID count')
        plt.xlabel('PMIDs')

    plt.title('Ontology:'+ontology+' & Species:Yeast')
    plt.xticks(x,xlabels,rotation='vertical')

    for i, v in enumerate(x):
        plt.text(v-0.125, y[i]+0.05, ylabels[i], color='black', fontweight='bold')
        
    plt.savefig(outfolder+'/'+term+'_'+ontology+'.png',bbox_inches='tight')
    return

"""This is database search option, you can find anything about GOID"""
def GetOBOData(GOID):
    file=open("go.obo","r")
    allGOterms=file.read().split("[Term]")
    for goterm in allGOterms:
        temp=goterm.split("\n")[1]
        if(GOID in temp):
            return goterm.split("\n")
        
"""This function will fetch you shanon value for entire database"""
def GetShannonIndex(D):
    shannon=0
    sumofD=sum(D.values())
    for key,value in D.items():
        freq=float(value)/sumofD
        shannon=shannon+freq*math.log(freq)
    return -shannon

"""This function is forShannon's Equitability calculation"""
def GetShannonEquitability(D):
    return GetShannonIndex(D)/math.log(len(D))

"""This function is for comparing two files with -cmpr argument"""
def CompareGAFs(files,term,ontology):
    file1=open(files[0])
    file2=open(files[1])
    data1=GOTermCounter(file1,ontology,term)
    data2=GOTermCounter(file2,ontology,term)
    topn1=GetTopNFromDictionary(data1,top,log10)
    topn2=GetTopNFromDictionary(data2,top,log10)
    
    outfile="COMPARE.txt"
    fh=open(outfile,'a')
    fh.write("\nFiles:"+'\t'.join(files)+"\nCommon term type:"+term+"\nCommon Ontology type:"+ontology+"\n"+term+"\t"+files[0]+"\t"+files[1]+"\n")
    for i,n1 in enumerate(topn1[2]):
        for j,n2 in enumerate(topn2[2]):
            if(n1==n2):
                fh.write(str(n1)+'\t'+str(topn1[3][i])+'\t'+str(topn2[3][j])+"\n")
    fh.close()
    return

"""Argument parser"""
def ArgParse():
    parser=argparse.ArgumentParser(description="GOFindBias is an analytical tool built to analyse the .gaf files. Please visit: https://github.com/Pranavkhade/GOFindBias for more details.")
    required = parser.add_mutually_exclusive_group(required=True)
    required.add_argument('-i','--input',nargs='+',metavar='GAF_FILE', type=str, help='Names of the input GAF file(s).')
    required.add_argument('-cmpr','--compare',nargs=2,metavar='FILENAME', type=str,help='Names of the two GAF files')
    
    parser.add_argument('-ls','--logscale',metavar='1 OR O', type=int, help='For graphs, 0: Counts in normal scale 1: Counts in log scale [default=0]', default=0)
    parser.add_argument('-ts','--topstat',metavar='TOP', type=int, help='Top n statistics sorted from highest to lowest [default=10]', default=10)
    parser.add_argument('-e','--evidence',nargs='+',metavar='EVIDENCE_CODE',type=str,help='Accepts Standard Evidence Codes outlined in (http://geneontology.org/page/guide-go-evidence-codes). All 3 letter code for each standard evidence is acceptable. In addition to that EXPEC is accepted which will pull out all annotations which are made experimentally. COMPEC will extract all annotations which have been done computationally. Similarly, AUTHEC and CUREC are also accepted.',default='EXPEC')
    
    args=parser.parse_args()
    return args


"""User selected EvidenceCode array setter"""
def ECSetter(ECARGS):
    flag=True
    if(ECARGS=='EXPEC'):
        EC=EXPEC
        flag=False
    elif(ECARGS=='COMPEC'):
        EC=COMPEC
        flag=False
    elif(ECARGS=='AUTHEC'):
        EC=AUTHEC
        flag=False
    elif(ECARGS=='CUREC'):
        EC=CUREC
        flag=False
    elif(ECARGS=='IEA'):
        EC=IEA
        flag=False
    
    if(flag==True):
        EC=ECARGS
        
    return EC
    
    
###############################################################################Variables and switches
ROOT_BPO='GO:0008150'
ROOT_CCO='GO:0005575'
ROOT_MFO='GO:0003674'
root_terms=[ROOT_BPO,ROOT_CCO,ROOT_MFO]

EXPEC = ["EXP","IDA","IPI","IMP","IGI","IEP"]

COMPEC = ["ISS","ISO","ISA","ISM","IGC","IBA","IBD","IKR","IRD","RCA"]

AUTHEC = ["TAS","NAS"]

CUREC = ["IC","ND"]

IEA = ["IEA"]
###############################################################################Defaults
ontologies=["F","C","P"]
terms=['GO term','PMID']
###############################################################################Set-up
def set_up(gafname):
    outfolder=re.sub('.gaf','',gafname)
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    
    fh=open(outfolder+"/Shannon's_Statistics.txt","w")
    fh.close()
    return

###############################################################################Main
def run(term,ontology,gafname):
    file=open(gafname,"r")
    D=GOTermCounter(file,ontology,term)
    file.close()
    
    global outfolder
    outfolder=re.sub('.gaf','',gafname)
    
    x,y,xlabels,ylabels,logtext=GetTopNFromDictionary(D,top,log10)
    GenerateBarPlot(x,y,xlabels,ylabels,ontology,term,logtext)
    
    fh=open(outfolder+"/Shannon's_Statistics.txt","a")
    fh.write("%r\t%r\t%r\n"%(term,ontology,round(GetShannonEquitability(D),4)))
    fh.close()
    return

###############################################################################CallMain

def main():
    commandLineArg = sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
        
    global log10,top,EC
    
    args=ArgParse()
    
    EC=ECSetter(args.evidence)
        
    gafnames=args.input
    log10=args.logscale
    top=args.topstat
    
    
    
    
    """Running main GOFindBias"""
    if args.input is not None:
        for i in gafnames:
            set_up(i)
        
        for i in terms:
            for j in ontologies:
                for k in gafnames:
                    run(i,j,k)
    
    """Running Optional Terms"""
    if args.compare is not None:
        fh=open("COMPARE.txt","w")
        fh.close()
        for term in terms:
            for ontology in ontologies:
                CompareGAFs(args.compare,term,ontology)
###############################################################################Post Cleaning
if __name__ == "__main__":
    main()
