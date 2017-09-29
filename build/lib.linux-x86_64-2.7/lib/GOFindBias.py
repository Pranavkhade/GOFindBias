#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:40:36 2017

@author: pranavk
"""
import argparse
import os
import re
import operator
import matplotlib.pyplot as plt
import math
import Bio.UniProt.GOA as goa
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
        
    plt.savefig('graph_output/'+term+'_'+ontology+'.png',bbox_inches='tight')
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

"""Do you even lift?"""
def GetShannonEquitability(D):
    return GetShannonIndex(D)/math.log(len(D))
    
###############################################################################Variables and switches
ROOT_BPO='GO:0008150'
ROOT_CCO='GO:0005575'
ROOT_MFO='GO:0003674'
root_terms=[ROOT_BPO,ROOT_CCO,ROOT_MFO]

InferredFromExperiment='EXP'
InferredFromDirectAssay='IDA'
InferredFromPhysicalInteraction='IPI'
InferredFromMutantPhenotype='IMP'
InferredFromGeneticInteraction='IGI'
InferredFromExpressionPattern='IEP'
EC=[InferredFromExperiment,InferredFromDirectAssay,InferredFromPhysicalInteraction,InferredFromMutantPhenotype,InferredFromGeneticInteraction,InferredFromExpressionPattern]


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
def run(ontology,term):
    file=open(gafname,"r")
    D=GOTermCounter(file,ontology,term)
    file.close()
    
    x,y,xlabels,ylabels,logtext=GetTopNFromDictionary(D,top,log10)
    GenerateBarPlot(x,y,xlabels,ylabels,ontology,term,logtext)
    
    fh.write("%r\t%r\t%r\n"%(term,ontology,round(GetShannonEquitability(D),4)))
    return

###############################################################################CallMain

def main():
    parser=argparse.ArgumentParser(description="This is X V0.1. Please use X -h for the help.")
    parser.add_argument('gafname', metavar='NameOfTheGAFFile', type=str, help='Please provide the name of the .gaf file')
    parser.add_argument('loginput', metavar='Log10[1 or 0]', type=int, help='Parse 1 if you want the scale to be Log to the base 10', default=0)
    parser.add_argument('topstat', metavar='TopStatistics', type=int, help='How many top entires you want to include?', default=10)
    args=parser.parse_args()
    
    global gafname,log10,top
    
    gafname=args.gafname
    log10=args.loginput
    top=args.topstat

    for i in terms:
        for j in ontologies:
            run(j,i)
###############################################################################Post Cleaning
if __name__ == "__main__":
    main()
#fh.close()
