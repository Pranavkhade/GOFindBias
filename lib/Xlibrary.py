#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:38:11 2017

@author: pranavk
"""

import re
import operator
import matplotlib.pyplot as plt
import math
import Bio.UniProt.GOA as goa

###############################################################################Functions

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
