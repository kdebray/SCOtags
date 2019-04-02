#!/usr/bin/env python
#title           :PhantomSpikesRemover.py
#description     :This script aims to automatically identify phantom spikes in Phylogenetic Informativeness profiles and remove unusual substitution rates causing these phantom spikes
#author          :Kevin Debray
#date            :20190401
#version         :0.1
#usage           :python PhantomSpikesRemover.py -h
#notes           :
#python_version  :3.6.3
#==============================================================================

print("\n\n   ***This is PhantomSpikesRemover.py, a program that automates the removal of unusual subtitution rates that cause phantom spikes in Phylogenetic Informativeness profiles***\n   ***Written by: Kevin Debray (kevin.debray@univ-angers.fr)***\n\n")

import os
import glob
import shutil
import argparse
from itertools import groupby, count

""" Requirements:
	biopython==1.70
        kneed==0.2.4
	numpy==1.15.3
	pandas==0.22.0
	scipy==1.0.0 """
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import Gapped, generic_dna
from Bio.Nexus import Nexus
from kneed import KneeLocator
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--datapoints", dest="datapoints", required=True, type=str, help="Path to reformated datapoints obtained through the PhyDesign application")
parser.add_argument("-r", "--ratefile", dest="ratefile", required=True, type=str, help="Path to rate file obtained through the PhyDesign application")
parser.add_argument("-i", "--infiles", dest="infiles", required=True, type=str, help="Path to the directory of individual alignments to clean")
parser.add_argument("-o", "--outdir", dest="outdir", required=False, type=str, default="./", help="Output file name")
args = parser.parse_args()

DataPoints = args.datapoints
RateFile = args.ratefile
InAln = args.infiles
OutDir = args.outdir

#==============================================================================
#	Section: Functions
#==============================================================================
def SpikesFinder(datapoints):
    """ Scan PI profiles to identify profiles with more than 1 maximu """
    LocList = []
    with open(datapoints, 'r') as infile:
        for l in infile:
            loc = l.strip().split()[0]
            DataPointsList = [float(i) for i in l.strip().split()[1:]] #Convert datapoints to list of floats
            DataPointsArray = np.pad(np.array(DataPointsList), 1, mode = 'reflect') #Convert list to array with the option 'reflect' so that first and last values can further be seen as local maxima
            MaxNb = len(DataPointsArray[argrelextrema(DataPointsArray, np.greater)[0] - 1]) #Get the number of local maxima
            if MaxNb > 1:
                LocList.append(loc)
    return LocList


def BreakPointFinder(ratefile, loc):
    """ Find an elbow in rate distribution that determines a rate threshold """
    with open(ratefile, 'r') as infile:
        for l in infile:
            if l.startswith(loc):
                SortedRateList = [float(i) for i in sorted(set(l.strip().split(':')[-1].split(',')), key=float)]
                BreakPoint = KneeLocator(range(len(SortedRateList)), SortedRateList, curve='convex', direction='increasing')
                RateThreshold = SortedRateList[BreakPoint.knee]               
    return(RateThreshold)


def SitesToRemove(ratefile, loc, thr):
    """ Create a file listing all sites with unusual substitution rate """
    with open(ratefile, 'r') as infile:
        for l in infile:
            if l.startswith(loc):
                RateList = [float(i) for i in l.strip().split(':')[-1].split(',')]
                for pos in range(len(RateList)):
                    if RateList[pos] > thr:
                        with open("SitesToRemove.txt", 'a') as outfile:
                            outfile.write("%s\t%i\n" % (loc, pos))
                
        
def CleanRateLine(ratefile, loc, thr):
    """ Remove sites with unusual substitution rate from line in RateFile.txt """
    with open(ratefile, 'r') as infile:
        for l in infile:
            if l.startswith(loc):
                NewRateLine = loc+':'
                RateList = [float(i) for i in l.strip().split(':')[-1].split(',')]
                for i in range(len(RateList)):
                    if RateList[i] <= thr and RateList[i] != 0:
                        NewRateLine += str(RateList[i])+','
                    elif RateList[i] == 0:
                        NewRateLine += '0,'
    return NewRateLine[:-1]


def CleanAln(aln, ratefile, loc, thr):
    """ Remove sites with unusual substitution rate from input alignment """
    Aln = AlignIO.read(aln, 'nexus')
    with open(ratefile, 'r') as infile:
        for l in infile:
            if l.startswith(loc):
                GoodSites = []
                RateList = [float(i) for i in l.strip().split(':')[-1].split(',')]
                for pos in range(len(RateList)):
                    if RateList[pos] <= thr:
                        GoodSites.append(pos)
                GoodGroups = groupby(GoodSites, key=lambda n, c=count():n-next(c))
                tmp = [list(g) for k, g in GoodGroups]
                GoodIntervals = [str(x[0]) if len(x) == 1 else "{},{}".format(x[0],x[-1]) for x in tmp]
                # Saving alignment
                CleanAln = MultipleSeqAlignment([])
                it = 0
                for interval in GoodIntervals:
                    SubStart = int(interval.split(',')[0])
                    SubEnd = int(interval.split(',')[-1])
                    SubAln = MultipleSeqAlignment(Aln[:, SubStart:SubEnd+1], alphabet=Gapped(generic_dna))
                    if it == 0:
                        NewAln = SubAln
                    else:
                        NewAln += SubAln
                    it = 1
                OutName = OutDir + loc + '.clean.nex'
                with open(OutName, 'w') as outhandle:
                    outhandle.write(NewAln.format('nexus'))


#==============================================================================
#	Section: Main Code
#==============================================================================

# Clean directory
print("Preparing analysis...")
if os.path.exists("SitesToRemove.txt"):
    os.remove("SitesToRemove.txt")

if not os.path.exists(OutDir):
    os.mkdir(OutDir)
else:
   shutil.rmtree(OutDir)
   os.mkdir(OutDir)

# Identify PI profiles with spikes
print("Identifying spikes in PI profiles...")
LocusWithSpikes = SpikesFinder(DataPoints)

# Initiate a clean rate file
NewRateFile = RateFile.split('/')[-1].rsplit('.',1)[0]+'.clean.txt'
if os.path.exists(NewRateFile):
    os.remove(NewRateFile)

# Loop through all nexus input alignments
print("Cleaning...")
for aln in glob.glob(InAln+'*.nex'):
    Locus = aln.split('/')[-1].split('.')[0]
    if Locus in LocusWithSpikes:
        print("...%s: a clean-up is required" % Locus)

        threshold = BreakPointFinder(RateFile, Locus)

        SitesToRemove(RateFile, Locus, threshold)

        NewRateLine = CleanRateLine(RateFile, Locus, threshold)
        with open(NewRateFile, 'a') as outfile:
            outfile.write(NewRateLine+'\n')

        CleanAln(aln, RateFile, Locus, threshold)

    else:
        print("...%s: no clean-up is required" % Locus)
        with open(RateFile, 'r') as infile:
            for l in infile:
                if l.startswith(Locus):
                    with open(NewRateFile, 'a') as outfile:
                        outfile.write(l)
        shutil.copy2(aln, OutDir)
