#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 11:33:03 2015

@author: Luis Cunha

usage: ./get-DIP-breakpoints-star.py prefix

prefix: the basename of the sequencing files. e.g. H112_03_curated.fastq.gz => prefix = H112_03_curated

this program reformats and annotates the 'splice' breakpoints of output of the STAR aligner with the start/end point's coverage (extracted from the plotcov file)

"""

import sys


def usage():
	return "usage:"+"\n"+"\tget-DIP-breakpoints-star.py <SJ.out.tab file> <plotcov file> <output file>"

if len(sys.argv)!=4:
	print usage()
	sys.exit(0)


def parseCoverage(plotcovfile):
    d={}
    with open(plotcovfile, "r") as plotfileHandler:
        plotlines=plotfileHandler.readlines() 
        for line in plotlines:
           l=line.strip().split()
           seg=l[0][0]
           if seg not in d:
               d[seg]={}
           d[seg][int(l[1])]=l[2]
           d[seg]['name']=l[0]
    return d


def annotateDIP(dipfile, outfile, d):
    with open(dipfile, "r") as dipfileHandler:
        out=open(outfile,"a")
        header="Segment"+"\t"+"Start"+"\t"+"End"+"\t"+"Count"+"\t"+"Start_cov"+"\t"+"End_cov"+"\n"
        out.write(header)
        lines=dipfileHandler.readlines()
	for line in lines:
	   l=line.strip().split()
	   seg=l[0][0]
           start=int(l[1])
           end=int(l[2])
           count = int(l[6])
	   covstart=d[seg][start]
           covend = d[seg][end]   
	   name=d[seg]['name']
	   print seg, start, covstart, end, covend
           s=name+"\t"+str(start)+"\t"+str(end)+"\t"+str(count)+"\t"+str(covstart)+"\t"+str(covend)+"\n"   
	   out.write(s)
	out.close()
    return 0

     
if __name__ == "__main__":

    dipfile = sys.argv[1]
    plotcovfile = sys.argv[2]
    outfile = sys.argv[3]

    try:
	coverageDictionary = parseCoverage(plotcovfile)
    except:
	print "Error reading plotcov file\n" 
        print usage()
	sys.exit(0)
    try:
    	annotateDIP(dipfile, outfile, coverageDictionary)
    except:
	print "Error reading dipfile \n" 
	print usage()
	sys.exit(0)

