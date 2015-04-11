#!/usr/bin/env python

"""
Author: Luis Cunha
Date: January 30th 2015
Description:
        script to evaluate influenza assembly pipeline.
        It checks for segment completeness, fragments, problems with the sequence, and present of high frequency variants
Usage:
        1) Run the script with .info.cov file after first assembly.
        2) curate (using the log file produce by this script as guidance,
           annotate the sequence
           generate report
        3) rerun this script from within the curated directory, to annotate the output log file with sequence errors
           and high frequency variants
"""

import os, sys, argparse
from operator import itemgetter


class Segments:
        """
        reads info.cov file and builds dictionary of each segment
        """
        def __init__(self, infocovfile):
                self.segments={}
                self.segmentList = ([str(i) for i in range(1, 9)])
		self.infocovfile = infocovfile
                for segment in self.segmentList:
                        self.segments[segment]=[]
                # read file and build segment dictionary
                try:
                        f = open(self.infocovfile)
                except IOError:
                        print('error opening .info.cov file')
                else:
                        with f as fileHandle:
                                fileHandle.readline()
                                lines = fileHandle.readlines()
                                for line in lines:
                                        line_ = line.strip().split()
                                        d = {"id":line_[0], "length":line_[1], "GC":line_[2], "fullcount":int(line_[3]), "bpcovered":line_[4], "fractioncovered":line_[5]}
                                        self.segments[line_[0][0]].append(d)
                """
                if more than one segment present, sort by the number of reads
                write back as a line string
                could (should) continue the code based on a dictionary structure
                but the code below was already written to handle line strings
                so the lazy approach involves less code to change from dictionary to string, than modify the code after
                """
                for segment in self.segmentList:
                        l = self.segments[segment]
                        self.segments[segment]=[]
                        sortedSegments = sorted(l, key=itemgetter('fullcount'), reverse=True)
                        for y in sortedSegments:
                                d2s = y["id"]+"\t"+y["length"]+"\t"+y["GC"]+"\t"+str(y["fullcount"])+"\t"+y["bpcovered"]+"\t"+y["fractioncovered"]
                                self.segments[segment].append(d2s)
	
	@property
        def segments(self):
                return self.segments
	
	@property
	def segmentList(self):
		return self.segmentList

class processSegments:
        """
        process the segment and write the output file
        """
        def __init__(self, segments_, outfile_):
                self.outfile=outfile_
                self.segments = segments_.segments
                self.segmentList = segments.segmentList
                self.proteins = {"1": "PB2", "2":"PB1", "3":"PA", "4":"HA", "5":"NP", "6":"NA", "7":"MP", "8":"NS"}

        def process(self):

                try:
                        open(self.outfile, 'w').close()
                        output = open(self.outfile, "a")
                except IOError:
                        print "error: can't create output file"
                        sys.exit(0)
                for x in self.segmentList:
                        # for each segment from 1..8
                        # write the segment number and protein name, e.g. "Segment 1 (PB2)"
                        l = self.segments[x]
                        output_string = "Segment " + x + " (" + self.proteins[x] + ")\n"
                        output.write(output_string)
                        output_string=""
                        # handle possible outcomes of the assembly
                        # 1) this segment is not present:
                        if len(l)<1:
                                output_string = "\t No sequence found\n"
                                output.write(output_string)
                        # 2) there's only only instance of the segment, which could be complete or fragment:
                        elif len(l)==1:
                                s=l[0].split()
                                meta=s[0].split("|")
                                if meta[2]=="complete":
                                        output_string = "\tAssembly complete\n"
                                        output.write(output_string)
                                else:
                                        output_string = "\tAssembly incomplete, with following fragment found:\n\t"
                                        m = "|".join(meta[:3])
                                        output_string += "\t\t"+ m + " :\t"+s[3] + " reads\n"
                                        output.write(output_string)
                        # 3) there are multiple instances of this segment
                        else:
                                # count the total number of reads, and the number of complete vs fragment instances of
                                # this segment
                                completecount = 0
                                fragmentcount = 0
                                reads = 0.0
                                for segment_ in l:
                                        s = segment_.split()
                                        reads += int(s[3])
                                segment_ = l[0]
                                s=segment_.split()
                                s1 = s[0].split("|")
                                output_string="\t\t"+"|".join(s1[:3]) + " :\t" + s[3] + " reads\n"
                                if s1[2] == "complete":
                                        completecount += 1
                                else: fragmentcount += 1
                                # handle each additional instance, be it a fragment or complete:
                                for segment_ in l[1:]:
                                        s = segment_.split()
                                        s1 = s[0].split("|")
                                        if int(s[3])*100.0/reads>10:
                                                output_string += "\t\t"+"|".join(s1[:3]) + " :\t" + s[3] + " reads\n"
                                                if s1[2] == "complete": completecount+=1
                                                else: fragmentcount += 1
                                # there's no complete instances of this segment. All are fragments
                                if completecount == 0:
                                        output.write("\tAssembly incomplete, with following fragments found:\n")
                                        output.write(output_string)
                                # in addition to a complete instance, there could be more complete, or fragments
                                # annotate only those with a proportion larger than 10% of the total reads
                                elif completecount == 1:
                                        if fragmentcount == 0: #none had fraction > 10%
                                                output.write("\tAssembly complete\n")
                                        else:
                                                output.write("\tMore than one segment detected:\n")
                                                output.write(output_string)
                                else:
                                                output.write("\tMore than one segment detected:\n")
                                                output.write(output_string)

                output.close()


class Annotations:
        def __init__(self, annotation_file):
                self.annotationFile = annotation_file
                self.annotations = {}

        def processannotations(self):
                try:
                        f = open(self.annotationFile)
                except IOError:
                        print('error opening .info.cov file')
                else:
                        with f as fileHandler:
                                lines=fileHandler.readlines()
                                #print lines
                                s=lines[0].split()[1][0]
                                notifications=[]
                                for line in lines[1:]:
                                        if line[0]==">":
                                                self.annotations[s]=notifications
                                                s=line.split()[1][0]
                                                notifications=[]
                                        if "ERROR" in line or "WARNING" in line:
                                                notifications.append(line)
                                self.annotations[s]=notifications
	@property
	def annotations(self):
        	return self.annotations


class Variants:
        def __init__(self, variants_file):
                self.variantsFile = variants_file
                self.variants = {}

        def processvariants(self):
                try:
                        f = open(self.variantsFile)
                except IOError:
                        print('error opening .info.cov file')
                else:
                        with f as fileHandler:
                                lines = fileHandler.readlines()
                                for line in lines[1:]:
                                    try:
                                        l=line.split()
                                        meta = l[0].split("|")
                                        if float(l[2]) > 0.5 and meta[2] == "complete":
                                                if meta[0] not in self.variants:
                                                        self.variants[meta[0]]=[]
                                                        str = "|".join(meta[0:3]) + "\t:\t" + l[1]+"\t" + l[2]
                                                self.variants[meta[0]].append(str)
                                    except:
                                            pass
        @property
	def variants(self):
                return self.variants


class WriteFinalOutput:
        def __init__(self, annotations_, variants_, outfile, assemblyfile):
                self.outfile = outfile
                self.assemblyfile = assemblyfile
                self.annotations = annotations_
                self.variants = variants_

        def write(self):
                try:
                        open(self.outfile, 'w').close() #empty the file if it exists
                        output= open(self.outfile, 'a')
                        f=open(self.assemblyfile, "r")
                except IOError:
                        print "error: can't create output file"
                        sys.exit(0)

                else:
                        with f as fileHandler:
                                lines=fileHandler.readlines()
                                segmentNumber=lines[0].split()[1]
                                output.write(lines[0])
                                content=[]
                                for line in lines[1:]:
                                        if line[0:7]=="Segment":
                                                output.write("".join(content))
                                                if len(self.annotations[segmentNumber])>0:
                                                        for x in self.annotations[segmentNumber]:
                                                                output.write("\t"+x)
                                                if segmentNumber in self.variants:
                                                        output.write("\tVariants present at a frequency greater than 50%:\n")
                                                        for x in self.variants[segmentNumber]:
                                                                output.write("\t\t"+x+"\n")
                                                segmentNumber = line.split()[1]  #reset
                                                content=[]
                                                output.write(line)
                                                continue
                                        content.append(line)

                                output.write("".join(content))
                        output.close()


if __name__ == "__main__":
        def usage(code="0"):
		print "error: " + str(code)
                print "\n\tusage: process_assembly --cov <example.info.cov>\n\tor"
		print "\tusage: process_assembly --anno <example.anno> --var <example.variant.calls.txt>\n"               
		sys.exit(0)
        parser = argparse.ArgumentParser(description='Parse the info.cov file and check for segment completeness and errors')
        parser.add_argument('--cov',  help="info.cov file - must be in the same path where script is called")
        parser.add_argument('--anno', help="annotation file")
	parser.add_argument('--var', help="called variants file")
	args = parser.parse_args()
        # first run of the script, after assembly
	if args.cov!=None:
		print args
		covfile = args.cov
        	try:
			if covfile.split("_assembly_")[1] != "cap3_cdhit.info.cov":
                		usage("1")
	     	except: usage("4")
	        basename = covfile.split("_assembly_")[0]
	        outfile = basename + "_assembly.out"
		segments = Segments(covfile)
                process = processSegments(segments, outfile)
                process.process()
        # second run of the script, after annotation step                                  
	else:
		if args.anno== None or args.var == None:
			usage("2")
		annofile = args.anno
		varfile = args.var
		try:
			if varfile.split("_curated")[1] != ".variants.calls.txt": usage()
			if annofile[-5:] != ".anno": usage()
			basename =  varfile.split("_curated")[0]
			outfile = basename + "_assembly.out"	
		except: usage("3")
		else: 
        		if "curated" in os.getcwd():
				annotation_file = args.anno
                		variant_file = args.var
                		assembly_file = "../"+basename + "_assembly.out"
                		processAnnotations = Annotations(annotation_file)
                		processVariants = Variants(variant_file)
                		processAnnotations.processannotations()
                		annotations = processAnnotations.annotations
                		processVariants.processvariants()
                		variants = processVariants.variants
                		writeFinalOutput = WriteFinalOutput(annotations, variants, outfile, assembly_file)
                		writeFinalOutput.write()
			else: 
				print "If the script is being run post-annotation, make sure to be called from within the curated directory"
				sys.exit(0)
            #   os.remove(assembly_file)
