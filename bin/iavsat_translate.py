#!/usr/bin/env python

import sys
import re


class TranslateSequence:
    '''
	translate the given DNA sequence to protein sequence
    '''
    def __init__(self):
        self.codon_table = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA":"T","ACC":"T","ACG":"T","ACT":"T","AGA":"R","AGC":"S","AGG":"R","AGT":"S","ATA":"I","ATC":"I","ATG":"M","ATT":"I","CAA":"Q","CAC":"H","CAG":"Q","CAT":"H","CCA":"P","CCC":"P","CCG":"P","CCT":"P","CGA":"R","CGC":"R","CGG":"R","CGT":"R","CTA":"L","CTC":"L","CTG":"L","CTT":"L","GAA":"E","GAC":"D","GAG":"E","GAT":"D","GCA":"A","GCC":"A","GCG":"A","GCT":"A","GGA":"G","GGC":"G","GGG":"G","GGT":"G","GTA":"V","GTC":"V","GTG":"V","GTT":"V","TAA":"","TAC":"Y","TAG":"","TAT":"Y","TCA":"S","TCC":"S","TCG":"S","TCT":"S","TGA":"","TGC":"C","TGG":"W","TGT":"C","TTA":"L","TTC":"F","TTG":"L","TTT":"F"}

    def translate(self, sequence_):
        g=sequence_
        splitted_g = [g[0+i:3+i] for i in range(0, len(g), 3)]
        if len(splitted_g[-1]) != 3:
            splitted_g = splitted_g[:-1]
        trans = [self.codon_table[c] for c in splitted_g]
        return "".join(trans)


class Proteins:
    '''
      Modify the sequence header to include the protein names instad of the segment
      input: annotation file (output of the iavsat_annotate.sh script"
      return: dictionary of DNA sequence headers : protein name
      note that different headers for same segment are unique based on the last column assembly pipeline identifier
    '''
    def __init__(self, annotationFile_):
        self.annotationFile = annotationFile_
        self.proteins = {}
        try:
                f = open(self.annotationFile, "r")
        except:
                print "Error: cannot open .anno file\n"
                sys.exit(0)
        else:
	        with f as file:
            		for line_ in file:
                		if "protein_id" in line_:
                    			header = line_.split()[1]
                    			n=next(file)
                    			proteinName = n.strip().split()[1]
                    			if "gene" in n:
                        			self.proteins[header] = proteinName
                    			else:
                        			proteinName = next(file).strip().split()[1]
                        			self.proteins[header] = proteinName

    @property
    def proteins(self):
        return self.proteins


class WriteProteinSequence:
    '''
      Read the coding sequence file, translate the sequence, and write it out to the .protein.fa file
      arguments: coding sequence file name, Headers modified to include protein name instead of segment, output file name
      return: none
      side effect: writes output file
    '''
    def __init__(self, cdsFile_, proteins_, outfile_ ):
        self.cdsFile = cdsFile_
        self.proteins = proteins_
        self.outfile = outfile_

    def write(self):
	try:
		cdsfile = open(self.cdsFile, "r")
	except:
		print "Error: cannot open .cds.fa file\n"
		sys.exit(0)
	else:
        	with cdsfile as cds:
            		seq=[]
           		translateSeq = TranslateSequence()
        		try:
				f=open(self.outfile, "w")
                                f.close()
                		outfile = open(self.outfile, "a")
        		except:
                		print "Error: cannot open output file for writing\n"
				print self.outfile
                		sys.exit(0)
        		else:
                		with outfile as g:
                			lines = cds.readlines()
                			line0 = lines[0].strip().split("|")[1:]
                			l0 = ">"+self.proteins[lines[0].strip()[1:]]+"|"+"|".join(line0)+"\n"
                			g.write(l0)
                			for line in lines[1:]:
                    				if line[0] == ">":
                        				line0 = line.strip().split("|")[1:]
                        				l0 = ">" + self.proteins[line.strip()[1:]] + "|" + "|".join(line0)+"\n"
                        				seq = translateSeq.translate("".join(seq)) + "\n"
                        				seq_formated = re.sub("(.{100})", "\\1\n", seq, 0, re.DOTALL)
                        				g.write(seq_formated)
                        				g.write(l0)
                        				seq = []
                    				else:
                        				seq.append(line.strip())
                			seq = translateSeq.translate("".join(seq)) + "\n"
                			seq_formated = re.sub("(.{100})", "\\1\n", seq, 0, re.DOTALL)
                			g.write(seq_formated)


if __name__ == "__main__":
    if len(sys.argv) != 2:
	print "Usage: iavsat_translate.py <.cds.fa file>\n"
    basename = sys.argv[1].split("_cds")[0]
    cdsFile = sys.argv[1]
    annotationFile = basename + ".anno"
    outfile = basename+"_protein.fa"
    proteins = Proteins(annotationFile)
    writeProteinSequence = WriteProteinSequence(cdsFile, proteins.proteins, outfile)
    writeProteinSequence.write()
