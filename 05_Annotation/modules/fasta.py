'''
Created on Oct 27, 2015

@author: sven.twardziok

Version 0.9

'''

from Bio import SeqIO
import subprocess, re, itertools, csv

class SplitSeqs(object):
	def __init__(self, sequences, outdir, nfiles=500):
		nseqs = 0
		with open(sequences, "r") as infile:
			for record in SeqIO.parse(infile, "fasta"):
				nseqs += 1
		seqsperfile = nseqs/nfiles
		for i in range(1, nfiles+1):
			tmpoutdir = "%s/part_%i" %(outdir, i)
			subprocess.call(["mkdir", "-p", tmpoutdir])
		self.fasta_parts = {}
		with open(sequences, "r") as infile:
			tmpcounter = 0
			nfile = 1
			for record in SeqIO.parse(infile, "fasta"):
				if nfile < nfiles:
					if tmpcounter==0:
						tmpfilename = "%s/part_%i/part_%i.fasta" % (outdir, nfile, nfile)
						outfasta = open(tmpfilename, "w")
						self.fasta_parts["part_%i" % (nfile)] = tmpfilename
						SeqIO.write(record, outfasta, "fasta")
						tmpcounter = tmpcounter+1
					else:
						SeqIO.write(record, outfasta, "fasta")
						tmpcounter = tmpcounter+1
					if tmpcounter>=(seqsperfile-1):
						outfasta.close()
						tmpcounter = 0
						nfile = nfile+1
				else:
					if tmpcounter==0:
						tmpfilename = "%s/part_%i/part_%i.fasta" % (outdir, nfile, nfile)
						outfasta = open(tmpfilename, "w")
						self.fasta_parts["part_%i" % (nfile)] = tmpfilename
						SeqIO.write(record, outfasta, "fasta")
						tmpcounter = tmpcounter+1
					else:
						SeqIO.write(record, outfasta, "fasta")
						tmpcounter = tmpcounter+1

class PrintCdsStats(object):
	def __init__(self, infasta, outstats):
		stop_codons = ["TGA", "TAG", "TAA"]
		start_codons = ["ATG"]
		with open(infasta, "r") as infile:
			with open(outstats, "w") as outfile:
				rowpattern = {"id":"none", "length": 0, "status": "fragment"}
				variables = ["id", "length", "status"]
				writer = csv.DictWriter(outfile, fieldnames=variables)
				writer.writeheader()
				for record in SeqIO.parse(infile, "fasta"):
					x = str(record.seq)
					outdata = dict(rowpattern)
					outdata["id"] = record.id
					outdata["length"] = len(x)
					if len(str(record.seq)) % 3 != 0:
						outdata["status"] = "no translation"
					elif any(x[i:i+3] in stop_codons for i in range(3,len(x)-3,3)):
						outdata["status"] = "internal stop"
					elif x[0:3] in start_codons and x[(len(x)-3):len(x)] in stop_codons:
						outdata["status"] = "complete"
					elif not x[0:3] in start_codons and x[(len(x)-3):len(x)] in stop_codons:
						outdata["status"] = "no start"
					elif x[0:3] in start_codons and not x[(len(x)-3):len(x)] in stop_codons:
						outdata["status"] = "no stop"
					writer.writerow(outdata)
