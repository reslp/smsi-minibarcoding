#!/usr/bin/env python3
import sys
#from Bio.Blast import NCBIWWW
from Bio import SeqIO
#from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import argparse

if sys.version_info[0] < 3:
	raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="reorient_sequences.py", description = """This script will reorder sequences based on the blast hit direction.""", epilog = """written by Philipp Resl""")
pars.add_argument('--blast_results', dest="blastres", required=True, help="NCBI BLAST results file in XML format.")
pars.add_argument('--sequence_file', dest="seqfile", required=True, help="Sequence file to be checked in FASTA format.")
args=pars.parse_args()

#file_name = "all_barcodes_max750.fa"

#blastx_cline = NcbiblastnCommandline(query=file_name, db="carabidae_blastdb", outfmt=5, out="carabidae_results.xml")
#stdout, stderr = blastx_cline()

sequences = SeqIO.parse(args.seqfile, "fasta")


def get_sequence(seqid):
	for sequence in sequences:
		if seqid == sequence.id:
			return sequence

i = 1
orientation = ""
with open(args.blastres, "r") as f:
	for blast_record in NCBIXML.parse(f):
		if blast_record.alignments:
			if blast_record.alignments[0].hsps:
				strand = blast_record.alignments[0].hsps[0].strand
				if strand[1] == "Minus":
					seq = get_sequence(blast_record.query)
					seq.seq = seq.seq.reverse_complement()
					print(seq.format("fasta").strip())
				else:
					seq = get_sequence(blast_record.query)
					print(seq.format("fasta").strip())	
		else:
			print("Something went wrong with blast for:", blast_record.query, "sequence will be written as-is", file=sys.stderr)
			seq = get_sequence(blast_record.query)
			print(seq.format("fasta"))
