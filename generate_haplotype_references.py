#!/usr/bin/env python3

###############################
# Harvey Orel
# Date: 27 Apr 2023
#
# Description: Splits sequences into separate gene folders to be used as references for read mapping.
#
# Note: Input is a directory of fasta alignments. Output is to "Phased_References/" in launch directory;
#
#       Arguments are: <input_alignment_directory>
###############################

import os
import sys
import re
import shutil
from Bio import SeqIO


# Get cwd
directory = str(os.getcwd())

# input directory containing alignment files
input_directory = str(sys.argv[1])
output_directory = directory + "/Phased_References/"

# Create output directory
if not os.path.exists(output_directory):
	os.makedirs(output_directory)
else:
	shutil.rmtree(output_directory)
	os.makedirs(output_directory)
	print("Warning: directory 'Phased_References' already exists... Overwriting")


for alignment in os.listdir(input_directory):
	if alignment.endswith(".fasta") or alignment.endswith(".fas") or alignment.endswith(".FNA"):
		alignment_basename = os.path.basename(alignment)
		gene_name = alignment_basename.split('.')[0]

		os.makedirs(output_directory + gene_name)

		# Parse the input FASTA file and extract sequences
		for record in SeqIO.parse(input_directory + alignment, "fasta"):
			# Create a file name based on the sequence ID
			file_name = record.id + ".fasta"

			# Write the sequence to a new file in the gene directory
			with open(os.path.join(output_directory, gene_name, file_name), "w") as output_file:
			  SeqIO.write(record, output_file, "fasta")

		print("Alignment '%s' sequences written to %s%s" %(alignment_basename, output_directory, gene_name))

print("\n\nDone!!! \n")
