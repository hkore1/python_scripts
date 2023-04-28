#!/usr/bin/env python3

###############################
# Harvey Orel
# Date: 27 Apr 2023
#
# Description: Runs bbsplit to test sample association on a directory of phased references, containing subfolders named by gene 
# with phased haplotype references for each sequence in the dataset inside. (designed for use after 'generate_haplotype_references.py')
#
# Note: Inputs are (1) a directory of phased references, (2) Paired read files for the samples of interest. 
# Output is to "Phased_References/" in launch directory;
#
#       Arguments are: <phased_reference_directory> <read_directory>
###############################

import os
import sys
import re
import shutil
import subprocess
import argparse
import numpy as np

helptext = '''This script runs bbsplit to test sample association on a directory of phased references, containing subfolders named by gene 
 with phased haplotype references for each sequence in the dataset inside.
'''


def createReferencePaths(input_gene_directory):
	'''Given an input <input_gene_directory>, create a list of paths to each reference stored in the directory'''

	reference_paths = []
	for ref in os.listdir(input_gene_directory):
		reference_paths.append(os.path.abspath(input_gene_directory + '/' + ref))
	return(reference_paths)

def createReadfilePaths(input_read_directory):
	'''Given an input <input_read_directory>, create a list of paths to each pair of reads stored in the directory'''

	readlist = []
	for filename in sorted(os.listdir(input_read_directory)):
		readlist.append(filename)

	if len(readlist) % 2 == 0: # checks that readlist is an even number (i.e. all reads are paired)
		ali_array = np.array_split(readlist, len(readlist)/2)
		return(ali_array)

def main():
	parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("phased_reference_directory",help="Directory containing folders for each gene with phased reference sequences from each sample")
	parser.add_argument("read_directory",help="Directory containing paired read files for the samples of interest")
	args = parser.parse_args()

	if len(vars(args)) != 2:
		print("Please supply exactly two arguments!\n")
		sys.exit(1)

	phased_reference_directory = args.phased_reference_directory
	read_directory = args.read_directory

	# Get cwd
	directory = str(os.getcwd())

	output_directory = directory + "/bbsplit_output/"

	# Create output directory
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)
	else:
		shutil.rmtree(output_directory)
		os.makedirs(output_directory)
		print("Warning: directory 'bbsplit_output' already exists... Overwriting")

	# Get read paths
	read_array = createReadfilePaths(read_directory)

	for read_pair in read_array:
		R1 = read_directory + read_pair[0]
		R2 = read_directory + read_pair[1]

		SampleName = read_pair[0].split("_R")[0]

		for folder in os.listdir(phased_reference_directory):
			if os.path.isdir(os.path.join(phased_reference_directory, folder)):
				ref_paths = ','.join(createReferencePaths(os.path.join(phased_reference_directory, folder)))

				p=subprocess.Popen("~/Downloads/bbmap/bbsplit.sh build=%s \
						ref=%s  \
						in=%s \
						in2=%s \
						ambiguous=random \
						ambiguous2=all \
						refstats=%s%s-%s_bbsplit-stats.txt" % (folder,ref_paths,R1,R2,output_directory,SampleName,folder), shell=True).wait()
			else:
				pass



if __name__ == "__main__":
	main()