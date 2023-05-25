#!/usr/bin/env python3

###########################################
# Harvey Orel
# Date: 15 Oct 2022
# 	Modified: 25 May 2023 - Major update to include running samples in parallel.
# Description: Batch runs GetOrganelle (Jin et al. 2020) on a directory of paired-end read files. 
# Usage: Run from the command line, in the following format:
#
#				python3 batch_run_get_organelle_from_reads.py -go <###> -i <###> -a "<###>" -o <###> -t <###>
#
#		Arguments: 	-go = Path to directory with get_organelle_from_reads.py
#					-i = Path to directory with input read files
#					-a = A string of settings for GetOrganelle assembly, as they would be entered for a single GetOrganelle run
#					-o = The name for the output parent directory (containing subfolders of assemblies)
#					-t = The number of samples to run at a time. The product of this and the -t value applied to GetOrganelle settings should equal
#						 the total number of available CPUs
#
# Note: If no assembly settings are provided, it will run with some default settings I've found work for me.
#		It will automatically exclude any non fastq or fastq.gz files, or any files that don't have corresponding
# 		R1 and R2 files, that are in the input directory. Outputs to the '-o' folder in the
#		working directory.
#
# Example commands to run (note quotations around assmebly settings): 
#	python3 batch_run_get_organelle_from_reads.py -go /Users/JohnSmith/opt/miniconda3/bin/ -i read_files/
#	python3 batch_run_get_organelle_from_reads.py -go /Users/JohnSmith/opt/miniconda3/bin/ -i read_files/ -a "-R 5 -k 21,45 -F embplant_pt -t 4 --out-per-round --memory-save" -o getorg_aseemblies -t 2
###########################################

import sys
import os
import argparse
import re
import concurrent.futures
import subprocess as sp

parser = argparse.ArgumentParser(description='Batch run GetOrganelle')

parser.add_argument('-go', metavar='--GetOrganelle_directory', type=str, nargs=1, help='path to directory with get_organelle_from_reads.py')
parser.add_argument('-i', metavar='--input_directory', type=str, nargs=1, help='path to directory with input read files')
parser.add_argument('-a', metavar='--assembly_settings', type=str, nargs='+', help='settings for GetOrganelle assembly, as a string')
parser.add_argument('-o', metavar='--output_directory', type=str, nargs=1, help='name for output directory, as a string')
parser.add_argument('-t', metavar='--threads', type=int, nargs=1, help='number of threads (i.e. samples to run at a time)')

args = parser.parse_args()

print("\nListing input arguments...\n")




###############################
### Logic for  arguments... ###
###############################

# Test path to getorganelle directory
if args.go is None:
	print("No direction to get_organelle_from_reads.py, EXITING...!!!")
	sys.exit()
else:
	print("Path to get_organelle_from_reads.py: '" + str(args.go[0]) + "'")
	pass

# Test input file directory
if args.i is None:
	print("No path to input files, EXITING...!!!")
	sys.exit()
else:
	print("Path to input files: '" + str(args.i[0]) + "'")
	pass

# Check assembly settings
if args.a is None:
	print("No assembly settings provided, using defaults: '-R 30 -k 21,45,65,85,105 -F embplant_pt -t 1 --out-per-round --memory-save -w 0.67'\n")
	args.a = ["-R","30","-k","21,45,65,85,105","-F","embplant_pt","-t","1","--out-per-round","--memory-save","-w","0.67"]
elif args.a == 'arg_was_not_given':
	print("Option given, but no command-line argument '-a', using defaults: '-R 30 -k 21,45,65,85,105 -F embplant_pt -t 1 --out-per-round --memory-save -w 0.67'\n")
	args.a = ["-R","30","-k","21,45,65,85,105","-F","embplant_pt","-t","1","--out-per-round","--memory-save","-w","0.67"]
else:
	print("Using assembly settings provided: '" + str(' '.join(args.a)) + "'\n")

# Test output directory
if args.o is None:
	print("No output directory name provided, using 'getorganelle_assemblies'")
	args.o = "getorganelle_assemblies"
else:
	print("Name of output directory: '" + str(args.o[0]) + "'")
	pass

# Test threads
if args.t is None:
	print("No thread number specified, using 1")
	args.t = 1
else:
	print("Number of threads: '" + str(args.t[0]) + "'")
	pass


### Confirm parameters for script ###
getorganelle_dir = str(args.go[0])
input_dir = str(args.i[0])
assembly_settings = ' '.join(args.a)
output_dir = str(os.getcwd() + '/' + args.o[0])
samplethreads=int(args.t[0])




########################################
### Load get_organelle_from_reads.py ###
########################################

sys.path.insert(1, getorganelle_dir)

try:
	import get_organelle_from_reads
	print("Successfully loaded 'get_organelle_from_reads.py'")
except ModuleNotFoundError:
	print("Error importing 'get_organelle_from_reads.py', check file path...")
	sys.exit()




##########################################################
###  Organise list of read files to run the program on ###
##########################################################

### Coerce read files in input directory to list of tuples ###
#create list of fastq files
filenames = []
for filename in os.listdir(input_dir):
	if filename.endswith(".fastq") or filename.endswith(".fastq.gz"):
		filenames.append(filename)

#get base file names (filename preceding R1 or R2) and put matches in tuple
tuples = []
for file in filenames:
    pattern = re.findall(r"(.+)(_R\d{1})", file)[0][0]
    matches = tuple([x for x in filenames if x.startswith(pattern)])
    tuples.append(matches)

#remove any items from list that arent a tuple with 2 elements
removed_files = []
tuples_cleaned = []
for i in tuples:
	if len(i) != 2:
		removed_files.append(i)
	else:
		tuples_cleaned.append(i)

#print some feedback (and sort tuples so that R1 is always first)
print("\nRemoved " +  str(len(removed_files)) + " files:")
for i in removed_files:
	print(i)

#sorting
tuples_cleaned_sorted = []
for i in tuples_cleaned:
	i = sorted(i)
	tuples_cleaned_sorted.append(tuple(i))

#finalise tuple list
tuples_cleaned_sorted = list(set(tuples_cleaned_sorted))
print("\nLocated " + str(len(tuples_cleaned_sorted)) + " sets of paired read files for assembly:")
print(tuples_cleaned_sorted)




###############################
### Create output directory ###
###############################

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
elif os.path.exists(output_dir):
	print("Warning: prexisting " + str(args.o[0]) + " folder in the working directory, new files will write to that directory")




########################
### Define functions ###
########################

def glob_re(pattern, strings):
    return filter(re.compile(pattern).match, strings)

def run_get_org(sample):
	# First, get the information for an individual sample...
	pattern = str(re.findall(r".*(?=_R\d{1})", sample[0])[0])
	readfile1 = input_dir + str(sample[0])
	readfile2 = input_dir + str(sample[1])
	output = output_dir + "/" + pattern + "_output"
	sample_name=re.split(r"_H.{8}_", pattern)[0]

	# Run get_organelle_from_reads.py
	p1=sp.Popen("python3 %s/get_organelle_from_reads.py -1 %s -2 %s -o %s %s --prefix %s_" %(getorganelle_dir,readfile1,readfile2,output,assembly_settings,sample_name), shell=True).wait()

	# Get output info from previous command, for input in next command
	getorg_fasta = glob_re(r'1.path_sequence.fasta', os.listdir(output))

	# Run evaluate_assembly_using_mapping.py
	p2=sp.Popen("%s/evaluate_assembly_using_mapping.py -f %s/%s_.fasta -1 %s/extended_1.fq -2 %s/extended_2.fq -o %s/assembly_evaluation -c yes --draw" %(getorganelle_dir,output,getorg_fasta,output,output,output), shell=True).wait()




######################
### Run the script ###
######################

if __name__ == '__main__':

	print("\n")
	print("******************************")
	print("*** Running GetOrganelle *****")
	print("******************************")
	print("\n")

	executor1 = concurrent.futures.ProcessPoolExecutor(samplethreads)
	futures1 = [executor1.submit(run_get_org,sample) for sample in tuples_cleaned_sorted]
	concurrent.futures.wait(futures1)


