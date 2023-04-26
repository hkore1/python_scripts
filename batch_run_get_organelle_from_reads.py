#!/Users/bin/python

###########################################
# Harvey Orel
# Date: 15 Oct 2022
# Description: Batch runs GetOrganelle (Jin et al. 2020) on a directory of paired-end read files. 
# Usage: Run from the command line, in the following format:
#
#				python3 batch_run_get_organelle_from_reads.py -go <###> -i <###> -a "<###>"
#
#		Arguments: 	-go = path to directory with get_organelle_from_reads.py
#					-i = path to directory with input read files
#					-a = a string of settings for GetOrganelle assembly, as they would be entered for a single GetOrganelle run
#
# Note: If no assembly settings are provided, it will run with some default settings I've found work for me.
#		It will automatically exclude any non fastq or fastq.gz files, or any files that don't have corresponding
# 		R1 and R2 files, that are in the input directory. Outputs to a 'getorganelle_assemblies' folder in the
#		working directory.
#
# Example commands to run (note quotations around assmebly settings): 
#	python3 batch_run_get_organelle_from_reads.py -go /Users/JohnSmith/opt/miniconda3/bin/ -i read_files/
#	python3 batch_run_get_organelle_from_reads.py -go /Users/JohnSmith/opt/miniconda3/bin/ -i read_files/ -a "-R 5 -k 21,45 -F embplant_pt -t 4 --out-per-round --memory-save"
###########################################

import sys
import os
import argparse
import re

parser = argparse.ArgumentParser(description='Batch run GetOrganelle')

parser.add_argument('-go', metavar='--GetOrganelle_directory', type=str, nargs=1, help='path to directory with get_organelle_from_reads.py')
parser.add_argument('-i', metavar='--input_directory', type=str, nargs=1, help='path to directory with input read files')
parser.add_argument('-a', metavar='--assembly_settings', type=str, nargs='+', help='settings for GetOrganelle assembly, as a string')

args = parser.parse_args()

print("\nListing input arguments...")

### Logic for parsing arguments... ###
#getorganelle directory
if args.go is None:
	print("No direction to get_organelle_from_reads.py, EXITING...!!!")
	sys.exit()
else:
	print("Path to get_organelle_from_reads.py: '" + str(args.go[0]) + "'")
	pass

#input file directory
if args.i is None:
	print("No path to input files, EXITING...!!!")
	sys.exit()
else:
	print("Path to input files: '" + str(args.i[0]) + "'")
	pass

#assembly settings
if args.a is None:
	print("No assembly settings provided, using defaults: '-R 30 -k 21,45,65,85,105 -F embplant_pt -t 1 --out-per-round --memory-save -w 0.67'\n")
	args.a = ["-R","30","-k","21,45,65,85,105","-F","embplant_pt","-t","1","--out-per-round","--memory-save","-w","0.67"]
elif args.a == 'arg_was_not_given':
	print("Option given, but no command-line argument '-a', using defaults: '-R 30 -k 21,45,65,85,105 -F embplant_pt -t 1 --out-per-round --memory-save -w 0.67'\n")
	args.a = ["-R","30","-k","21,45,65,85,105","-F","embplant_pt","-t","1","--out-per-round","--memory-save","-w","0.67"]
else:
	print("Using assembly settings provided: '" + str(' '.join(args.a)) + "'\n")


### Set parameters for script ###
getorganelle_dir = str(args.go[0])
input_dir = str(args.i[0])
assembly_settings = ' '.join(args.a)
output_dir = str(os.getcwd()) + "/getorganelle_assemblies"


### Load get_organelle_from_reads.py ###
sys.path.insert(1, getorganelle_dir)

try:
	import get_organelle_from_reads
	print("Successfully loaded 'get_organelle_from_reads.py'")
except ModuleNotFoundError:
	print("Error importing 'get_organelle_from_reads.py', check file path...")
	sys.exit()


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

print("\nLocated " + str(len(tuples_cleaned)) + " sets of paired read files for assembly:")
for i in tuples_cleaned:
	i = sorted(i)
	print(i)


### Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
else:
	print("Warning: prexisting 'getorganelle_assemblies' folder in working directory, new files will write to that directory")
	pass


### Function to run GetOrg ###
def run_get_org(readfile1, readfile2, output, settings):
	os.system(getorganelle_dir + "get_organelle_from_reads.py" + readfile1 + readfile2 + output + settings)


print("\n")
print("******************************")
print("*** Running GetOrganelle *****")
print("******************************")
print("\n")


for sample in tuples_cleaned:
	pattern = str(re.findall(r".*(?=_R\d{1})", sample[0])[0])
	R1 = str(sample[0])
	R2 = str(sample[1])
	O = output_dir + "/" + pattern + "_output"

	run_get_org(readfile1 = " -1 " + input_dir + R1,
		readfile2 = " -2 " + input_dir + R2,
		output = " -o " + O,
		settings = " " + assembly_settings)

