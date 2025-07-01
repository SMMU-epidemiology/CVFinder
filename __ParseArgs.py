#!/usr/bin/env python

import argparse
import os

###########################################################
#
# Auxilary functions
#
###########################################################

# Valid file path check (Does not check file formatting, but checks if given
# path exists and is readable)
def valid_infile(in_file):
    if not os.path.isfile(in_file):
        raise argparse.ArgumentTypeError("{0} is not a valid input file path".format(in_file))

    if os.access(in_file, os.R_OK):
        return in_file
    else:
        raise argparse.ArgumentTypeError("{0} is not a readable input file".format(in_file))

# Checking valid MAPQ values (for all values that must be >= 0)
def valid_int(x):
    x = int(x)

    if x < 0:
        raise argparse.ArgumentTypeError("%s must be a non-zero integer" % x)

    return x

# Checking valid integer values (for all values that must be >0)
def positive_int(x):
    x = int(x)

    if x <= 0:
        raise argparse.ArgumentTypeError("%s must be a positive integer" % x)

    return x


###########################################################
#
# main function
#
###########################################################

def parseArgs():
	# 
	parser = argparse.ArgumentParser(
	    description = "="*25 + "  This is the python version of CVFinder.  " + "="*25
	)

	required_inputs = parser.add_argument_group("Mandatory parameter")

	# Required input file from user
	required_inputs.add_argument(
		"-f",
		"--infile",
		type = valid_infile,
		help = "The input BAM/SAM file containing chimeric reads. Output by BWA-MEM."
	)

	# Other optional run parameters
	parser.add_argument(
		"-o",
		"--out",
		default = "Chimeric_Cluster.txt",
		type = str,
		help = "The name of output file. Default: Chimeric_Cluster.txt"
	)
	parser.add_argument(
		"-q",
		"--mapq",
		default = 30,
		type = valid_int,
		help = "Threshold of MAPQ. Default: 30."
	)
	parser.add_argument(
		"-l",
		"--linear",
		action = "store_true",
		help = "If the virus genome is linear or not. Default: circular."
	)
	parser.add_argument(
		"-m",
		"--method",
		default = 0,
		type = int,
		choices = [0, 1],
		help = "Method chosen to estimate the complexity of \
		        sequence. 0 for DUST, 1 for entropy. Default: 0"
	)
	parser.add_argument(
		"-v",
		"--virus",
		default = "chrVirus",
		type = str,
		help = "The name of virus used by aligner. Note that \
		        it must be identical with virus name listed \
		        in the header of input SAM file. Default: chrVirus"
	)
	parser.add_argument(
		"-n",
		"--vlen",
		default = 3215,
		type = positive_int,
		help = "Length of virus genome, default 3215 (which \
		        is the typical length of HBV). This values \
		        will be overwriten by virus length listed \
		        in SAM header."
	)
	parser.add_argument(
		"-s",
		"--cluster_size",
		default = 1,
		type = positive_int,
		help = "Threshold of cluster size, default 1, i.e., \
		        outputing all clusters."
	)
	parser.add_argument(
		"-u",
		"--duplicates",
		action = "store_false",
		help = "Filtering duplicates or not.\
		        Fragments with exactly identical coordinates \
		        (starts, stops, and gaps, etc) were considered \
		        duplicates. Default: filtering"
	)
	parser.add_argument(
		"-d",
		"--distance",
		default = 200,
		type = positive_int,
		help = "Threshold of distance between the fragment and \
		        the current cluster. The distance determines \
		        if a fragment should be added to an existing \
		        cluster or not. Default: 200"
	)
	parser.add_argument(
		"-r",
		"--refine",
		default = 1,
		type = int,
		choices = [0, 1],
		help = "Refining the cluster or not. By default, after making \
	 	        the clusters, those with divergent \
	 	        viral coordinates will be split. Default: 1"
	)
	parser.add_argument(
		"-c",
		"--virus_thre",
		default = 800,
		type = positive_int,
		help = "Threshold during the process of refinement. Default: 800"
	)
	parser.add_argument(
		"-b",
		"--bound_human",
		default = 1000,
		type = positive_int,
		help = "Threshold of human part used in the process of \
		        merging two clusters. Default: 1000"
	)
	parser.add_argument(
		"-w",
		"--bound_virus",
		default = 250,
		type = positive_int,
		help = "Threshold of virus part used in the process of \
	            merging two clusters. Default: 250"
	)
	parser.add_argument(
		"-p",
		"--bp_thre",
		default = 10,
		type = positive_int,
		help = "Threshold of distance when merging breakpoints. Default: 10"
	)


	args = parser.parse_args()
	return args


# Test
if __name__ == "__main__":
	args = parseArgs()

	print(args.infile)
	print(args.out)
	print(args.mapq)

	if args.linear:
		print("The virus is linear.")

	print(args.method)
	print(args.virus)
	print(args.vlen)
	print(args.cluster_size)

	if args.duplicates:
		print("Not filtering the duplicates.")

	print(args.distance)

	if args.refine == 1:
		print("Refining clusters.")

	print(args.virus_thre)
	print(args.bound_human)
	print(args.bound_virus)
	print(args.bp_thre)
