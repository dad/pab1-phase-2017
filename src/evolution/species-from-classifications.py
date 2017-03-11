#! python

import sys, os, math, random, argparse
import util, biofile, phyloutil, na
from io import StringIO
from Bio import Phylo
from Bio.Phylo import Newick
from Bio.Phylo import NewickIO


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Build phylogeny from NCBI classification information")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="guide filename")
	# Optional arguments
	parser.add_argument("--binomial", dest="binomial", action="store_true", default=False, help="output binomial species names (like Homo sapiens)")
	parser.add_argument("--verbose", dest="verbose", action="store_true", default=False, help="write a header and other information?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-g", "--debug", dest="debug", action="store_true", default=False, help="debugging?")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = open(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)

	# Write out parameters
	if options.verbose:
		data_outs.write("# Run started {}\n".format(util.timestamp()))
		data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
		data_outs.write("# Parameters:\n")
		optdict = vars(options)
		for (k,v) in optdict.items():
			data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Read input
	fname =os.path.expanduser(options.in_fname)
	if not os.path.isfile(fname):
	 	raise IOError("# Error: file {} does not exist".format(fname))

	# Get directory of guide file
	path = os.path.dirname(fname)
	curwd = os.getcwd()

	species_names = []
	with open(fname,'r') as inf:
		os.chdir(path)
		tab = util.readTable(inf, header=True)
		rows = tab.dictrows
		if options.debug:
			# Just work on 
			rows = [x for x in tab.dictrows][:2]
		just_started = True
		for row in rows:
			spec_fname = row['filename']
			#print(spec_fname)
			if not na.isNA(spec_fname):
				spec_inf = util.readTable(open(spec_fname,'r'), header=True)
				species = None
				if options.binomial:
					row = phyloutil.getClassificationEntryByRank(spec_inf,'species')
					species = row['name']
				else:
					species = spec_inf.row(spec_inf.nrows-1)['name']
				if not species is None:
					species_names.append(species)

	n_written = 0
	if options.verbose:
		dout = util.DelimitedOutput()
		dout.addHeader('species','Species name','s')
		dout.addHeader('hydrophobicity','Mean hydrophobicity','f')
		if options.randomize:
			dout.addHeader('species','Number of randomizations'.format(options.reps),'d')
			dout.addHeader('hydrophobicity.rand','Mean hydrophobicity after randomization {} times'.format(options.reps),'f')
			dout.addHeader('hydrophobicity.rand.sd','Standard deviation of hydrophobicity after randomization {} times'.format(options.reps),'f')
		dout.describeHeader(data_outs)
	for species in species_names:
		data_outs.write("{}\n".format(species))
		n_written += 1

	# Write out stopping time
	if options.verbose:
		data_outs.write("# Run finished {}\n".format(util.timestamp()))
	# Re-set the working director to initial value
	os.chdir(curwd)

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

