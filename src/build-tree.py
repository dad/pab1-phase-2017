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

	tree_root = Newick.Clade()
	tree_root.parent = None
	tree_root.name = "cellular organisms"


	# Get directory of guide file
	path = os.path.dirname(fname)
	curwd = os.getcwd()

	species_names = []
	with open(fname,'r') as inf:
		os.chdir(path)
		tab = util.readTable(inf, header=True)
		rows = tab.dictrows
		if options.debug:
			rows = [x for x in tab.dictrows][:2]
		just_started = True
		for row in rows:
			spec_fname = row['filename']
			#print(spec_fname)
			if not na.isNA(spec_fname):
				spec_inf = util.readTable(open(spec_fname,'r'), header=True)
				twig = phyloutil.treeFromClassificationTable(spec_inf)
				added = phyloutil.mergeTrees(tree_root, twig, add_to_leaf=just_started)
				if added:
					just_started = False
					species_names.append(row['updated.species'])
					#print(spec_fname)
				else:
					info_outs.write("# Didn't add {}\n".format(spec_fname))
				#phyloutil.printTree(tree_root)
	# Testing
	# Write tree
	# Read it back in
	# Extract leaf species
	# Check to make sure they're all the ones we expect
	#stream = StringIO()
	#Phylo.write(tree_root, stream, 'newick')
	#print(stream.getvalue())
	#tmpfile = tempfile.NamedTemporaryFile(mode="w")
	#tmpfile.write(stream.getvalue())
	#print(tmpfile.name)
	#in_tree = NewickIO.parse(stream)
	in_tree = tree_root
	in_species = [n.name for n in in_tree.get_terminals()]
	in_species.sort()
	species_names.sort()
	for (s1, s2) in zip(species_names, in_species):
		#s2_ = s2.replace("\\",'')
		s2_ = s2
		print(s1, s2_)
		assert s1 == s2_

	# Write out tree
	Phylo.write(tree_root, data_outs, 'newick')
	n_written = len(species_names)

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))
	# Re-set the working director to initial value
	os.chdir(curwd)

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

