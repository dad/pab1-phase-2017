#! python

import sys, os, math, random, argparse
import util, biofile, na
from Bio import Phylo
from Bio.Phylo import NewickIO
from io import StringIO

def extractSpeciesName(text):
	res = h.split('[')[1].split(']')[0]
	return res

def makeShortSpeciesName(text):
	flds = text.split()
	if len(flds)>1:
		res = "{}.{}".format(flds[0][0], flds[1])
	else:
		res = flds[0]
	if len(flds)>2:
		res = "{}_{}".format(res, '_'.join(flds[2:]))
	return res

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Update FASTA file with Newick tree, keyed by species names")
	# Required arguments
	parser.add_argument(dest="fasta_in_fname", default=None, help="FASTA input filename")
	parser.add_argument(dest="tree_in_fname", default=None, help="tree (Newick format) filename")
	parser.add_argument(dest="mapping_in_fname", default=None, help="species-name mapping filename")
	# Optional arguments
	parser.add_argument("--generate-short-ids", dest="generate_short_ids", action="store_true", default=False, help="generate short species IDs (e.g. S.cerevisiae) for each sequence?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--fasta-out", dest="fasta_out_fname", default=None, help="FASTA output filename")
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
	fname =os.path.expanduser(options.fasta_in_fname)
	if not os.path.isfile(fname):
	 	raise IOError("# Error: file {} does not exist".format(fname))
	with open(fname,'r') as inf:
		# Read a FASTA file?
		(headers, seqs) = biofile.readFASTA(inf)
	# Read tree
	tree_fname = os.path.expanduser(options.tree_in_fname)
	if not os.path.isfile(tree_fname):
	 	raise IOError("# Error: file {} does not exist".format(tree_fname))
	tree_string = ""
	with open(tree_fname,'r') as inf:
		lines = inf.readlines()
		for line in lines:
			if not line.strip()[0] == '#':
				tree_string += line.strip()
	trees = NewickIO.parse(StringIO(tree_string))
	tree = next(trees)
	# Read mapping file
	map_fname =os.path.expanduser(options.mapping_in_fname)
	if not os.path.isfile(map_fname):
	 	raise IOError("# Error: file {} does not exist".format(map_fname))
	with open(map_fname,'r') as inf:
		map_table = util.readTable(inf, header=True)

	# Create mapping 
	mapping_dict = dict(zip(map_table['species'], map_table['updated.species']))

	# Update the FASTA headers
	#new_headers = []
	#new_seqs = []
	seq_dict = {}
	header_dict = {}
	short_species_names = {}
	for (i,h) in enumerate(headers):
		species_name = extractSpeciesName(h)
		short_name = makeShortSpeciesName(species_name)
		try:
			updated_species_name = mapping_dict[species_name]
			if not na.isNA(updated_species_name):
				new_header = "{}[{}]{}".format(h.split('[')[0], updated_species_name, h.split(']')[1])
				#new_headers.append(new_header)
				#new_seqs.append(seqs[i])
				seq_dict[updated_species_name] = seqs[i]
				header_dict[updated_species_name] = new_header
				short_species_names[updated_species_name] = short_name
		except KeyError as ke:
			print(ke)

	# Iterate over tree and write out FASTA in tree-sorted order
	n_written = 0
	sorted_headers = []
	sorted_seqs = []
	for indiv in tree.get_terminals():
		#print(indiv.name)
		#spec = extractSpeciesName(h)
		spec = indiv.name.strip()
		spec.replace('\\','')
		#print("^{}$\n".format(spec))
		if spec in header_dict:
			hdr = header_dict[spec]
			if options.generate_short_ids:
				hdr = "{} {}".format(short_species_names[spec], hdr)
			seq = seq_dict[spec]
			sorted_headers.append(hdr)
			sorted_seqs.append(seq)
		else:
			#print(dir(indiv))
			info_outs.write("# Can't find {}\n".format(spec))

	if not options.fasta_out_fname is None:
		fasta_outs = util.OutStreams(open(os.path.expanduser(options.fasta_out_fname),'w'))
		biofile.writeFASTA(sorted_seqs, fasta_outs, headers=sorted_headers)
		info_outs.write("# Wrote {} entries to {}\n".format(len(sorted_seqs), options.fasta_out_fname))

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

