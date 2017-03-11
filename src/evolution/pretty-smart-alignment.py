#! python

import sys, os, math, random, argparse, collections #, biopython
import util, biofile, muscle, translate, orthodbutil

'''
Make an alignment of SMART sequences pretty by:
1) Sorting by phylogeny
2) Making header "X.yyyyy" from [Genus species]
'''

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Prettify a SMART alignment")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input FASTA filename")
	# Optional arguments
	parser.add_argument("--tree", dest="tree_fname", default=None, help="filename of Newick-format species tree")
	parser.add_argument("--one-per-species", dest="one_per_species", action="store_true", default=False, help="one sequence per species?")
	parser.add_argument("--gap-threshold", dest="gap_threshold", type=float, default=None, help="fraction of the total possible alignment gaps")
	parser.add_argument("--insertion-group-threshold", dest="gap_threshold", type=float, default=None, help="one sequence per species?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--fasta-out", dest="fasta_out_fname", default=None, help="output FASTA filename")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	fasta_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = open(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)
	if not options.fasta_out_fname is None:
		outf = open(options.fasta_out_fname,'w')
		fasta_outs.addStream(outf)
	else:
		# By default, write to stdout
		fasta_outs.addStream(sys.stdout)

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	with open(options.in_fname,'r') as inf:
		# Read a FASTA file?
		(headers, seqs) = biofile.readFASTA(inf)
	info_outs.write("# Read {:d} sequences\n".format(len(seqs)))

	def getSpeciesName(hdr):
		# [Cyanidioschyzon merolae strain 10D]
		return hdr.split('[')[-1].split(']')[0].strip()

	def getShortSpeciesName(hdr):
		"""Turn 'Homo sapiens sapiens' into 'H.sapiens'"""
		d = getSpeciesName(hdr).split()
		genus = d[0]
		species = d[1]
		return "{:s}.{:s}".format(genus[0], species)

	# Prune by species?
	if options.one_per_species:
		new_seqs = []
		new_headers = []
		species_dict = {}
		for (xi, h) in enumerate(headers):
			species = getSpeciesName(h)
			species_dict[species] = xi
		new_indices = sorted(species_dict.values())
		new_seqs = [seqs[xi] for xi in new_indices]
		new_headers = [headers[xi] for xi in new_indices]
		seqs = new_seqs
		headers = new_headers
		info_outs.write("# Retained {:d} sequences, enforcing one per species\n ".format(len(headers)))


	# Goal: find sequences responsible for the production of the most gaps
	# Quantify, for each sequence S, the number of sites T at which S is one of only M with a residue
	# Each sequence gets a list of (M,T) pairs.
	# Sequences that have a unique insertion of length I get (1,I) which is terrible.
	# The total number of gaps induced is G = (N-M)*T, out of a total possible number of N*L
	# Enforce a cutoff, eliminating sequences with fraction of gaps G/(N*L) larger than some threshold?
	if not options.gap_threshold is None:
		gap = '-'
		seqs_to_remove = []
		L = len(seqs[0])
		N = len(seqs)
		gap_dist = [None]*L
		# At each site, identify the set of non-gapped sequences
		total_gaps = 0
		for ci in range(L):
			#n_aligned = [x for x[ci] in seqs].count(gap)
			gap_dist[ci] = [xi for xi in range(N) if seqs[xi][ci] != gap]
			total_gaps += N-len(gap_dist[ci])

		# Find the adjacent gaps.
		# Ignore sites with X% 
		gap_chunks = {}
		for xi in range(N):
			gap_chunks[xi] = []
			cur_chunk = 0
			for ci in range(L):
				seqids = gap_dist[ci]
				M = len(seqids)
				frac_gap = (N-M)/N
				if frac_gap>0.9 and xi in seqids:
					cur_chunk += (N-M)
				else:
					gap_chunks[xi].append(cur_chunk)
					cur_chunk = 0
			gap_chunks[xi] = sorted(gap_chunks[xi], reverse=True)

		def biggest_nth_chunk(n):
			return dict([(xi,gap_chunks[xi][n]/total_gaps) for xi in range(N) if len(gap_chunks[xi])>n])

		offenders_to_remove = []
		for i in range(50):
			bc = biggest_nth_chunk(i)
			sl = sorted(bc.items(), key=lambda x: x[1], reverse=True)
			#print([(round(v,4),x) for (x,v) in sl[:10]])
			for (xi, v) in sl:
				if v > options.gap_threshold:
					offenders_to_remove.append(xi)
		seqs_to_keep = sorted(list(set(list(range(N))).difference(set(offenders_to_remove))))
		seqs = [seqs[xi] for xi in seqs_to_keep]
		headers = [headers[xi] for xi in seqs_to_keep]


	# Then sort by length

	# Put nice headers in

	new_headers = []
	for (xi, h) in enumerate(headers):
		short_spec = getShortSpeciesName(h)
		new_headers.append("{:s} {:s}".format(short_spec, h))
	headers = new_headers

	# Degap the remaining alignment
	gap = '-'
	new_seqs = []
	n = len(seqs)
	# Find positions with at least one non-gap character
	nongap = [xi for xi in range(len(seqs[0])) if (''.join([s[xi] for s in seqs])).count(gap)<n]
	#print nongap
	for s in seqs:
		new_s = ''.join([s[xi] for xi in nongap])
		new_seqs.append(new_s)
	seqs = new_seqs

	# Write output
	n_written = len(seqs)
	biofile.writeFASTA(seqs, fasta_outs, headers)

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.fasta_out_fname is None:
		info_outs.write("# Wrote {} entries to {}\n".format(n_written, options.fasta_out_fname))
		outf.close()

