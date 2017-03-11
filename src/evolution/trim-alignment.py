#! python

import sys, os, math, random, argparse, itertools
import util, biofile, muscle, translate, geneutil
import scipy as sp

'''
From http://web.expasy.org/decrease_redundancy/
Outline of algorithm by Cedric Notredame
--
Computes all the pairwise alignments (PAM250, gop=-10, gep=-1) or use a multiple alignment.
Measure the %id (number id/number matches) of each pair
if a minimum identity min% is set: all the sequences with less than min% identity with ANY sequence in the set will be removed so that in the remaining set ALL the pairs of sequences have more than min% identity. The removal will stop uncompleted if the set becomes smaller than n.
Remove one of the two closest sequences until either n is reached or until all the sequences have less than max% identity.
return the new set.

My algorithm:
- sort by length
- starting with batches of size n sequences
	- compute pairwise similarity
	- for any sequences with sim > threshold, remove the one with highest mean similarity to others
'''

# Get species name from SMART FASTA header
def getSpeciesName(hdr):
	# [Cyanidioschyzon merolae strain 10D]
	return hdr.split('[')[-1].split(']')[0].strip()

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Prune alignment to remove similar sequences")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input FASTA filename")
	# Optional arguments
	parser.add_argument("--identity", dest="identity_threshold", type=float, default=0.99, help="maximum proportion identity allowed in alignment")
	parser.add_argument("--fraction-aligned", dest="fraction_aligned_threshold", type=float, default=0.9, help="maximum proportion identity allowed in alignment")
	parser.add_argument("--chunksize", dest="chunksize", type=int, default=None, help="number of sequences to initially consider in a group when assessing pairwise identity (defaults to all)")
	parser.add_argument("--anchor", dest="anchor", action='append', default=[], help="identifiers of sequences that must remain in the alignment")
	parser.add_argument("--one-per-species", dest="one_per_species", action="store_true", default=False, help="one sequence per species?")
	parser.add_argument("-g", "--debug", dest="debug", action='store_true', default=False, help="debug?")
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

	# Sort by length
	# Save indices so we can back them out later
	lengths = []
	for (xi, s) in enumerate(seqs):
		lengths.append((geneutil.truelen(s), xi))
	lengths.sort()
	original_seqs = seqs[:]
	original_headers = headers[:]
	seqs = [seqs[xi] for (L,xi) in lengths]
	headers = [headers[xi] for (L,xi) in lengths]

	# Accumulate sequences that must remain in alignment
	anchor_ids = options.anchor
	anchor_indices = []
	for (xi, h) in enumerate(headers):
		for a in anchor_ids:
			if a in h:
				anchor_indices.append(xi)
				info_outs.write("# Anchoring sequence {:s}\n".format(h))

	if options.debug:
		seqs = seqs[:20]
		headers = headers[:20]
	# Split up and 
	N = len(seqs)
	chunksize = N
	if not options.chunksize is None:
		chunksize = options.chunksize
	# similarities
	sims = sp.zeros((N,N))
	sims_hash = {}
	sequences_to_remove = []

	done = False
	while not done:
		num_chunks = int(math.floor(N/chunksize))
		chunks = [range(i*chunksize,min(i*chunksize+chunksize,N)) for i in range(num_chunks)]
		if N % chunksize > 0:
			chunks.append(range(chunks[-1][-1]+1, N))
		assert list(itertools.chain(*chunks)) == list(range(N))

		# With chunks in place, compute similarity
		max_similarity = 0.0
		max_pair = None
		for chunk in chunks:
			for xi in chunk[:-2]:
				for xj in range((xi+1),chunk[-1]):
					si = seqs[xi]
					sj = seqs[xj]
					# Don't bother to compute similarities for sequences whose lengths are different enough to violate
					# the requirement for fraction aligned
					li = lengths[xi][0]
					lj = lengths[xj][0]
					if options.debug:
						assert li == geneutil.truelen(si)
						assert lj == geneutil.truelen(sj)
					if min(li/float(lj), lj/float(li)) < options.fraction_aligned_threshold:
						continue
					# Don't bother to compute similarities for sequences we've already decided to remove
					if (xi not in sequences_to_remove) and (xj not in sequences_to_remove):
						if sims[xi][xj] == 0.0:
							sim = translate.compare(si,sj)
							sims[xi][xj] = sim.identity
							sims_hash[(xi,xj)] = sim
						else:
							sim = sims_hash[(xi,xj)]
						if min(sim.fraction_aligned_x, sim.fraction_aligned_y) >= options.fraction_aligned_threshold and max_similarity < sim.identity:
							max_similarity = sim.identity
							max_pair = (xi, xj)
			# We've identified the most-similar pair, and all similarities
			# Now remove any above a threshold identity
			sorted_ids = sorted([(v, k) for (k,v) in sims_hash.items()], reverse=True, key=lambda x: x[0].identity)
			for si in sorted_ids:
				(sim, (xi, xj)) = si
				if (xi not in sequences_to_remove) and (xj not in sequences_to_remove) and (sim.identity >= options.identity_threshold):
					# Try to remove one of these sequences
					# Prioritize
					# First check if one or more can be removed -- not anchored
					remaining_options = [x for x in [xi, xj] if not x in anchor_indices]
					if len(remaining_options)==1:
						# Remove this one sequence
						sequences_to_remove.append(remaining_options[0])
					elif len(remaining_options)==2:
						# Remove the sequence most similar to other sequences
						sim_i = [x for x in sims[xi] if x>0.0]
						avg_sim_i = sp.mean(sim_i)
						sim_j = [x for x in sims[xj] if x>0.0]
						avg_sim_j = sp.mean(sim_j)
						if avg_sim_j > avg_sim_i:
							sequences_to_remove.append(xj)
						else:
							sequences_to_remove.append(xi)

			info_outs.write("# Removing {} sequences at chunksize {} of {}\n".format(len(sequences_to_remove), chunksize, N))
		# Expand the chunksize and continue
		if chunksize == N:
			done = True
		if chunksize < N:
			chunksize = min(2*chunksize, N)

	# Grab only the sequences we're not removing; restore order
	new_indices = lengths[:]
	headers_to_remove = []
	for xi in sequences_to_remove:
		new_indices[xi] = None
		headers_to_remove.append(headers[xi])
	new_indices = [x for x in new_indices if not x is None]
	new_indices.sort(key=lambda x: x[1])
	seqs = []
	headers = []
	for (xi,entry) in enumerate(new_indices):
		(L, original_index) = entry
		seqs.append(original_seqs[original_index])
		headers.append(original_headers[original_index])
	# Did this work?
	for h in headers_to_remove:
		#print orthodbutil.headerDict(h)['organism_name']
		assert not (h in headers)

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

	# Write output
	biofile.writeFASTA(seqs, fasta_outs, headers)

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.fasta_out_fname is None:
		info_outs.write("# Wrote {} entries to {}\n".format(len(seqs), options.fasta_out_fname))
		outf.close()

