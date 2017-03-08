#! python

import sys, os, math, random, argparse, time
import util, biofile, muscle, translate, geneutil

'''
Prune out:
1) Sequences with X
2) Sequences that don't begin
Choose 
'''

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Prune sequences based on completeness")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input FASTA filename")
	# Optional arguments
	parser.add_argument("--remove-x", dest="remove_x", action='store_true', default=False, help="remove sequences containing X's?")
	parser.add_argument("--remove-non-met-start", dest="remove_non_met_start", action='store_true', default=False, help="remove sequences that don't start with Met?")
	parser.add_argument("--fasta-out", dest="fasta_out_fname", default=None, help="output FASTA filename")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	fasta_outs = util.OutStreams()

	# Start up output
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

	new_headers = []
	new_seqs = []
	for (h,s) in zip(headers,seqs):
		remove = False
		if options.remove_x:
			remove = 'X' in s
		if options.remove_non_met_start:
			remove = remove or s.replace('-','')[0] != 'M'
		if not remove:
			new_headers.append(h)
			new_seqs.append(s)


	# Write output
	n_written = len(new_seqs)
	biofile.writeFASTA(new_seqs, fasta_outs, new_headers)

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.fasta_out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.fasta_out_fname))
		outf.close()

