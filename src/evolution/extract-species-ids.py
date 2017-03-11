#! python

import sys, os, math, random, argparse
import util, biofile

def extractSpeciesName(text):
	res = h.split('[')[1].split(']')[0]
	return res

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Extract species names, and optionally FASTA file with species names as identifiers, from input FASTA")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename")
	# Optional arguments
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
	fname =os.path.expanduser(options.in_fname)
	if not os.path.isfile(fname):
	 	raise IOError("# Error: file {} does not exist".format(fname))
	with open(fname,'r') as inf:
	 	# Read a FASTA file?
	 	(headers, seqs) = biofile.readFASTA(inf)
	 	# Read a tab-delimited file?

	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('species','Species name','s')
	dout.describeHeader(data_outs)
	n_written = 0
	new_headers = []
	for h in headers:
		spec = extractSpeciesName(h)
		new_headers.append(spec)
		line = spec + '\n'
		data_outs.write(line)
		n_written += 1

	n_fasta_written = 0
	if not options.fasta_out_fname is None:
		fasta_outs = util.OutStreams(open(os.path.expanduser(options.fasta_out_fname),'w'))
		biofile.writeFASTA(seqs, fasta_outs, headers=new_headers)
		info_outs.write("# Wrote {} entries to {}\n".format(len(seqs), options.fasta_out_fname))

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

