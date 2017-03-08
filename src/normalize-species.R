# Goal: turn species names pulled from a FASTA file into
# species names valid for NCBI Taxonomy

library(taxize)


fetchClassifications <- function(species_names) {
	res <- classification(species_names, db = 'ncbi')
	res
}

saveClassifications <- function(classifications, outdir) {
	species_names <- names(classifications)
	# Build filenames -- we'll use them later.
	fnames <- sapply(seq_along(species_names), function(xi) {
		species <- species_names[xi]
		elems <- strsplit(species, ' ', fixed=T)[[1]]
		fname <- paste(paste(elems[1], elems[2], xi,sep='-'), ".txt", sep='')
		if (is.null(nrow(classifications[[xi]]))){ fname <- NA }
		fname
		})

	if (!dir.exists(outdir)) {
		dir.create(outdir, recursive=TRUE)
	}

	out.res <- lapply(seq_along(classifications), function(xi, n){
		# Each entry is a data frame
		classif <- classifications[[xi]]
		species <- n[[xi]]
		fname <- fnames[[xi]]
		if (!is.null(nrow(classif))) { # NA returned if no hit found
			outname <- paste(outdir, fname, sep='/');
			cat("# Writing to", outname, '\n');
			write.table(classif, file=outname, sep='\t', row.names=F, quote=F);
		}
		species
		}, n=names(classifications))
	updated.species <- lapply(seq_along(classifications), function(xi){
		# Each entry is a data frame
		classif <- classifications[[xi]]
		res <- NA
		if (!is.null(nrow(classif))) { # NA returned if no hit found
			res <- classif[nrow(classif),'name']
		}
		res
		})
	guide_table <- data.frame(id=seq_along(classifications), species=species_names, updated.species=unlist(updated.species), filename=fnames)
	guide.fname <- paste(outdir, "guide.txt", sep='/')
	write.table(guide_table, file=guide.fname, row.names=F, quote=F, sep='\t')
	cat("# Wrote guide information to", guide.fname, '\n');
}

args <- commandArgs(TRUE)
species_fname <- args[1]
outdir <- args[2]

x <- read.table(species_fname, header=F, sep='\t')
cls <- fetchClassifications(x$V1)
saveClassifications(cls, outdir)
