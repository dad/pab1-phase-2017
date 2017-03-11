PY = python
DATA = ../../data
LIB = ../../../base/src

PHYLO_DATA = $(DATA)
PABP_FINAL = $(PHYLO_DATA)/alignment/pabp-smart-orthologs-aligned-trimmed.fa
BASE_FASTA = $(PHYLO_DATA)/pabp-smart-orthologs-aligned-trimmed.fa
SORTED_FASTA = $(PHYLO_DATA)/pabp-sorted.fa

all: filter-and-align create-sorted extract-region calculate-frequencies

# This command takes a raw unaligned FASTA file, removes partial sequences,
# aligns the sequences, then trims to make an alignment with one sequence per species.
filter-and-align: prune align trim
prune:
	$(PY) prune-fasta.py $(PHYLO_DATA)/smart-pabp-orthologs.fa --remove-x --remove-non-met-start \
	--fasta-out $(PHYLO_DATA)/smart-pabp-orthologs-pruned.fa

align:
	$(PY) $(LIB)/muscle.py $(PHYLO_DATA)/smart-pabp-orthologs-pruned.fa --out $(PHYLO_DATA)/pabp-smart-orthologs-aligned.fa

trim:
	$(PY) trim-alignment.py $(PHYLO_DATA)/pabp-smart-orthologs-aligned.fa --anchor S288c --anchor PABP1_HUMAN \
	--anchor ENSP00000313007 --identity 0.95 --fraction-aligned 0.95 --one-per-species --out $(PHYLO_DATA)/tmp.log \
	--fasta-out $(BASE_FASTA)

# This command takes an alignment in FASTA format with species names
# given in the headers in [species name] format, and builds a phylogenetic
# tree (without distances), then prunes and sorts the FASTA alignment
# according to the tree.
create-sorted: extract-species fetch-classifications build-tree update-fasta

extract-species:
	$(PY) extract-species-ids.py $(BASE_FASTA) --out $(PHYLO_DATA)/species-ids.txt

fetch-classifications:
	Rscript normalize-species.R $(PHYLO_DATA)/species-ids.txt $(PHYLO_DATA)/class

build-tree:
	$(PY) build-tree.py $(PHYLO_DATA)/class/guide.txt --out $(PHYLO_DATA)/pabp-tree.txt

reextract-species:
	$(PY) species-from-classifications.py $(PHYLO_DATA)/class/guide.txt --binomial --out $(PHYLO_DATA)/species-binomial-only.txt

update-fasta:
	$(PY) update-fasta-with-tree.py $(BASE_FASTA) $(PHYLO_DATA)/pabp-tree.txt $(PHYLO_DATA)/class/guide.txt \
	--generate-short-ids --fasta-out $(SORTED_FASTA)

# This command produces an alignment pruned of sequences that produce large alignment-spanning gaps.
# Useful for display.
pabpretty:
	$(PY) pretty-smart-alignment.py $(SORTED_FASTA) --gap-threshold 0.004 \
		--out $(PHYLO_DATA)/pretty.log --fasta-out $(PHYLO_DATA)/pabp-sorted-display.fa


# Given an alignment, extract a chunk of columns where the query species starts and ends
# with the given sequence.
extract-region:
	$(PY) $(LIB)/extract-aligned-region.py $(SORTED_FASTA) --query cerevisiae \
		--start-sequence YQQATAAAAAAAAGMP --end-sequence ANDNNQFYQ \
		--fasta-out $(PHYLO_DATA)/pab-pdomain.fa
	$(PY) $(LIB)/extract-aligned-region.py $(SORTED_FASTA) --query cerevisiae \
		--start-sequence YQQATAAAAAAAAGMP --end-sequence ANDNNQFYQ --exclude \
		--fasta-out $(PHYLO_DATA)/pab-non-pdomain.fa

# Calculate amino acid frequencies.
calculate-frequencies:
	$(PY) $(LIB)/protprop.py --in $(PHYLO_DATA)/pab-pdomain.fa --aas all --out $(PHYLO_DATA)/pab-pdomain-freqs.txt
	$(PY) $(LIB)/protprop.py --in $(PHYLO_DATA)/pab-non-pdomain.fa --aas all --out $(PHYLO_DATA)/pab-non-pdomain-freqs.txt
