# In shell:
# > echo $PYTHONEXE
# python
PY = $$PYTHONEXE
RSCRIPT = $$RSCRIPT
DATA = ../data
FIG = ../figures
TMP = ../../../../tmp
LIB = ../../base/src

PABP_FINAL = $(DATA)/alignment/pabp-smart-orthologs-aligned-trimmed.fa
BASE_FASTA = $(DATA)/pabp-smart-orthologs-aligned-trimmed.fa
PHYLO_DATA = $(DATA)
SORTED_FASTA = $(PHYLO_DATA)/pabp-sorted.fa

# This command takes a raw unaligned FASTA file, removes partial sequences,
# aligns the sequences, then trims to make an alignment with 
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
	"$$RSCRIPT" normalize-species.R $(PHYLO_DATA)/species-ids.txt $(PHYLO_DATA)/class

build-tree:
	$(PY) build-tree.py $(PHYLO_DATA)/class/guide.txt --out $(PHYLO_DATA)/pabp-tree.txt

reextract-species:
	$(PY) species-from-classifications.py $(PHYLO_DATA)/class/guide.txt --binomial --out $(PHYLO_DATA)/species-binomial-only.txt

update-fasta:
	$(PY) update-fasta-with-tree.py $(BASE_FASTA) $(PHYLO_DATA)/pabp-tree.txt $(PHYLO_DATA)/class/guide.txt \
	--generate-short-ids --fasta-out $(SORTED_FASTA)

pabpretty:
#	$(PY) pretty-orthodb-alignment.py $(SGAGG_DATA)/alignment/pab1-orthologs-orthodb-trimmed.fa --one-per-species --out $(SGAGG_DATA)/alignment/tmp.log --fasta-out $(SGAGG_DATA)/alignment/pab1-orthologs-orthodb-pretty.fa
#	$(PY) pretty-orthodb-alignment.py $(SGAGG_DATA)/alignment/pab1-orthologs-orthodb-trimmed.fa --one-per-species --gap-threshold 2 --out $(SGAGG_DATA)/alignment/tmp.log --fasta-out $(SGAGG_DATA)/alignment/pab1-orthologs-display.fa
	$(PY) pretty-smart-alignment.py $(SORTED_FASTA) --gap-threshold 0.004 \
		--out $(PHYLO_DATA)/pretty.log --fasta-out $(PHYLO_DATA)/pabp-sorted-display.fa


# Given a tree 
extract-region:
	$(PY) $(LIB)/extract-aligned-region.py $(SORTED_FASTA) --query cerevisiae \
		--start-sequence YQQATAAAAAAAAGMP --end-sequence ANDNNQFYQ \
		--fasta-out $(PHYLO_DATA)/hypr.fa
	$(PY) $(LIB)/extract-aligned-region.py $(SORTED_FASTA) --query cerevisiae \
		--start-sequence YQQATAAAAAAAAGMP --end-sequence ANDNNQFYQ --exclude \
		--fasta-out $(PHYLO_DATA)/non-hypr.fa

filter-extracted:
	$(PY) filter-jointly.py --

compute-hydrophobicity:
	$(PY) compute-hydrophobicity.py $(PHYLO_DATA)/hypr.fa --tree $(PHYLO_DATA)/pabp-tree.txt \
		--scale Hopp-Woods --randomize --reps 1 --out $(PHYLO_DATA)/pabp-hypr-all-hydro.txt
	$(PY) compute-hydrophobicity.py $(PHYLO_DATA)/hypr.fa --tree $(PHYLO_DATA)/pabp-tree.txt \
		--aas ILMVA --scale Hopp-Woods --randomize --reps 1 --out $(PHYLO_DATA)/pabp-hypr-ILMVA-hydro.txt
	$(PY) compute-hydrophobicity.py $(PHYLO_DATA)/hypr.fa --tree $(PHYLO_DATA)/pabp-tree.txt \
		--phylo-node Chordata --scale Hopp-Woods --randomize --reps 1 --out $(PHYLO_DATA)/pabp-hypr-chordata-all-hydro.txt
	$(PY) compute-hydrophobicity.py $(PHYLO_DATA)/hypr.fa --tree $(PHYLO_DATA)/pabp-tree.txt \
		--phylo-node Fungi --scale Hopp-Woods --randomize --reps 1 --out $(PHYLO_DATA)/pabp-hypr-fungi-all-hydro.txt
	$(PY) compute-hydrophobicity.py $(PHYLO_DATA)/hypr.fa --tree $(PHYLO_DATA)/pabp-tree.txt \
		--phylo-node Chordata --aas ILMVA --scale Hopp-Woods --randomize --reps 1 --out $(PHYLO_DATA)/pabp-hypr-chordata-ILMVA-hydro.txt
	$(PY) compute-hydrophobicity.py $(PHYLO_DATA)/hypr.fa --tree $(PHYLO_DATA)/pabp-tree.txt \
		--phylo-node Fungi --aas ILMVA --scale Hopp-Woods --randomize --reps 1 --out $(PHYLO_DATA)/pabp-hypr-fungi-ILMVA-hydro.txt

prop:
	$(PY) $(LIB)/protprop.py --in $(PHYLO_DATA)/hypr.fa --aas ILMVA --out $(PHYLO_DATA)/ILMVA-freqs.txt
	$(PY) $(LIB)/protprop.py --in $(PHYLO_DATA)/hypr.fa --aas all --out $(PHYLO_DATA)/hypr-freqs.txt
	$(PY) $(LIB)/protprop.py --in $(PHYLO_DATA)/non-hypr.fa --aas all --out $(PHYLO_DATA)/non-hypr-freqs.txt
