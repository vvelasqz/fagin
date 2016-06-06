# Cadmium

# Input

 - For each species
   - GFF file (must at least include gene models)
   - Full genome (GFF reference)
 - Search space file (Synder output)

# Pipeline

 - Identify target genes that overlap the search space.
 - Search the query protein against the overlapping target gene's coding sequences

# Output

 - Classification of each gene

# TODO

 - [x] Script to retrieve sample data for species of interest
 - [ ] Script to check input data
 - [ ] Script to prepare summarizing datasets for included genomes
 - [ ] R code for importing required data
 - [ ] Classify focal species gene against extra-species targets
