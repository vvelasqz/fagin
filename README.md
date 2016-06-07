# Cadmium

# Input

 For information on inputs, call `1_prepare-inputs.sh -H`

 In short, the following is required

 - Phylogeny for all included species
 - Name of the focal species
 - Synteny map for the focal species and each other species
 - For each species
   - GFF file (must at least include gene models)
   - Full genome (GFF reference)

# Pipeline

 - Identify target genes that overlap the search space.
 - Search the query protein against the overlapping target gene's coding sequences

# Output

 The output of this pipeline is a set of text files and an HTML report
 including the following:

 - The phylogenetic tree of all species in the pipeline
 - Genome summary for each species 
 - Overall statistics for classifications
 - Visualizations of overall statistics
 - Origin page for each orphan gene

# TODO

 - [x] Script to retrieve sample data for species of interest
 - [x] Script to check input data
 - [ ] Script to extract relevant fasta files from GFF and genomes 
 - [ ] Script to get search intervals using Synder
 - [ ] Script to prepare summarizing datasets for included genomes
 - [ ] R code for importing and iterating through all data
 - [ ] R code to classify focal species genes against each search interval
 - [ ] R code to merge search interval classifications
 - [ ] R code to get final gene class from species class vector
 - [ ] R code to statistically summarize results
 - [ ] R code to visualize results
 - [ ] R code to prepare report for each gene
 - [ ] R code to prepare HTML webpage, merging all results
