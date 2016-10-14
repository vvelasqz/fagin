
# ===================================================================
# REQUIRED INPUT
# ===================================================================

FOCAL_SPECIES <- 'Arabidopsis_thaliana'




# ===================================================================
# The defaults below are usable for a first pass
# ===================================================================

# Base p-value cutoffs (these will be adjusted for multiple testing
# query protein versus target gene alignments
PROT2PROT_PVAL <- 0.05

# query protein versus all SI translated ORFs
PROT2ALLORF_PVAL <- 0.05

# query protein versus translated ORFs from spliced transcripts
PROT2TRANSORF_PVAL <- 0.05

# query genes versus entire SI (nucleotide match)
DNA2DNA_PVAL <- 0.05

# Number of simulations
PROT2PROT_NSIMS <- 1e3
PROT2ALLORF_NSIMS <- 1e3
PROT2TRANSORF_NSIMS <- 1e3
DNA2DNA_MAXSPACE <- 1e8
INDEL_THRESHOLD <- 0.25




# ===================================================================
# You shouldn't change anything below
# The correct configuration will automatically be produced
# ===================================================================

HOME          = "/home/shoggoth/src/git/fagin"
FAA_DIR       = "/home/shoggoth/src/git/fagin/input/faa"
GFF_DIR       = "/home/shoggoth/src/git/fagin/input/gff"
SYN_DIR       = "/home/shoggoth/src/git/fagin/input/syn"
GENE_DIR      = "/home/shoggoth/src/git/fagin/input/gene"
GENOME_DIR    = "/home/shoggoth/src/git/fagin/input/fna"
SI_DIR        = "/home/shoggoth/src/git/fagin/input/maps"
SPECIES_FILE  = "/home/shoggoth/src/git/fagin/input/species"
SCAFLEN       = '/home/shoggoth/src/git/fagin/input/stat/scaffold-lengths.tab'
NSTRINGS      = '/home/shoggoth/src/git/fagin/input/stat/nstrings.tab'
ORPHAN_LIST   = "/home/shoggoth/src/git/fagin/input/orphan-list.txt"
KBCOMP        = '/home/shoggoth/src/git/fagin/input/stat/kb-composition.tab'
ORFGFF        = "/home/shoggoth/src/git/fagin/input/orf-gff"
ORFFAA        = "/home/shoggoth/src/git/fagin/input/orf-faa"
TRANS_ORF     = "/home/shoggoth/src/git/fagin/input/trans-orf"
CACHE         = "/home/shoggoth/src/git/fagin/cache"
DECISION_TREE = "/home/shoggoth/src/git/fagin/etc/cds.yaml"
TREE          = "/home/shoggoth/src/git/fagin/input/tree"
