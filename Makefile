# NOTE: This is not a very reliable Makefile. It is convenient for me, but
# needs a major overhaul before general application.

all:
	src/prologue/1_prepare-inputs.sh
	[[ -d input/faa ]]  || src/prologue/2_extract-fasta.sh
	[[ -d input/stat ]] || src/prologue/3_gather-genome-summary-data.sh
	src/prologue/4_get-search-intervals.sh
	src/prologue/5_prepare-report.sh

.PHONY: archive
archive:
	src/epilogue/archive.sh

.PHONY: clean
clean:
	rm -f *log log report.pdf
	rm -rf cache
	rm -f gmon.out
	rm -f preconfig.sh
	rm -f runconfig.R

.PHONY: distclean
distclean:
	rm -rf input
	rm -f log report.pdf
	rm -rf cache
	rm -f gmon.out
