# NOTE: This is not a very reliable Makefile. It is convenient for me, but
# needs a major overhaul before general application.

all:
	./1_prepare-inputs.sh
	[[ -d input/faa ]] || ./2_extract-fasta.sh
	[[ -d stat ]]      || ./3_gather-genome-summary-data.sh
	./4_get-search-intervals.sh
	./5_prepare-report.sh

.PHONY: clean
clean:
	rm -f log report.pdf
	rm -rf cache
	rm -rf input
