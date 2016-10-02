# NOTE: This is not a very reliable Makefile. It is convenient for me, but
# needs a major overhaul before general application.

all:
	./1_prepare-inputs.sh
	[[ -d input/faa ]]  || ./2_extract-fasta.sh
	[[ -d input/stat ]] || ./3_gather-genome-summary-data.sh
	./4_get-search-intervals.sh
	./5_prepare-report.sh

# Builds synder, recalculates search intervals, cleans house, makes tags
.PHONY: refresh
refresh:
	cd ../synder && make && make install PREFIX=${HOME} && make test
	rm -rf input/maps/db
	./4_get-search-intervals.sh
	make clean
	cd src/report && make deepclean && ctags .

.PHONY: archive
archive:
	./archive.sh

.PHONY: clean
clean:
	rm -f *log log report.pdf
	rm -rf cache
	rm -f gmon.out

.PHONY: deepclean
deepclean:
	rm -rf input
	rm -f log report.pdf
	rm -rf cache
	rm -f gmon.out
