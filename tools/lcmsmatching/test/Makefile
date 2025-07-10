export BIODB_CACHE_DIRECTORY=$(CURDIR)/../lcms.biodb.cache

all: clean
	../bash-testthat/testthat.sh test-lcmsmatching

clean:
	$(RM) *-output.tsv
	$(RM) *-output-peaks.tsv
	$(RM) *.html *-main.tsv *-peak.tsv *.txt *.log
	$(RM) -r output

.PHONY: all clean
