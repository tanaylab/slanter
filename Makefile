.PHONY: all docs build check

all: docs

for-commit: docs build check

docs:
	rm -rf NAMESPACE man
	R -e 'library(devtools); document()'

build:
	rm -f ../slanter_*.tar.gz
	cd .. && R CMD build slanter

check:
	R CMD check --as-cran ../slanter_*.tar.gz
