.PHONY: all docs build check

all: for-commit

for-commit: docs build check

docs:
	rm -rf NAMESPACE man docs
	R -e 'library(devtools); document()'
	R -e 'pkgdown::build_site()'

build:
	rm -f ../slanter_*.tar.gz
	cd .. && R CMD build --resave-data slanter

check:
	R CMD check --as-cran ../slanter_*.tar.gz
