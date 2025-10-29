PACKAGE=$(shell awk '/^Package: / { print $$2 }' DESCRIPTION)
VERSION=$(shell awk '/^Version: / { print $$2 }' DESCRIPTION)
TARBALL=$(PACKAGE)_$(VERSION).tar.gz

build:
	R CMD build .

build-lite:
	R CMD build --no-build-vignettes .

check:
	R CMD check "$(TARBALL)"

install:
	R CMD INSTALL --install-tests --html --example "$(TARBALL)"

test:
	parallel -j 8 --halt now,fail=1 Rscript ::: tests/test*.R

build-and-test: build-lite test

coverage:
	Rscript -e "library(covr); coverage <- package_coverage(path='.', type='tests'); print(coverage); to_cobertura(coverage, file='coverage.xml'); report(coverage, 'coverage.html', FALSE)"

clean:
	$(RM) $(TARBALL)
	$(RM) -r $(PACKAGE).Rcheck
	$(RM) -r coverage.* lib/

all: build check install

.PHONY: all clean build build-lite check install test build-and-test coverage
