all:
	Rscript -e 'roxygen2::roxygenise(".")'
	R CMD INSTALL .
	Rscript -e 'tinytest::test_package("stochgen")'

clean:
	rm -f src/*.o src/*.so *.tar.gz