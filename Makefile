all: build web

build:
	Rscript -e 'devtools::document()'
	Rscript -e 'devtools::build()'

check:
	Rscript -e 'devtools::check()'

web:
	Rscript -e 'pkgdown::build_site()'
