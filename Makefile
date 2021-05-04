VERSION = $(shell cat DESCRIPTION | grep Version | sed 's/Version: //')
R       = R -q -e

test:
  $(R) "devtools::test()"

check:
  $(R) "devtools::check()"

build_src:
  $(R)"devtools::build()"

build_doc:
  $(R) "pkgdown::build_site()"

document:
  $(R) "devtools::document(roclets = c('rd', 'collate', 'namespace', 'vignette'))"

install: document
  $(R) "devtools::install(upgrade = 'never')"
