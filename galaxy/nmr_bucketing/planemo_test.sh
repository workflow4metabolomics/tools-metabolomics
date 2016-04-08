planemo conda_init
planemo conda_install .
planemo test --install_galaxy --conda_dependency_resolution --galaxy_branch "dev"

#All 1 test(s) executed passed.
#nmr_bucketing[0]: passed


planemo shed_test -t testtoolshed --install_galaxy --galaxy_branch "dev"

#All 1 test(s) executed passed.
#testtoolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_bucketing/NmrBucketing/1.0.1[0]: passed
