export BIODB_CACHE_DIRECTORY=$(CURDIR)/lcms.biodb.cache
CONDA_PREFIX=$(HOME)/.plncnd
ifndef GALAXY_VERSION
	GALAXY_VERSION=18.09
endif

all:

test:
	$(MAKE) -C $@

planemo-venv/bin/planemo: planemo-venv
	. planemo-venv/bin/activate && pip install "pip>=7"
	. planemo-venv/bin/activate && pip install planemo

planemo-venv:
	virtualenv planemo-venv

planemolint: planemo-venv/bin/planemo
	. planemo-venv/bin/activate && planemo lint --no_xsd

planemotest: R_LIBS_USER=
planemotest: planemo-venv/bin/planemo
	. planemo-venv/bin/activate && planemo test --conda_dependency_resolution --conda_prefix "$(CONDA_PREFIX)" --galaxy_branch release_$(GALAXY_VERSION)

planemo-testtoolshed-diff: dist/lcmsmatching/ planemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_diff --shed_target testtoolshed

planemo-testtoolshed-update: dist/lcmsmatching/ planemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_update --check_diff --shed_target testtoolshed

planemo-toolshed-diff: dist/lcmsmatching/ planemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_diff --shed_target toolshed

planemo-toolshed-update: dist/lcmsmatching/ planemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_update --check_diff --shed_target toolshed

dist/lcmsmatching/:
	mkdir -p $@
	cp -r README.md lcmsmatching lcmsmatching.xml test-data $@

clean:
	$(MAKE) -C test $@
	$(RM) -r $(BIODB_CACHE_DIRECTORY)
	$(RM) -r dist
	$(RM) -r planemo-venv
	$(RM) -r planemotest.log
	$(RM) -r $(HOME)/.planemo
	$(RM) -r $(CONDA_PREFIX)
	$(RM) test_*-output.tsv
	$(RM) tool_test_output.*
	$(RM) *.zip

.PHONY: all clean test planemo-lint planemo-test planemon-install planemo-toolshed-diff planemo-toolshed-update planemo-testtoolshed-diff planemo-testtoolshed-update pkg
