SCRIPT_NAME=extract-from-isa
REPOS_NAME=isaextractor
TOOLS_XML=$(wildcard isa2*.xml)

all:

test:
	make -C $@

planemo-venv/bin/planemo: planemo-venv
	. planemo-venv/bin/activate && pip install --upgrade pip setuptools
	. planemo-venv/bin/activate && pip install planemo

planemo-venv:
	virtualenv -p python2.7 $@

plint: planemo-venv/bin/planemo
	. planemo-venv/bin/activate && planemo lint $(TOOLS_XML)

ptest: planemo-venv/bin/planemo
	. planemo-venv/bin/activate && planemo test --conda_dependency_resolution --galaxy_branch release_19.01 $(TOOLS_XML)

dist/$(REPOS_NAME)/:
	mkdir -p $@
	cp -Lr README.md $(SCRIPT_NAME) $(TOOLS_XML) test-data $@

ptesttoolshed_diff: dist/$(REPOS_NAME)/ planemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_diff --shed_target testtoolshed

ptesttoolshed_update: dist/$(REPOS_NAME)/ planemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_update --check_diff --shed_target testtoolshed

ptoolshed_diff: dist/$(REPOS_NAME)/ planemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_diff --shed_target toolshed

ptoolshed_update: dist/$(REPOS_NAME)/ planemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_update --check_diff --shed_target toolshed

clean:
	$(RM) -r $(HOME)/.planemo
	$(RM) -r planemo-venv
	make -C test $@
	$(RM) tool_test_output.*
	$(RM) -r dist

.PHONY:	all clean test planemolint planemotest
