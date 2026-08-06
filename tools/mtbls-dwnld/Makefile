TOOL_NAME=mtbls-dwnld
REPOS_NAME=mtblsdwnld
TOOL_XML=mtbls-dwnld.xml
PYTHON_VERSION=3.6.10

all: isaslicer.py

isaslicer.py: isatools-galaxy/tools/isatools/isaslicer.py isaslicer.deps
	ln -sf $<

test: test-venv isaslicer.py
	. test-venv/bin/activate && $(MAKE) -C $@ ; deactivate

isaslicer.deps: test-venv
	. test-venv/bin/activate && pip install pandas isatools ; deactivate
	. test-venv/bin/activate && python -c 'import pandas;import isatools' ; deactivate

%-venv: install_python
	PYENV_VERSION=$(PYTHON_VERSION) python3 -m venv $@
	. $@/bin/activate && pip install --upgrade pip ; deactivate

test-venv/bin/planemo: test-venv
	. test-venv/bin/activate && pip install planemo ; deactivate

install_python:
	pyenv install -s $(PYTHON_VERSION)

plint: test-venv/bin/planemo
	. test-venv/bin/activate && planemo lint $(TOOL_XML) ; deactivate

ptest: test-venv/bin/planemo isaslicer.py
	. test-venv/bin/activate && planemo test --conda_dependency_resolution --galaxy_branch release_20.09 $(TOOL_XML) ; deactivate

dist/$(REPOS_NAME)/: isaslicer.py
	mkdir -p $@
	cp -Lr README.md $(TOOL_NAME) $(TOOL_NAME).xml test-data isaslicer.py $@

ptesttoolshed_diff: dist/$(REPOS_NAME)/ test-venv/bin/planemo
	. test-venv/bin/activate && cd $< && planemo shed_diff --shed_target testtoolshed ; deactivate

ptesttoolshed_update: dist/$(REPOS_NAME)/ test-venv/bin/planemo
	. test-venv/bin/activate && cd $< && planemo shed_update --check_diff --skip_metadata --shed_target testtoolshed ; deactivate

ptoolshed_diff: dist/$(REPOS_NAME)/ test-venv/bin/planemo
	. test-venv/bin/activate && cd $< && planemo shed_diff --shed_target toolshed ; deactivate

ptoolshed_update: dist/$(REPOS_NAME)/ test-venv/bin/planemo
	. test-venv/bin/activate && cd $< && planemo shed_update --check_diff --skip_metadata --shed_target toolshed ; deactivate

clean:
	$(MAKE) -C test $@
	$(RM) -r $(HOME)/.planemo
	$(RM) -r test-venv
	$(RM) tool_test_output.*
	$(RM) -r dist
	$(RM) isaslicer.py

.PHONY:	all clean test plint ptest ptesttoolshed_diff ptesttoolshed_update ptoolshed_diff ptoolshed_update install_python
