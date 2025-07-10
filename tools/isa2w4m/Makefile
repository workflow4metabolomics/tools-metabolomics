SCRIPT_NAME=isa2w4m.py
REPOS_NAME=isa2w4m
TOOL_XML=isa2w4m.xml
PYVER37=3.7.13
PYVER38=3.8.12
PYVER39=3.9.9
PYVER=$(PYVER38) $(PYVER39)
PYENV:=$(shell command -v pyenv)
ifeq (,$(PYENV))
PYTHON=python3 -s
else
PYTHON=$(PYENV) exec python3 -s
endif
VM = isa2w4m-archlinux
CONDA_DIR=$(HOME)/plnmconda
export TMPDIR=$(HOME)/plnmtmp
PLANEMO_DIR=$(HOME)/plnmws

all:

test: $(addprefix test,37 38 39)

testvm:
	teston -rD -f Makefile -f isa2w4m.py -f isa2w4m.xml -f requirements.txt -f tests/test-isa2w4m -f test-data -c gnumake -c pyenv -c unzip -c git -t test39 archlinux/archlinux

pyver:
ifneq (,$(PYENV))
	for ver in $(PYVER) ; do pyenv install -s $$ver ; done
	pyenv local $(PYVER)
endif
	echo "Using $$($(PYTHON) --version)."

define test_on_pyver
test$(1): export PYENV_VERSION=$$(PYVER$(1))
test$(1): pyver requirements.txt
	$(PYTHON) -m venv test$(1)-venv
	. test$(1)-venv/bin/activate && test "Python $$(PYENV_VERSION)" = "$$$$(python3 --version)"
	. test$(1)-venv/bin/activate && python3 -s -m pip install --no-cache-dir --upgrade pip
	. test$(1)-venv/bin/activate && python3 -s -m pip install --no-cache-dir -r requirements.txt
	. test$(1)-venv/bin/activate && python3 -s -c 'import isatools'
	. test$(1)-venv/bin/activate && tests/test-isa2w4m
endef

$(foreach ver,37 38 39,$(eval $(call test_on_pyver,$(ver))))

$(TMPDIR):
	mkdir -p "$@"

# TODO Check for presence of ~/.local/lib/python* folders before running planemo.
define ptest_on_pyver
planemo$(1)-venv: export PYENV_VERSION=$$(PYVER$(1))
planemo$(1)-venv: pyver
	$(PYTHON) -m venv $$@
	. planemo$(1)-venv/bin/activate && test "Python $$(PYENV_VERSION)" = "$$$$(python3 --version)"

planemo$(1)-venv/bin/planemo: planemo$(1)-venv
	. planemo$(1)-venv/bin/activate && python3 -s -m pip install planemo

plint$(1): planemo-venv/bin/planemo
	. planemo$(1)-venv/bin/activate && planemo --directory $(PLANEMO_DIR)$(1) lint $(TOOL_XML)

ptest$(1): planemo$(1)-venv/bin/planemo $(TMPDIR)
	. planemo$(1)-venv/bin/activate && planemo --directory $(PLANEMO_DIR)$(1) test --conda_prefix "$(CONDA_DIR)$(1)" --conda_dependency_resolution --galaxy_python_version $(subst 3,3.,$(1)) $(TOOL_XML)

pttsdiff$(1): dist/$(REPOS_NAME)/ planemo$(1)-venv/bin/planemo $(HOME)/.planemo.yml
	@echo "Check difference with testtoolshed version."
	. planemo$(1)-venv/bin/activate && cd $$< && planemo shed_diff --shed_target testtoolshed

pttsupdate$(1): dist/$(REPOS_NAME)/ planemo$(1)-venv/bin/planemo $(HOME)/.planemo.yml
	@echo "Push new version to testtoolshed."
	. planemo$(1)-venv/bin/activate && cd $$< && planemo shed_update --check_diff --skip_metadata --shed_target testtoolshed

ptsdiff$(1): dist/$(REPOS_NAME)/ planemo$(1)-venv/bin/planemo $(HOME)/.planemo.yml
	@echo "Check difference with toolshed version."
	. planemo$(1)-venv/bin/activate && cd $$< && planemo shed_diff --shed_target toolshed

ptsupdate$(1): dist/$(REPOS_NAME)/ planemo$(1)-venv/bin/planemo $(HOME)/.planemo.yml
	@echo "Push new version to toolshed."
	. planemo$(1)-venv/bin/activate && cd $$< && planemo shed_update --check_diff --skip_metadata --shed_target toolshed
endef

$(foreach ver,37 38 39,$(eval $(call ptest_on_pyver,$(ver))))

plint: plint39

ptest: ptest39

dist/$(REPOS_NAME)/:
	mkdir -p $@
	cp -Lr README.md $(SCRIPT_NAME) $(TOOL_XML) test-data $@

pttsdiff: pttsdiff39

ptsdiff: ptsdiff39

pttsupdate: pttsupdate39

ptsupdate: ptsupdate39

clean:
	$(RM) -r *-venv $(HOME)/plnm*
	$(RM) tool_test_output.*
	$(RM) -r dist
	$(RM) -r tests/output
	$(RM) .python-version
	$(RM) -r $(TMPDIR)

.PHONY: all clean test plint ptest ptesttoolshed_diff ptesttoolshed_update ptoolshed_diff ptoolshed_update
