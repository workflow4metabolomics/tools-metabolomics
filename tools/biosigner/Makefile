all:

test:
	$(MAKE) -C $@

planemo-venv/bin/planemo: planemo-venv
	. planemo-venv/bin/activate && pip install --upgrade pip setuptools
	. planemo-venv/bin/activate && pip install planemo

planemo-venv:
	virtualenv planemo-venv

planemolint: planemo-venv/bin/planemo
	. planemo-venv/bin/activate && planemo lint --no_xsd

planemotest: R_LIBS_USER=
planemotest: planemo-venv/bin/planemo
	. planemo-venv/bin/activate && planemo test --conda_dependency_resolution --galaxy_branch release_17.09

planemo-testtoolshed-diff: distplanemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_diff --shed_target testtoolshed

planemo-testtoolshed-update: distplanemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_update --check_diff --shed_target testtoolshed

planemo-toolshed-diff: distplanemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_diff --shed_target toolshed

planemo-toolshed-update: distplanemo-venv/bin/planemo
	. planemo-venv/bin/activate && cd $< && planemo shed_update --check_diff --shed_target toolshed

clean:
	$(MAKE) -C test $@
	$(RM) -r planemo-venv
	$(RM) -r planemotest.log

.PHONY: all clean test planemolint planemotest planemon-install planemo-toolshed-diff planemo-toolshed-update planemo-testtoolshed-diff planemo-testtoolshed-update testvenv
