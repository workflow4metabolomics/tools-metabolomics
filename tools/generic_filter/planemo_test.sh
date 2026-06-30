planemo conda_init --conda_prefix /tmp/mc/
planemo conda_install --conda_prefix /tmp/mc/ .
planemo test --install_galaxy --conda_dependency_resolution --conda_prefix /tmp/mc/ --no_cleanup