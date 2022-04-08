
if [ ! -d ".venv" ];then
  python3 -m virtualenv .venv ;
  . .venv/bin/activate \
  && pip install planemo \
  ;
fi


. .venv/bin/activate \
&& planemo conda_install \
  --conda_prefix /tmp/mixmodel/ \
  . \
;

. .venv/bin/activate \
&& planemo test \
  --install_galaxy \
  --conda_dependency_resolution \
  --conda_prefix /tmp/mixmodel/ \
  --no_cleanup \
;
