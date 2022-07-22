
TOOLNAME=MS2snoop

######### BEGIN  SCRIPT INFO #########

# SCRIPT_RUNNER=perl
# SCRIPT_RUNNER=python
# SCRIPT_RUNNER=python3
SCRIPT_RUNNER=Rscript

ifeq ($(SCRIPT_RUNNER),Rscript)
	SCRIPT_FILE=$(TOOLNAME).R
endif
ifeq ($(SCRIPT_RUNNER),python)
	SCRIPT_FILE=$(TOOLNAME).py
endif
ifeq ($(SCRIPT_RUNNER),python3)
	SCRIPT_FILE=$(TOOLNAME).py
endif
ifeq ($(SCRIPT_RUNNER),perl)
	SCRIPT_FILE=$(TOOLNAME).pl
endif

override SCRIPT_OPTIONS := \
  -c $(TMP_DIR)/test-data/compounds_pos.csv \
  -f $(TMP_DIR)/test-data/peaklist_fragments.tsv \
  -p $(TMP_DIR)/test-data/peaklist_precursors.tsv \
  --ionization="[M+H]+" \
  -o $(TMP_DIR)/out.tsv \
  $(SCRIPT_OPTIONS)
#   --mzdecimal=6 \

######### ENDOF SCRIPT INFO #########


XML_NAME=$(TOOLNAME).xml
CONDA_ENV_NAME=/tmp/$(TOOLNAME)-envs/$(TOOLNAME)-conda-env
VIRTUAL_ENV_NAME=/tmp/$(TOOLNAME)-envs/$(TOOLNAME)-venv
CONDA_LOCKFILE=.conda-lock.yml
TOOL_ENV=$(CONDA_ENV_NAME)/envs/$(TOOLNAME)


CONDA_ENV_PATH := $(CONDA_ENV_NAME)
VIRTUAL_ENV_PATH := $(VIRTUAL_ENV_NAME)
# CONDA_ENV_PATH := $(TMP_DIR)/$(CONDA_ENV_NAME)
# VIRTUAL_ENV_PATH := $(TMP_DIR)/$(VIRTUAL_ENV_NAME)

CONDA_BIN := $(CONDA_ENV_PATH)/bin/conda
CONDA_ACTIVATE := $(CONDA_ENV_PATH)/bin/activate


CONDA_DOWNLOAD_DOMAIN:=repo.anaconda.com
CONDA_DOWNLOAD_FILE:=Miniconda3-py37_4.12.0-Linux-x86_64.sh
CONDA_DOWNLOAD_PATH:=miniconda/$(CONDA_DOWNLOAD_FILE)
CONDA_DOWNLOAD_URL:=https://$(CONDA_DOWNLOAD_DOMAIN)/$(CONDA_DOWNLOAD_PATH)
CONDA_INSTALL_SCRIPT=/tmp/install_conda.sh


## Always creates a temp directory ; it will be removed at the end
ifeq ($(TMP_DIR),)
location = $(CURDIR)
self := $(location)/Makefile
%:
	@TMP_DIR=$(shell mktemp -d); \
	echo $${TMP_DIR} ; \
	trap 'rm -rf "$${TMP_DIR}"' EXIT; \
	$(MAKE) \
		-f $(self) \
		--no-print-directory \
		TMP_DIR=$${TMP_DIR} \
		CWD=$(pwd) \
		$@ \
	;
else


%: $(TOOL_ENV)


test:$(TOOL_ENV) link_files
	export PATH=$(TOOL_ENV):$(TOOL_ENV)/bin:${PATH} \
	&& cd $(TMP_DIR) \
	&& $(SCRIPT_RUNNER) ./$(SCRIPT_FILE) $(SCRIPT_OPTIONS) \
	&& echo "" \
	&& echo "-------------------" \
	&& echo "Find the results in $(TMP_DIR)." \
	&& echo "Press enter to terminate the script and erase temp files..." \
	&& read x \
	;

## The only good way to iterate on files is to not rely on `find` nor `ls`!
link_files:
	@for f in ./*; do \
		ln -s $(shell pwd)/$$f $(TMP_DIR)/ ; \
	done \
	;

$(TOOL_ENV): $(CONDA_LOCKFILE)
	@$(CONDA_ACTIVATE) && $(CONDA_BIN)-lock \
		install \
		--name $(TOOLNAME) \
		--conda $(CONDA_BIN) \
		$(CONDA_LOCKFILE) \
	;

$(CONDA_LOCKFILE): $(VIRTUAL_ENV_PATH) $(CONDA_ENV_PATH)
	@if [ ! -e "$(CONDA_LOCKFILE)" ];then \
		"$(VIRTUAL_ENV_PATH)"/bin/pip3 install conda-lock \
		&& if [ "$(shell which conda-lock)" = "" ];then \
			$(CONDA_BIN) install \
				-y \
				-c conda-forge \
				conda-lock \
			; \
		fi \
		&& $(CONDA_BIN)-lock \
			-f environment.yml \
			-p linux-64 \
			--lockfile "$(CONDA_LOCKFILE)" \
		; \
	fi \
	;

$(CONDA_ENV_PATH):
	@if [ ! -d "$(CONDA_ENV_PATH)" ];then \
		echo "conda env not created yet, creating..." \
		&& if [ ! -e $(CONDA_INSTALL_SCRIPT) ];then \
			wget \
				-O $(CONDA_INSTALL_SCRIPT) \
				$(CONDA_DOWNLOAD_URL) \
			; \
		fi \
		&& bash $(CONDA_INSTALL_SCRIPT) -b -p "$(CONDA_ENV_PATH)" \
		&& $(CONDA_BIN) install -y -n base -c conda-forge conda-lock \
		; \
	fi \
	;

$(VIRTUAL_ENV_PATH):
	@if [ ! -d "$(VIRTUAL_ENV_PATH)" ];then \
		python3 -m virtualenv "$(VIRTUAL_ENV_PATH)" ; \
	fi \
	;


endif ## ends temp directory rules