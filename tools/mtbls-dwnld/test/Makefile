ifeq ($(CONTINUOUS_INTEGRATION),true)
TESTTHAT_FLAGS=-gp
endif

all: mtbls-dwnld

mtbls-dwnld:
	../bash-testthat/testthat.sh $(TESTTHAT_FLAGS) test-mtbls-dwnld.sh

clean:
	$(RM) -r MTBLS* *.zip outputs

.PHONY: clean all
