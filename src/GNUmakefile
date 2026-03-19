ifneq ($(MAKELEVEL),0)

.DEFAULT_GOAL := all
.PHONY: force

%: force
	@exec $(MAKE) -f Makefile $@

force:

else

.PHONY: all clean mostlyclean distclean maintainer-clean check install python \
	RNAfold RNAcofold RNAsubopt

BUILDROOT ?= .build
SRC_BUILD := $(BUILDROOT)/src
BIN_BUILD := $(SRC_BUILD)/bin
PY_BUILD := $(BUILDROOT)/interfaces/Python

define require_configure
	@test -f "$(SRC_BUILD)/Makefile" || { \
		echo "Run ./configure in /Users/howardsalis/ViennaRNA/src first."; \
		exit 1; \
	}
endef

all clean mostlyclean check install:
	$(call require_configure)
	$(MAKE) -C $(SRC_BUILD) $@

python:
	$(call require_configure)
	$(MAKE) -C $(PY_BUILD) all

distclean maintainer-clean:
	@test -f "$(BUILDROOT)/Makefile" || { \
		echo "Nothing configured under $(BUILDROOT)"; \
		exit 0; \
	}
	$(MAKE) -C $(BUILDROOT) $@

RNAfold RNAcofold RNAsubopt:
	$(call require_configure)
	$(MAKE) -C $(BIN_BUILD) libhelpers.la $@

endif
