#
# Written by Jonas Juselius <jonas@iki.fi> Tue Dec  4 13:38:16 EET 2001
#


srcdir:=.
top_srcdir:=.
builddir:=.

include $(builddir)/make.config

basdir=${prefix}/basis

auxdir=$(top_srcdir)/config
mkinstalldirs=$(auxdir)/mkinstalldirs

libdirs:=lib_mod libsrint lib_cc
		 
progdirs:=getHF CPSd.SF EE-CIS.SF EE-CPSd.SF grCPSd.SF CISdoo.SD.SF EE-grCPSd.SF


ifneq ($(MPI_TYPE), none)
	progdirs+=pararun 
endif
ifeq ($(ENABLE_GIMIC), yes)
	progdirs+=cpdens
endif

progdirs+=pyscf2cc

dirs:=$(libdirs) $(progdirs)

.PHONY: all lib prog clean distclean 

all: lib prog

lib: 
	@mkdir -p lib bin
	@set -e; \
   for i in $(libdirs); do \
		if [ -d $$i ]; then \
		  if [ -f $$i/Makefile ]; then \
		    $(MAKE) -C $$i all ;\
		  fi; \
		fi; \
   done;

prog: lib
	@set -e; \
   for i in $(progdirs); do \
		if [ -d $$i ]; then \
		  if [ -f $$i/Makefile ]; then \
		    $(MAKE) -C $$i all ;\
		  fi; \
		fi; \
   done;

clean:
	@set -e; \
	 for i in $(dirs); do \
		if [ -d $$i ]; then \
		  if [ -f $$i/Makefile ]; then \
			$(MAKE) -C $$i $@ || exit 1; \
		  fi; \
		fi; \
	 done; 

realclean:
	@set -e; \
	for i in $(dirs); do \
		if [ -d $$i ]; then \
		  if [ -f $$i/Makefile ]; then \
			$(MAKE) -C $$i $@ || exit 1; \
		  fi; \
		fi; \
	done; 

distclean: realclean
	-rm config.h config.log config.status Makefile

basis_sets:=$(wildcard $(top_srcdir)/basis/*)
basis_sets:=$(filter-out $(top_srcdir)/basis/CVS, $(basis_sets))

install: install_bin install_basis

install_all: install_bin install_lib install_basis

install_bin:
	@$(mkinstalldirs) $(bindir)
	@set -e; \
		for i in $(dirs); do \
			if [ -d $$i ]; then \
			  if [ -f $$i/Makefile ]; then \
			    $(MAKE) -C $$i install || exit 1; \
			  fi; \
			fi; \
		done; 
	
install_lib:
	@$(mkinstalldirs) $(libdir)
	@set -e; \
		for i in $(dirs); do \
			if [ -d $$i ]; then \
			  if [ -f $$i/Makefile ]; then \
			    $(MAKE) -C $$i install || exit 1; \
			  fi; \
			fi; \
		done; 

install_basis:
	@$(mkinstalldirs) $(basdir)
ifneq ($(basis_sets),)
	@set -e; \
	for i in $(basis_sets); do \
		$(INSTALL_PROGRAM) $$i $(basdir); \
		echo "Installed $$i in $(basdir)"; \
	done
endif

PACKAGE_NAME:=CCSO
PACKAGE_TARNAME:=ccso
PACKAGE_BUGREPORT:=404227815@qq.com, wangf44@gmail.com
PACKAGE_VERSION:=1.0
ARCHIVE_NAME=$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)
ARCHIVE_TARNAME=$(PACKAGE_TARNAME)_$(PACKAGE_VERSION)

.PHONY: archive
archive:
	git-archive --format=tar --prefix=$(ARCHIVE_NAME)/ HEAD | gzip -c >\
		$(ARCHIVE_TARNAME).tgz
