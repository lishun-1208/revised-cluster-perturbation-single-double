#
# Written by Jonas Juselius <jonas@iki.fi> Tue Dec  4 11:41:04 EET 2001
#
#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
builddir:=..

include $(builddir)/make.config

INST_PROGS:=
INST_LIBS:=lib_mod.a
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:= 

include $(builddir)/make.rules
 
VPATH=$(topdir)/librlinux

.PHONY: all
all: $(INST_PROGS) $(INST_LIBS)

src+=$(wildcard $(topdir)/librlinux/*.f) 
src+=$(wildcard $(topdir)/librlinux/*.f90) 
src+=$(wildcard $(topdir)/librlinux/*.c)
objs:=$(addsuffix .o, $(basename $(src)))
objs:=$(notdir $(objs))

special:=expnd4.o  # pgf77 has problems... see below.

$(libdir)/lib_mod.a: $(objs) 
	$(AR) rv $@ $?
	$(RAN) $@

ifeq ($(FC),pgf77)
$(special): %.o: %.f
	$(FC) -O2 -pc 64 -tp p6  -c $<
endif

