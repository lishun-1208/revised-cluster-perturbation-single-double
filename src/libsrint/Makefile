#
# Written by Jonas Juselius <jonas@iki.fi> Tue Dec  4 11:41:04 EET 2001
#
#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
builddir:=..
moddir:=../lib_mod
include $(builddir)/make.config

INST_PROGS:=
INST_LIBS:=libsrint.a
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:= 

include $(builddir)/make.rules
 
.PHONY: all
all: $(INST_PROGS) $(INST_LIBS)
objs:=$(addsuffix .o, $(basename $(src)))
objs:=$(notdir $(objs))

#$(shell cp ../lib_mod/mod* .)

special:=expnd4.o  # pgf77 has problems... see below.


$(libdir)/libsrint.a: $(objs) 
	$(AR) rv $@ $? 
	$(RAN) $@ 


