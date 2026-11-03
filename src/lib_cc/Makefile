#
# Written by Jonas Juselius <jonas@iki.fi> Tue Dec  4 11:41:04 EET 2001
#
#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
builddir:=..
#moddir:=../lib_mod
moddir:=.
include $(builddir)/make.config

INST_PROGS:=
INST_LIBS:=lib_cc.a
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=mod_mult mod_contract mod_common  mod_operator mod_davidson

include $(builddir)/make.rules
 
.PHONY: all
all: $(INST_PROGS) $(INST_LIBS)
objs:=$(addsuffix .o, $(basename $(src)))
objs:=$(notdir $(objs))

#$(shell cp ../lib_mod/mod* .)

special:=expnd4.o  # pgf77 has problems... see below.

#libs:=-lsrint -l_mod
$(libdir)/lib_cc.a: $(objs)
	$(AR) rv $@ $?
	$(RAN) $@ 

# ==============================================
# 强制模块优先编译：所有目标文件都依赖模块先生成
# ==============================================
MODULE_SRC := mod_davidson.f90 mod_mult.f90 mod_contract.f90 mod_common.f90  mod_operator.f90
MODULE_OBJ := mod_davidson.o mod_mult.o mod_contract.o mod_common.o  mod_operator.o

# 让所有 .o 文件都必须等模块编译完成才开始编译
$(filter-out $(MODULE_OBJ),$(objs)): $(MODULE_OBJ)
mod_contract.o mod_mult.o mod_common.o: mod_operator.o
mod_davidson.o: mod_common.o

# 额外确保模块本身先被编译
all: $(MODULE_OBJ)
