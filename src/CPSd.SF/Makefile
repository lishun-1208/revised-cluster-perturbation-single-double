#
# Written by Jonas Juselius <jonas@iki.fi> Tue Dec  4 11:41:04 EET 2001
#
#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
builddir:=..
installdir:=../..
 
include $(builddir)/make.config

INST_PROGS:=xCPSd.SF
INST_LIBS:=
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:= mod_variable

include $(builddir)/make.rules
 
all: $(INST_PROGS) $(INST_LIBS)

# l_cc must before lsrint and l_mod
libs:=-l_cc -lsrint -l_mod

$(bindir)/xCPSd.SF: $(objs) $(libs)
	$(FC) $(FLDFLAGS) -o $@ $(objs) $(libs) $(BLAS_LIBS)

# ==============================================
# 强制模块优先编译：所有目标文件都依赖模块先生成
# ==============================================
MODULE_SRC := mod_variable.f90
MODULE_OBJ := mod_variable.o

# 让所有 .o 文件都必须等模块编译完成才开始编译
$(filter-out $(MODULE_OBJ),$(objs)): $(MODULE_OBJ)

# 额外确保模块本身先被编译
all: $(MODULE_OBJ)

