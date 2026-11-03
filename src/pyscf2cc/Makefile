topdir := ..
builddir := ..
include ../make.config

FC := $(GNU_FC)
GCC_SUPPORT := $(shell gcc -print-libgcc-file-name | xargs dirname)
FFLAGS := -B$(GCC_SUPPORT)/ -fdefault-integer-8 -fPIC -O2 -fopenmp
F2PY_F90FLAGS := $(FFLAGS) -I$(CURDIR)
EXT_SUFFIX := $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
SO_WRITE_INT := write_int$(EXT_SUFFIX)
SO_INTERFACE := interface$(EXT_SUFFIX)
# Python development headers are extracted locally, without system installation.
export CFLAGS := -I$(DEPS_ROOT)/include/python3.10 -I$(DEPS_ROOT)/include
export LDFLAGS := -B$(GCC_SUPPORT)/ -L$(DEPS_LIB) -Wl,-rpath,$(DEPS_LIB) -Wl,--disable-new-dtags
export FC F90 := $(FC)

.PHONY: all clean install
all: $(SO_WRITE_INT) $(SO_INTERFACE)

mod_py2cc.o: mod_py2cc.f90 Makefile ../make.config
	$(FC) $(FFLAGS) -c $< -o $@

mod_interface.o: mod_interface.f90 Makefile ../make.config
	$(FC) $(FFLAGS) -c $< -o $@

libmod.a: mod_py2cc.o mod_interface.o
	$(AR) rcs $@ $^

$(SO_WRITE_INT): libmod.a write_int.f90 Makefile ../make.config
	$(PYTHON) -m numpy.f2py -m write_int -c write_int.f90 \
	  --fcompiler=gnu95 --f90exec=$(FC) --f77exec=$(FC) \
	  --f90flags="$(F2PY_F90FLAGS)" -L$(CURDIR) -lmod -lgomp

$(SO_INTERFACE): libmod.a interface.f90 Makefile ../make.config
	$(PYTHON) -m numpy.f2py -m interface -c interface.f90 \
	  --fcompiler=gnu95 --f90exec=$(FC) --f77exec=$(FC) \
	  --f90flags="$(F2PY_F90FLAGS)" -L$(CURDIR) -lmod -lgomp

# This source directory is a Python namespace package (see build-env.sh).
install: all

clean:
	rm -f *.o *.mod *.so *.a
