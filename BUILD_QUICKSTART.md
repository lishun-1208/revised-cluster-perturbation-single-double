# Building CPSd on another machine

This short guide covers the seven programs and the `pyscf2cc` Python extensions, excluding `cfour2cc`. Use the current sources and Makefiles together. The tested platform is Ubuntu 22.04 x86_64; other Linux systems may need different library and header paths. Detailed installation commands are available in [BUILD_DEPENDENCIES.md](BUILD_DEPENDENCIES.md).

## 1. Dependencies

Install these dependencies or load equivalent compiler/library modules on your machine:

| Dependency | Tested version / requirement |
| --- | --- |
| Build tools and system development files | GCC/G++ 11, GNU Make, binutils, C/C++ headers, glibc development files, and OpenMP runtime. |
| NVIDIA HPC SDK | 24.5, including `nvfortran`, CUDA 12.4, and cuBLAS. Required by the main programs even when running the CPU example. |
| GNU Fortran | GFortran 11, including its development files, `libgfortran`, and `libquadmath`; used for the Python extensions. |
| OpenBLAS with LAPACK | 0.3.20, **ILP64**, with unsuffixed Fortran symbols such as `dgemm_`. Ordinary LP64 OpenBLAS is incompatible. |
| Intel oneMKL | 2024.2.2 runtime; ILP64 interface, sequential layer, core library, and VML kernels. Needed for `vdmul`. |
| Python | Python 3.10 and matching development headers. |
| Python packages | NumPy 1.26.4, setuptools 69.5.1, wheel 0.48.0, packaging 26.3; SciPy 1.15.3, h5py 3.16.0, and PySCF 2.14.0 for integral generation. |

A virtual environment is recommended for Python packages. Keep Python 3.10 and NumPy 1.26 for the existing f2py build. The installed dependencies can live under `~/lish/soft`, system directories, or a cluster software directory, provided the paths below match.

## 2. Adjust the build configuration

Edit **`src/make.config`** to match the new machine:

| Variable | What to set |
| --- | --- |
| `SOFT_ROOT` | Base directory for locally installed dependencies; default: `$(HOME)/lish/soft`. |
| `NVHPC_ROOT` | SDK version directory containing `compilers/bin/nvfortran`; default: `$(SOFT_ROOT)/nvidia/hpc_sdk/Linux_x86_64/24.5`. |
| `DEPS_ROOT` | Dependency prefix containing `bin`, `include`, and `lib`; default: `$(SOFT_ROOT)/cpsd-deps/usr`. |
| `DEPS_LIB` | Directory containing `libgfortran.so.5` and `libquadmath.so.0`. The current Debian multiarch path may differ on your machine. |
| `OPENBLAS_LIB` | Directory containing `libopenblas64.so`. Adapt the library filename in `BLAS_LIBS` if your ILP64 installation uses another name. |
| `MKL_LIB` | Directory containing `libmkl_intel_ilp64.so.2`, `libmkl_sequential.so.2`, `libmkl_core.so.2`, and the other MKL runtime libraries. Adapt filenames in `VML_LIBS` if needed. |
| `PYTHON` | Absolute path to the Python 3.10 interpreter, preferably in a virtual environment. |
| `GNU_FC` | Absolute path to GFortran or its local wrapper; the main `FC` remains NVIDIA `nvfortran`. |
| `GPUFLAGS` | CUDA version and GPU targets supported by your SDK. Current setting: `-acc -cuda -cudalib=cublas -gpu=cc70,cc80,cc90,cuda12.4`. |

Also check these two files:

- **`src/pyscf2cc/Makefile`**: adjust `CFLAGS` if Python headers are outside `$(DEPS_ROOT)/include/python3.10`. Include the directory containing `Python.h` and any directory needed for the architecture-specific `pyconfig.h`.
- **`build-env.sh`**: adjust the software prefix, SDK version, and library paths for runtime use. Its current prefix is `$HOME/lish/soft`.

Keep the main programs' `-i8 -DINT64` and the extensions' `-fdefault-integer-8` consistent with ILP64 BLAS. Keep OpenBLAS before MKL in the link line and retain the MKL `--push-state` / `--pop-state` group. Do not pass a global `FC=...` to make, since the programs and extensions use different Fortran compilers. If reusing a relocated NVIDIA SDK, regenerate its machine configuration with `makelocalrc` as described in the detailed guide.

## 3. Compile and check

After installing dependencies and updating the paths, run in Bash:

```bash
cd /path/to/CPSd
source ./build-env.sh

# Remove objects from another machine or an earlier compiler configuration.
make -C src clean
make -C src -j8
```

The top-level build creates its output directories and builds internal libraries before programs. It skips `src/cfour2cc`. Use `clean` when changing compilers, library paths, or GPU flags; avoid `distclean`, which removes configuration files.

The executables appear in `src/bin/`:

```text
xgetHF          xCPSd.SF         xEE-CIS.SF       xEE-CPSd.SF
xgrCPSd.SF      xCISdoo.SD.SF    xEE-grCPSd.SF
```

The two Python extensions appear in `src/pyscf2cc/`. Check their imports and the executable's runtime libraries:

```bash
python -c 'import numpy, pyscf, pyscf2cc.write_int, pyscf2cc.interface; print("Imports OK")'
ldd src/bin/xgetHF
```

Resolve any `not found` libraries by correcting `build-env.sh` and sourcing it again. Keep `src/bin` on `PATH` and the repository's `src` directory on `PYTHONPATH`. For a numerical check, follow the H2O/STO-3G smoke test in [the detailed guide](BUILD_DEPENDENCIES.md#6-set-the-runtime-environment-and-verify). CPU execution passed on the original machine; GPU execution still needs a suitable GPU and working driver.
