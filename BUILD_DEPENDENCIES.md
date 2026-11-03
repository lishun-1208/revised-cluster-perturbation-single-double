# CPSd build dependencies and reproducible build guide

For a short overview of dependencies, configuration changes, and build commands, see [BUILD_QUICKSTART.md](BUILD_QUICKSTART.md).

This guide builds the seven programs and two Python extensions in this checkout. It intentionally excludes `src/cfour2cc`. It describes the successful build on **Ubuntu 22.04.5 LTS, x86_64**, tested on 2026-11-03. The checked-in Makefiles and the current `getHF` source fix are required; an older checkout may need the same build-file changes. The commands below keep locally added software under `~/lish/soft` and do not install Ubuntu packages system-wide, except for the administrator's base-tool step.

## What is required

| Component | Tested version | Origin and purpose |
| --- | --- | --- |
| Ubuntu, glibc | 22.04.5; glibc 2.35 | Existing system and runtime ABI. |
| GCC and G++ | 11.4.0 | Existing C/C++ compilers; GCC support files also serve local GFortran. |
| GNU Make, binutils | 4.3; 2.38 | Existing build driver, assembler, linker, `ar`, and `ranlib`. |
| C/C++ development files | `build-essential` 12.9ubuntu3; `libc6-dev` 2.35; `libgcc-11-dev` and `libstdc++-11-dev` 11.4.0 | Existing headers, start files, and GCC linker support. `libc6-dev` also provides `math-vector-fortran.h`. |
| Other system libraries | `libgomp1`, `libatomic1`, `libnuma1`, `libncursesw6`, `libtinfo6`, `zlib1g` | Existing compiler/runtime dependencies. |
| Python interpreter | CPython 3.10.12 | Existing `/usr/bin/python3`; used as the base of a local virtual environment. System `pip` and `python3-venv` were **not** present on the tested machine. |
| Shell and download tools | Bash 5.1, coreutils, findutils (`xargs`), `dpkg-deb`, `tar`, `gzip`, `wget`, CA certificates | Existing tools used by the commands below. |
| NVIDIA HPC SDK | 24.5-1, CUDA 12.4 / cuBLAS 12.4.2.65 | Installed locally; `nvfortran` builds the main programs, including CUDA Fortran code. |
| GFortran | 11.4.0 | Ubuntu packages extracted locally; builds the f2py extensions. |
| OpenBLAS with LAPACK | 0.3.20, pthread, **ILP64** | Ubuntu packages extracted locally; supplies BLAS/LAPACK with 64-bit integers and unsuffixed Fortran symbols. |
| GNU Fortran runtime | `libgfortran5` and `libquadmath0` 12.3.0 | Ubuntu packages extracted locally. |
| Python 3.10 headers | 3.10.12 | Ubuntu package extracted locally; needed by f2py. |
| Intel oneMKL runtime | `mkl` wheel 2024.2.2 | Installed locally; supplies the `vdmul` VML routine via ILP64, sequential, and core libraries. |
| Python packages | NumPy 1.26.4, SciPy 1.15.3, h5py 3.16.0, PySCF 2.14.0 | Local virtual environment. NumPy builds the extensions; SciPy, h5py, and PySCF support integral generation and Python runtime. |
| Python build tools | pip 26.3, setuptools 69.5.1, wheel 0.48.0, packaging 26.3 | Local virtual environment. |

The pre-existing Ubuntu package revisions on the test host were:

| Packages already installed | Observed Ubuntu release |
| --- | --- |
| `gcc-11`, `g++-11`, `libgcc-11-dev`, `libstdc++-11-dev` | `11.4.0-1ubuntu1~22.04.3` |
| `build-essential`; `make`; `binutils` | `12.9ubuntu3`; `4.3-4.1build1`; `2.38-4ubuntu2.12` |
| `libc6`, `libc6-dev` | `2.35-0ubuntu3.15` |
| `linux-libc-dev` | `5.15.0-202.212` |
| `libgcc-s1`, `libstdc++6`, `libgomp1`, `libatomic1` | `12.3.0-1ubuntu1~22.04.3` |
| `libnuma1` | `2.0.14-3ubuntu2` |
| `libncursesw6`, `libtinfo6` | `6.3-2ubuntu0.3` |
| `zlib1g` | `1:1.2.11.dfsg-2ubuntu9.2` |
| `python3.10`, `python3.10-minimal`, `libpython3.10-stdlib` | `3.10.12-1~22.04.18` |
| `bash`; `coreutils`; `findutils`; `dpkg` | `5.1-6ubuntu1.1`; `8.32-4.1ubuntu1.4`; `4.8.0-1ubuntu3`; `1.21.1ubuntu2.6` |
| `tar`; `gzip`; `wget`; `ca-certificates` | `1.34+dfsg-1ubuntu0.1.22.04.6`; `1.10-4ubuntu4.2`; `1.21.2-2ubuntu1.5`; `20260816ubuntu0.22.04.1` |

Versions above are a **tested combination**, not a guarantee for every Linux distribution. A working NVIDIA GPU driver was unavailable on the test host; CPU execution of the example passed. The CUDA compiler and cuBLAS are still required by the sources and build flags.

The main Fortran build uses `-i8 -DINT64`; the Python extensions use GFortran `-fdefault-integer-8`. Use the ILP64 OpenBLAS library with normal names such as `dgemm_`. An LP64 OpenBLAS library, or an ILP64 build that exports only suffixed names such as `dgemm_64_`, is not interchangeable. OpenBLAS appears first in the link line for BLAS/LAPACK, while oneMKL's sequential VML layer supplies `vdmul`. [Intel documents the oneMKL library layers](https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2023-0/dynamic-libraries-in-the-lib-intel64-lin-directory.html). The build uses NVIDIA's OpenMP runtime for the programs and GNU `libgomp` for Python extensions; sequential oneMKL avoids another OpenMP runtime.

## 1. Prepare the system and paths

On a fresh Ubuntu 22.04 x86_64 machine, an administrator can install the base tools once. The Python venv package is convenient for a new machine, although it was absent from the original test host and is not needed by the `--without-pip` method below.

```bash
sudo apt update
sudo apt install build-essential gcc-11 g++-11 binutils libc6-dev linux-libc-dev \
  libgcc-11-dev libstdc++-11-dev libgomp1 libatomic1 libnuma1 \
  libncursesw6 libtinfo6 zlib1g python3.10 python3.10-venv \
  ca-certificates wget tar gzip dpkg file coreutils findutils grep sed mawk
```

Then, as an ordinary user, start in the repository root and keep these variables in the **same Bash session** for the remaining commands. Change `SOFT_ROOT` if desired. Paths with spaces are not supported by the current Makefiles.

```bash
cd /path/to/CPSd
export CPSD_ROOT="$(pwd -P)"
export SOFT_ROOT="$HOME/lish/soft"
export DEPS_ROOT="$SOFT_ROOT/cpsd-deps/usr"
export DEPS_LIB="$DEPS_ROOT/lib/x86_64-linux-gnu"
export OPENBLAS_LIB="$DEPS_LIB/openblas64-pthread"
export MKL_LIB="$SOFT_ROOT/mkl/lib"
export NVHPC_ROOT="$SOFT_ROOT/nvidia/hpc_sdk/Linux_x86_64/24.5"
export CPSD_PYTHON="$SOFT_ROOT/cpsd-python/bin/python"
export CPSD_DOWNLOADS="$SOFT_ROOT/cpsd-downloads"
mkdir -p "$CPSD_DOWNLOADS" "$SOFT_ROOT/cpsd-deps"
gcc --version | head -1
g++ --version | head -1
python3.10 --version
```

On Ubuntu 22.04, `gcc` and `g++` should resolve to GCC 11 for this recipe. If the first two commands show another major release, arrange for GCC/G++ 11 to be first on `PATH` before installing the SDK or building; the local GFortran wrapper discovers GCC support files from the active `gcc`.

The recipes assume the Ubuntu Jammy repositories are available to `apt-get download`, PyPI and the NVIDIA download server are reachable, and at least **40–45 GiB** of temporary and installed space is available for the large SDK archive, extraction, and final install. The tested archive was about 6.77 GB and the installed SDK about 14 GiB. For offline reproducibility, retain the `.deb` files, SDK archive, and Python wheels in a separate artifact store.

## 2. Extract local Ubuntu packages

`apt-get download` fetches packages from the configured Ubuntu 22.04 repositories; `dpkg-deb -x` extracts them into the user prefix. Nothing in this step calls `dpkg -i` or changes the system package database. The seven tested package releases were:

| Package | Tested Ubuntu release |
| --- | --- |
| `gfortran-11`, `libgfortran-11-dev` | `11.4.0-1ubuntu1~22.04.3` |
| `libgfortran5`, `libquadmath0` | `12.3.0-1ubuntu1~22.04.3` |
| `libopenblas64-0-pthread`, `libopenblas64-pthread-dev` | `0.3.20+ds-1` |
| `libpython3.10-dev` | `3.10.12-1~22.04.18` |

Jammy security updates may provide later compatible revisions. Save the downloaded archives if you need to reproduce these exact revisions later.

```bash
DEB_WORK="$(mktemp -d)"
cd "$DEB_WORK"
apt-get download gfortran-11 libgfortran-11-dev libgfortran5 libquadmath0 \
  libopenblas64-pthread-dev libopenblas64-0-pthread libpython3.10-dev
for package in ./*.deb; do
  dpkg-deb -x "$package" "$SOFT_ROOT/cpsd-deps"
  cp "$package" "$CPSD_DOWNLOADS/"
done
```

The extracted GFortran driver needs the system GCC 11 linker-support directory and the local runtime directory. Create a wrapper; `install` replaces the extracted `gfortran` symlink without writing through it. This wrapper computes the prefix at runtime, so the recipe does not contain another user's home directory.

```bash
cat > "$CPSD_DOWNLOADS/gfortran-wrapper" <<'SH'
#!/bin/sh
cpsd_deps_usr=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cpsd_gcc_support=$(gcc -print-libgcc-file-name)
cpsd_gcc_support=${cpsd_gcc_support%/*}
exec "$cpsd_deps_usr/bin/gfortran-11" \
  "-B$cpsd_gcc_support/" \
  "-L$cpsd_deps_usr/lib/x86_64-linux-gnu" \
  "-Wl,-rpath,$cpsd_deps_usr/lib/x86_64-linux-gnu" "$@"
SH
install -m 755 "$CPSD_DOWNLOADS/gfortran-wrapper" "$DEPS_ROOT/bin/gfortran"
export PATH="$DEPS_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$DEPS_LIB:$OPENBLAS_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$DEPS_ROOT/bin/gfortran" --version | head -1
test -f "$DEPS_ROOT/include/python3.10/Python.h"
test -f "$DEPS_ROOT/include/x86_64-linux-gnu/python3.10/pyconfig.h"
test -f "$OPENBLAS_LIB/libopenblas64.so"
```

The `libpython3.10-dev` headers need **both** `$DEPS_ROOT/include/python3.10` and `$DEPS_ROOT/include`: the first provides `Python.h`, and the second resolves the architecture-specific `pyconfig.h`. The repository's `src/pyscf2cc/Makefile` already supplies both. Its include path is specific to Python 3.10.

## 3. Install NVIDIA HPC SDK 24.5

Use the [NVIDIA HPC SDK 24.5 installation guide](https://docs.nvidia.com/hpc-sdk/archive/24.5/hpc-sdk-install-guide/index.html) and [release notes](https://docs.nvidia.com/hpc-sdk/archive/24.5/hpc-sdk-release-notes/index.html) for platform and installer details. A functioning GCC, G++, and GFortran setup must be visible before running the installer. This repository requires `nvfortran` and its CUDA/cuBLAS libraries even for the CPU smoke test because it contains CUDA Fortran declarations and cuBLAS calls.

```bash
cd "$CPSD_DOWNLOADS"
wget -c -O nvhpc_2024_245_Linux_x86_64_cuda_12.4.tar.gz \
  https://developer.download.nvidia.com/hpc-sdk/24.5/nvhpc_2024_245_Linux_x86_64_cuda_12.4.tar.gz
sha256sum nvhpc_2024_245_Linux_x86_64_cuda_12.4.tar.gz
tar -xzf nvhpc_2024_245_Linux_x86_64_cuda_12.4.tar.gz
cd nvhpc_2024_245_Linux_x86_64_cuda_12.4
NVHPC_SILENT=true NVHPC_INSTALL_TYPE=single NVHPC_DEFAULT_CUDA=12.4 \
  NVHPC_INSTALL_DIR="$SOFT_ROOT/nvidia/hpc_sdk" ./install
"$NVHPC_ROOT/compilers/bin/makelocalrc" \
  -x "$NVHPC_ROOT/compilers/bin" \
  -gcc /usr/bin/gcc -gpp /usr/bin/g++ \
  -g77 "$DEPS_ROOT/bin/gfortran" -cuda 12.4
"$NVHPC_ROOT/compilers/bin/nvfortran" -V
```

The observed SHA-256 of the **tested local archive** was `392c3315dfe4c8bbd6d77ac9efff7a1af837447cca564cfbeff2d10682d365f5`. This is a local reproducibility record, not a vendor-published checksum. `makelocalrc` writes machine-specific compiler configuration; rerun it when moving the installation or changing the GNU toolchain. A missing `liblto_plugin.so` during this step usually means the GFortran wrapper or GCC 11 support files are unavailable.

## 4. Create the Python environment and install oneMKL

The tested environment used system Python 3.10 without system pip. Bootstrapping a local venv works whether or not `python3.10-venv` and pip were initially installed system-wide. The versions below match [`build-logs/python-packages.txt`](build-logs/python-packages.txt).

```bash
python3.10 -m venv --without-pip "$SOFT_ROOT/cpsd-python"
wget -O "$CPSD_DOWNLOADS/get-pip.py" https://bootstrap.pypa.io/get-pip.py
"$CPSD_PYTHON" "$CPSD_DOWNLOADS/get-pip.py" --no-cache-dir
"$CPSD_PYTHON" -m pip install --no-cache-dir \
  pip==26.3 numpy==1.26.4 setuptools==69.5.1 wheel==0.48.0 \
  packaging==26.3 scipy==1.15.3 h5py==3.16.0 pyscf==2.14.0
"$CPSD_PYTHON" -m pip freeze > "$CPSD_DOWNLOADS/python-packages.txt"
"$CPSD_PYTHON" -m pip install --no-cache-dir --no-deps \
  --prefix "$SOFT_ROOT/mkl" mkl==2024.2.2
test -f "$MKL_LIB/libmkl_intel_ilp64.so.2"
test -f "$MKL_LIB/libmkl_sequential.so.2"
test -f "$MKL_LIB/libmkl_core.so.2"
```

Check the OpenBLAS integer interface before building:

```bash
"$CPSD_PYTHON" - <<'PY'
import ctypes
import os

library = ctypes.CDLL(os.path.join(os.environ['OPENBLAS_LIB'], 'libopenblas64.so'))
library.openblas_get_config.restype = ctypes.c_char_p
config = library.openblas_get_config().decode()
print(config)
assert 'USE64BITINT' in config, 'OpenBLAS must use 64-bit integers (ILP64)'
PY
```

The `mkl` wheel has its own prefix and is **not** included in the virtual environment's `pip freeze`. `--no-deps` applies only to that MKL installation: the build uses its sequential libraries and does not need the wheel's optional threading dependencies. Keep the entire `mkl/lib` directory because VML may load additional dispatched kernel libraries. The Python packages need normal dependency resolution.

NumPy 1.26 retains `numpy.distutils`, which this f2py build uses. [NumPy documents its removal for Python 3.12](https://numpy.org/doc/1.26/reference/distutils_status_migration.html). The tested setuptools 69.5.1 worked in this checkout; for a reproducible build, keep the tested versions instead of assuming that a newer Python/NumPy/setuptools stack is interchangeable. Python development headers must match the 3.10 interpreter.

## 5. Build

Return to the repository. The local prefix can be overridden at make time. Specify the GFortran wrapper through `GNU_FC`; **do not pass `FC=...` globally**, since the main programs and f2py extensions intentionally use different Fortran compilers.

```bash
cd "$CPSD_ROOT"
make -C src -j8 \
  SOFT_ROOT="$SOFT_ROOT" DEPS_ROOT="$DEPS_ROOT" \
  NVHPC_ROOT="$NVHPC_ROOT" MKL_LIB="$MKL_LIB" \
  PYTHON="$CPSD_PYTHON" GNU_FC="$DEPS_ROOT/bin/gfortran"
```

Expected products:

| Location | Files |
| --- | --- |
| `src/bin/` | `xgetHF`, `xCPSd.SF`, `xEE-CIS.SF`, `xEE-CPSd.SF`, `xgrCPSd.SF`, `xCISdoo.SD.SF`, `xEE-grCPSd.SF` |
| `src/lib/` | `lib_mod.a`, `libsrint.a`, `lib_cc.a` |
| `src/pyscf2cc/` | `write_int*.so`, `interface*.so` with Python 3.10 extension suffixes |

The top-level `src/Makefile` builds the internal libraries in dependency order, then all seven programs and the two extensions. It does not enter `src/cfour2cc`. No CFOUR or HDF5 Fortran installation is needed for this build; Python HDF5 support comes from h5py. The build uses `-gpu=cc70,cc80,cc90,cuda12.4`. If you change `GPUFLAGS`, toolchain, ABI, or command-line prefix variables after a build, first run `make -C "$CPSD_ROOT/src" clean` and then rerun the build with the same overrides. Use `clean`, not `distclean`; this is an in-tree build, and `distclean` removes generated configuration files.

## 6. Set the runtime environment and verify

For the default prefix (`$HOME/lish/soft`), source the repository's existing script:

```bash
source "$CPSD_ROOT/build-env.sh"
command -v python
command -v xgetHF
python -c 'import numpy, pyscf, pyscf2cc.write_int, pyscf2cc.interface; print("Python imports: OK")'
```

For a custom `SOFT_ROOT`, use these exports instead. `build-env.sh` currently selects `$HOME/lish/soft` by design.

```bash
export PATH="$CPSD_ROOT/src/bin:$SOFT_ROOT/cpsd-python/bin:$DEPS_ROOT/bin:$NVHPC_ROOT/compilers/bin:$PATH"
export PYTHONPATH="$CPSD_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export LD_LIBRARY_PATH="$MKL_LIB:$DEPS_LIB:$OPENBLAS_LIB:$NVHPC_ROOT/compilers/lib:$SOFT_ROOT/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

The `pyscf2cc` modules import through `PYTHONPATH=$CPSD_ROOT/src`; there is no separate `pip install` step for this repository. Source the environment explicitly in noninteractive job scripts. If desired, add a guarded `source /path/to/CPSd/build-env.sh` to your shell startup file for the default prefix.

The following small H₂O/STO-3G example checks PySCF integral generation and then runs all seven programs in the order used for the successful CPU test. It uses its own temporary working directory so it does not overwrite repository examples.

```bash
(
set -e
SMOKE_DIR="$(mktemp -d)"
cd "$SMOKE_DIR"
cat > CCSOINP <<'EOF'
ref_state=rhf
shell=close
estate_sym=0/1/0/0/0/0/0/0
ex_mult=1
iostate=0
spin_free=t
memory=1
reserved_size=100000
EOF
cat > mol.py <<'PY'
from pyscf import gto, scf
from pyscf2cc import py2cc

mol = gto.M(
    atom='O 0 0 -0.13209669; H 0 1.43152878 0.97970006; H 0 -1.43152878 0.97970006',
    unit='B', basis='sto-3g', symmetry='C2v', verbose=3,
)
mf = scf.RHF(mol).run()
assert mf.converged
py2cc.py2cc(mf, mol, 'RHF')
print('PASS: H2O/STO-3G integral generation', mf.e_tot)
PY
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python mol.py > integrals.log 2>&1
for program in xgetHF xCPSd.SF xEE-CIS.SF xEE-CPSd.SF \
               xgrCPSd.SF xCISdoo.SD.SF xEE-grCPSd.SF; do
  "$program" > "$program.log" 2>&1 || { tail -50 "$program.log"; exit 1; }
  printf 'Completed %s; inspect %s.log\n' "$program" "$program"
done
printf 'Smoke-test logs: %s\n' "$SMOKE_DIR"
)
```

Inspect the end of every log for normal completion: a Fortran `STOP` can sometimes yield exit status zero even after an application error. On the tested checkout, all eight steps completed; the RHF energy was about `-74.9632607411` Hartree, CPS(D) fourth-order total energy about `-75.0115089187` Hartree, and grCPS(D) fourth-order total energy about `-75.0115985386` Hartree. Small numerical differences across machines are expected. The detailed [successful CPU-test summary](build-logs/retest2/summary.json) and logs are retained. This smoke test does not validate GPU execution or every scientific mode.

## Common failures and portability limits

| Symptom | Check |
| --- | --- |
| `liblto_plugin.so` missing during `makelocalrc` or linking | GCC 11, `libgcc-11-dev`, and the local GFortran wrapper; then rerun `makelocalrc`. |
| `Python.h` or `pyconfig.h` missing | Local `libpython3.10-dev` extraction and both include directories in `src/pyscf2cc/Makefile`; interpreter must be Python 3.10. |
| `numpy.distutils`/f2py fails | Activate the pinned Python 3.10/NumPy 1.26 stack and use the `GNU_FC` wrapper. Do not mix NVIDIA and GNU `.mod` files. |
| `vdmul_` unresolved | Install the complete oneMKL wheel and pass `MKL_LIB` to make. |
| BLAS/LAPACK linkage succeeds but results are wrong or crash | Check OpenBLAS is ILP64 with unsuffixed Fortran symbols; do not substitute LP64. |
| cuBLAS symbols unresolved | Check NVIDIA SDK/CUDA 12.4 paths and retain the Makefile's `--push-state,--no-as-needed ... --pop-state` MKL link group. |
| Shared library missing at runtime | Source `build-env.sh` or export the custom-prefix library paths; do not move a compiled prefix without rebuilding or fixing embedded runtime paths. |
| `getHF` rejects reference state | Use lowercase `ref_state=rhf` or `uhf` in `CCSOINP`. |

The provided Makefiles target Ubuntu 22.04 x86_64: they contain Debian multiarch library paths, Python **3.10** header paths, NVIDIA's `Linux_x86_64/24.5` layout, and an x86_64 oneMKL wheel. Other Linux distributions or architectures require corresponding Makefile and dependency-layout changes. NVIDIA's local compiler configuration and default CPU target are machine-specific, so build on the destination machine when possible. GPU execution additionally requires a supported NVIDIA GPU and driver; it was not tested here.
