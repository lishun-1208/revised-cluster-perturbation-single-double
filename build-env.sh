# Source this file before building or running CPSd.
_cpsd_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
_cpsd_soft="$HOME/lish/soft"
export PATH="$_cpsd_root/src/bin:$_cpsd_soft/cpsd-python/bin:$_cpsd_soft/cpsd-deps/usr/bin:$_cpsd_soft/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/bin:$PATH"
export PYTHONPATH="$_cpsd_root/src${PYTHONPATH:+:$PYTHONPATH}"
export LD_LIBRARY_PATH="$_cpsd_soft/mkl/lib:$_cpsd_soft/cpsd-deps/usr/lib/x86_64-linux-gnu:$_cpsd_soft/cpsd-deps/usr/lib/x86_64-linux-gnu/openblas64-pthread:$_cpsd_soft/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/lib:$_cpsd_soft/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset _cpsd_root _cpsd_soft
