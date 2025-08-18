#!/usr/bin/env bash

# Build Loci third‑party deps from ext/ and emit a helpful ./configure line.
# This builds all the necessary dependencies for Loci, and gives the user the option
# to either build using Scotch or ParMETIS.
#
# Mirrors Loci's configure options: --with-metis, --with-parmetis, --with-scotch.
#
# Notes:
# - If --with-parmetis is set, METIS and GKlib will also be built.
# - If only --with-metis is set, GKlib will be built automatically.
# - If --with-scotch is set, both METIS, ParMETIS, and GKlib are disabled.
# - Build directories are ext/<dep>_build, installs are ext/<dep>_install by default.
set -euo pipefail

# --- globals ---
deps="hdf5,cgns,scotch,petsc" # default dependencies
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # script root directory
ext="$root/ext" # external dependencies directory
jobs=$(command -v nproc >/dev/null 2>&1 && nproc || echo 4) # default parallel jobs
clean=0 # for removing build/install directories and starting fresh
clean_ext=0 # for cleaning out the git submodule directories completely (repos themselves)
with_scotch=1
with_metis=0
with_parmetis=0
prefix="" # installation prefix

# --- helpers ---
say() { echo -e "\033[1;34m[build-deps]\033[0m $*"; }
warn(){ echo -e "\033[1;33m[warn]\033[0m $*" >&2; }
die() { echo -e "\033[1;31m[error]\033[0m $*" >&2; exit 1; }
run_cmd() { echo -e "\033[1;36m$*\033[0m"; "$@"; }
has_dep() { [[ ",$deps," == *",$1,"* ]]; } # membership test for comma list like "hdf5,cgns,..."
need() { command -v "$1" >/dev/null 2>&1 || die "Missing tool: $1"; }
abspath(){ python3 -c 'import os,sys;print(os.path.abspath(sys.argv[1]))' "$1"; }
norm() { echo "$1" | tr 'A-Z' 'a-z' | tr -s ',' ',' | sed -e 's/^,//;s/,$//'; } # normalize comma-separated list

need git
need make
need gcc
need g++
need gfortran || warn "gfortran not found; some deps may fail if Fortran is required."
need cmake


# Check for the submodule
have_submodule() { [[ -d "$ext/$1" && -n "$(ls -A "$ext/$1" 2>/dev/null || true)" ]]; }
ensure_submodule() {
  local name="$1"
  if ! have_submodule "$name"; then
    say "Initializing submodule: ext/$name"
    (cd "$root" && git submodule update --init "ext/$name")
  fi
}

print_help() {
  cat <<EOF
Usage: $(basename "$0") [options]

Builds Loci's external dependencies from submodules.
Dependencies are built in ext/<dep>_build and installed to <prefix>/<dep>_install.

Options:
  --prefix        Path where all dependencies will be installed (default: ext/)
  --with-metis    Build with METIS support
  --with-parmetis Build with ParMETIS support
  --with-scotch   Build with Scotch support
  --jobs N        Number of parallel build jobs (default: nproc or 4)
  --clean         Remove <dep>_build and <dep>_install directories (fresh builds)
  --clean-ext     ***DESTRUCTIVE*** Wipe ALL ext/*_build, deinit all
                  submodules under ext/, and remove their working trees.
  -h, --help      Show this help and exit

Examples:
  $(basename "$0")   # Builds all default dependencies(Scotch instead of ParMETIS)
  $(basename "$0") --with-parmetis --jobs 8 # Build with ParMETIS and METIS support, using 8 jobs
  $(basename "$0") --with-scotch # Build with Scotch support

EOF
}


# Deep clean of ext/: wipe *_build, *_install, and leave ext/* as empty dirs (pre-init state)
deep_clean_ext() {
  say "==== Cleaning ext/ (builds, installs, and submodule trees) ===="

  # 1) Remove build/install dirs
  say "Removing *_build from $ext and *_install from $prefix"
  [[ -d "$ext"    ]] && find "$ext"    -maxdepth 1 -type d -name '*_build'   -print -exec rm -rf {} +
  [[ -d "$prefix" ]] && find "$prefix" -maxdepth 1 -type d -name '*_install' -print -exec rm -rf {} +

  # 2) Deinitialize submodules under ext/, but only if registered; then ensure empty dirs
  if git -C "$root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    say "Deinitializing submodules under ext/ (conditionally)"
    # enumerate submodule paths under ext/ from .gitmodules
    mapfile -t SUBS < <(git -C "$root" config --file "$root/.gitmodules" \
                         --get-regexp '^submodule\..*\.path$' \
                       | awk '{print $2}' | grep '^ext/')

    gitdir="$(git -C "$root" rev-parse --git-dir)"

    for path in "${SUBS[@]}"; do
      # derive the submodule "name" used in .git/config from .gitmodules
      name="$(git -C "$root" config --file "$root/.gitmodules" \
               --get-regexp "^submodule\\..*\\.path$" \
             | awk -v p="$path" '$2==p{print $1}' \
             | sed -E 's/^submodule\.|\.[^.]+$//g')"

      # is this submodule registered in .git/config?
      if git -C "$root" config --get "submodule.$name.url" >/dev/null 2>&1; then
        # if the per-module config exists, deinit; otherwise it was never initialized—skip
        if [[ -d "$gitdir/modules/$path" ]] || [[ -d "$gitdir/modules/${path#ext/}" ]]; then
          git -C "$root" submodule deinit -f --quiet -- "$path" || true
        fi
      fi

      # Ensure the directory exists, then empty it (including dotfiles) without removing the dir
      mkdir -p "$root/$path"
      find "$root/$path" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} + 2>/dev/null || true

      say "  - cleaned $path"
    done

    # Optional: prune any leftover cached module clones under .git/modules/ext
    rm -rf "$gitdir/modules/ext" 2>/dev/null || true
  else
    warn "Not a git repo; skipped submodule deinit. Only wiped *_build/_install."
  fi

  say "ext/ cleanup complete."
}



# ---------- HDF5 ----------
build_hdf5() {
  ensure_submodule hdf5
  say "==== HDF5 ===="
  local src="$ext/hdf5"
  local bld="$ext/hdf5_build"
  local inst="$prefix/hdf5_install"
  (( clean )) && rm -rf "$bld" "$inst"; mkdir -p "$bld" "$inst"

  run_cmd cmake -S "$src" -B "$bld" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_TESTING=OFF \
    -DHDF5_BUILD_TOOLS=ON \
    -DCMAKE_INSTALL_PREFIX="$inst"

  run_cmd cmake --build "$bld" --parallel "$jobs"
  run_cmd cmake --install "$bld"
}

# ---------- CGNS ----------
build_cgns() {
  ensure_submodule cgns
  say "==== CGNS ===="

  local src="$ext/cgns"
  local bld="$ext/cgns_build"
  local inst="$prefix/cgns_install"
  local h5="$prefix/hdf5_install"

  [[ -d "$h5" ]] || die "CGNS needs HDF5. Build HDF5 first (not found: $h5)."
  (( clean )) && rm -rf "$bld" "$inst"; mkdir -p "$bld" "$inst"

  run_cmd cmake -S "$src" -B "$bld" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$inst" \
    -DHDF5_ROOT="$h5" \
    -DCGNS_ENABLE_PARALLEL=OFF \
    -DHDF5_NEED_MPI=OFF

  run_cmd cmake --build "$bld" --parallel "$jobs"
  run_cmd cmake --install "$bld"
}

# ---------- Scotch / PT-Scotch ----------
build_scotch() {
  ensure_submodule scotch
  say "==== Scotch / PT-Scotch ===="

  # Require MPI for PT-Scotch
  local cc_mpi=
  if command -v mpicc >/dev/null 2>&1; then
    cc_mpi="mpicc"
  fi
  if [[ -z "$cc_mpi" ]]; then
    die "MPI compiler not found (mpicc). Install MPI or build with a non-PTScotch path."
  fi

  local src="$ext/scotch"
  local bld="$ext/scotch_build"
  local inst="$prefix/scotch_install"
  (( clean )) && rm -rf "$bld" "$inst"; mkdir -p "$bld" "$inst"

  # Configure: PT-Scotch ON, shared libs, use MPI wrappers
  run_cmd cmake -S "$src" -B "$bld" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$inst" \
    -DCMAKE_C_COMPILER="$cc_mpi" \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_PTSCOTCH=ON \
    -DSCOTCH_DETERMINISTIC=FULL \
    -DTHREADS=OFF \
    -DCMAKE_INSTALL_RPATH="$inst/lib" \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON

  run_cmd cmake --build "$bld" --parallel "$jobs"
  run_cmd cmake --install "$bld"

}

# ---------- PETSc ----------
build_petsc() {
  ensure_submodule petsc
  say "==== PETSc ===="

  local src="$ext/petsc"            # PETSC_DIR (source)
  local inst="$prefix/petsc_install"   # install prefix OUTSIDE submodule
  (( clean )) && rm -rf "$inst"; mkdir -p "$inst"
  (
    cd "$src"
    run_cmd ./configure --prefix="$inst" --with-debugging=0 --with-clanguage=C --with-precision=double --with-scalar-type=real COPTFLAGS="-O3" --with-x=0 --with-shared-libraries=1
    run_cmd make PETSC_DIR=. PETSC_ARCH=arch-linux-c-opt all -j"$jobs"
    run_cmd make PETSC_DIR=. PETSC_ARCH=arch-linux-c-opt install
  )
}

# ---------- GKlib (Needed only for METIS/ParMETIS) ----------
build_gklib() {
  ensure_submodule GKlib
  say "==== GKlib ===="

  local src="$ext/GKlib"
  local bld="$ext/GKlib_build"
  local inst="$prefix/GKlib_install"
  (( clean )) && rm -rf "$bld" "$inst"; mkdir -p "$bld" "$inst"

  run_cmd cmake -S "$src" -B "$bld" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DCMAKE_INSTALL_PREFIX="$inst"

  run_cmd cmake --build "$bld" --parallel "$jobs"
  run_cmd cmake --install "$bld"
}

# ---------- METIS ----------
build_metis() {
  ensure_submodule METIS
  say "==== METIS ===="

  local gk_inst="$ext/GKlib_install"
  [[ -d "$gk_inst/include" && -d "$gk_inst/lib" ]] || die "METIS needs GKlib installed."

  local src="$ext/METIS"
  local bld="$src/build"
  local inst="$prefix/METIS_install"
  (( clean )) && rm -rf "$bld" "$inst"; mkdir -p "$bld" "$inst"

  (
    cd "$src"
    run_cmd make config prefix="$inst" gklib_path="$gk_inst" shared=1
    run_cmd make -C build -j"$jobs"
    run_cmd make -C build install

  )
}

# ---------- ParMETIS ----------
build_parmetis() {
  ensure_submodule ParMETIS
  say "==== ParMETIS ===="

  # Require MPI compilers (ParMETIS is MPI)
  need mpicc
  need mpicxx

  local metis_inst="$prefix/METIS_install"
  local gklib_inst="$prefix/GKlib_install"
  [[ -d "$metis_inst" ]] || die "ParMETIS: METIS not installed at $metis_inst"
  [[ -d "$gklib_inst" ]] || die "ParMETIS: GKlib not installed at $gklib_inst"

  local src="$ext/ParMETIS"
  local bld="$ext/ParMETIS_build"
  local inst="$prefix/ParMETIS_install"
  (( clean )) && rm -rf "$bld" "$inst"; mkdir -p "$bld" "$inst"

  run_cmd cmake -S "$src" -B "$bld" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$inst" \
    -DSHARED=ON \
    -DMETIS_PATH="$metis_inst" \
    -DGKLIB_PATH="$gklib_inst" \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_INSTALL_RPATH="$metis_inst/lib:$gklib_inst/lib"

  run_cmd cmake --build "$bld" --parallel "$jobs"
  run_cmd cmake --install "$bld"
}

# ---------- driver ----------

# Parse command line args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      prefix="$(abspath "$2")"; shift 2 ;;
    --with-scotch)
      with_scotch=1; shift ;;
    --with-metis)
      with_metis=1; with_scotch=0; shift ;;
    --with-parmetis)
      with_parmetis=1; with_metis=1; with_scotch=0; shift ;;
    --jobs)
      [[ $# -ge 2 ]] || { echo "Missing value for --jobs" >&2; print_help; exit 2; }
      jobs="$2"; shift 2 ;;
    --clean)
      clean=1; shift ;;
    --clean-ext)
      clean_ext=1; shift ;;
    -h|--help)
      print_help; exit 0 ;;
    *)
      echo "Unknown arg: $1" >&2
      print_help
      exit 2 ;;
  esac
done

if [[ -z "$prefix" ]]; then
  prefix="$ext"   # default ext/ folder
fi

if (( clean_ext )); then
  say "Cleaning external dependencies..."
  deep_clean_ext
  echo "ext/ wiped. Exiting."
  exit 0
fi

# Apply boolean feature toggles to the build plan:
if (( with_parmetis )); then
  # ensure gklib+metis+parmetis present; remove scotch
  deps="gklib,metis,parmetis,${deps//scotch/}"
fi
if (( with_metis )); then
  # ensure gklib+metis present; remove scotch
  deps="gklib,metis,${deps//scotch/}"
fi

echo "Final dependencies: $deps"


for d in ${deps//,/ }; do
  case "$d" in
    hdf5)     build_hdf5 ;;
    cgns)     build_cgns ;;
    scotch)   build_scotch ;;
    petsc)    build_petsc ;;
    gklib)    build_gklib ;;
    metis)    build_metis ;;
    parmetis) build_parmetis ;;
    *)        die "Unknown dep: $d" ;;
  esac
done



# Build absolute install paths for deps actually requested
cfg_prefix="$(abspath "$root/loci_install")"

h5_inst="";  cgns_inst="";  petsc_inst="";  scotch_inst=""; metis_inst=""; parmetis_inst=""
if has_dep hdf5     && [[ -d "$prefix/hdf5_install"     ]]; then h5_inst="$(abspath "$prefix/hdf5_install")"; fi
if has_dep cgns     && [[ -d "$prefix/cgns_install"     ]]; then cgns_inst="$(abspath "$prefix/cgns_install")"; fi
if has_dep petsc    && [[ -d "$prefix/petsc_install"    ]]; then petsc_inst="$(abspath "$prefix/petsc_install")"; fi
if has_dep scotch   && [[ -d "$prefix/scotch_install"   ]]; then scotch_inst="$(abspath "$prefix/scotch_install")"; fi
if has_dep metis    && [[ -d "$prefix/METIS_install"    ]]; then metis_inst="$(abspath "$prefix/METIS_install")"; fi
if has_dep parmetis && [[ -d "$prefix/ParMETIS_install" ]]; then parmetis_inst="$(abspath "$prefix/ParMETIS_install")"; fi

echo
say "Done. Install locations:"
[[ -n "$h5_inst"     ]] && echo "  HDF5   : $h5_inst"
[[ -n "$cgns_inst"   ]] && echo "  CGNS   : $cgns_inst"
[[ -n "$petsc_inst"  ]] && echo "  PETSc  : $petsc_inst"
[[ -n "$scotch_inst" ]] && ((with_scotch )) && echo "  Scotch   : $scotch_inst"
[[ -n "$metis_inst"  ]] && ((with_metis )) && echo "  METIS   : $metis_inst"
[[ -n "$parmetis_inst" ]] && ((with_parmetis )) && echo "  ParMETIS   : $parmetis_inst"

echo
echo "Configure Loci with:"

printf '  ./configure --prefix "%s" ' "$cfg_prefix"
[[ -n "$h5_inst"    ]] && printf ' --with-hdf5  "%s" ' "$h5_inst"
[[ -n "$cgns_inst"  ]] && printf ' --with-cgns  "%s" ' "$cgns_inst"
[[ -n "$petsc_inst" ]] && printf ' --with-petsc "%s" ' "$petsc_inst"
[[ -n "$scotch_inst" ]] && ((with_scotch )) && printf ' --with-scotch "%s" ' "$scotch_inst"
[[ -n "$metis_inst"  ]] && ((with_metis )) && printf ' --with-metis "%s" ' "$metis_inst"
[[ -n "$parmetis_inst" ]] && ((with_parmetis )) && printf ' --with-parmetis "%s" ' "$parmetis_inst"
echo
echo "  make -j$(command -v nproc >/dev/null 2>&1 && nproc || echo 4)"
echo "  make install"
echo


# Show LD_LIBRARY_PATH helper for runtime linking
llp_parts=()
[[ -n "$h5_inst"    ]] && llp_parts+=("$h5_inst/lib")
[[ -n "$cgns_inst"  ]] && llp_parts+=("$cgns_inst/lib")
[[ -n "$petsc_inst" ]] && llp_parts+=("$petsc_inst/lib")
[[ -n "$scotch_inst" ]] && ((with_scotch )) && llp_parts+=("$scotch_inst/lib")
[[ -n "$metis_inst"  ]] && ((with_metis )) && llp_parts+=("$metis_inst/lib")
[[ -n "$parmetis_inst" ]] && ((with_parmetis )) && llp_parts+=("$parmetis_inst/lib")

if (( ${#llp_parts[@]} )); then
  echo "If needed, set runtime library path:"
  echo "  export LD_LIBRARY_PATH=\"$(IFS=:; echo "${llp_parts[*]}"):\$LD_LIBRARY_PATH\""
fi
