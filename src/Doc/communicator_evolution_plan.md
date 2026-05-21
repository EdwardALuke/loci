# Communicator Evolution Plan

## Motivation

The Loci framework currently hardcodes `MPI_COMM_WORLD` as the
communicator for all parallel operations.  This prevents building
parallel schedules over processor subsets --- a capability needed for
coupled multi-physics, sub-domain solvers, and ensemble runs.

The framework already accepts an `MPI_Comm` parameter in many
lower-level APIs (e.g. `redistribute`, `GLOBAL_OR`, `hdf5CreateFile`),
but those parameters default to `MPI_COMM_WORLD`.  The `fact_db`
class, which is the central data manager, has no communicator member at
all.  The global variables `MPI_processes` and `MPI_rank` are likewise
derived from `MPI_COMM_WORLD` during `Loci::Init()`.

The plan below describes a phased, backward-compatible evolution that:

1. Stores a communicator in `fact_db` and derives rank/size from it.
2. Threads the communicator through APIs that currently hardcode
   `MPI_COMM_WORLD`.
3. Introduces a compile-time feature flag (`LOCI_STRICT_COMM`) that
   removes the `MPI_COMM_WORLD` defaults, forcing solver codes to
   pass communicators explicitly.

Each phase is designed so that **existing solver codes compile and run
without modification** until a project opts in to strict mode.

---

## Phase 0 — Communicator in `fact_db`

### 0.1  Add an `MPI_Comm` member to `fact_db`

Add a private member and public accessors:

```cpp
private:
  MPI_Comm comm_ ;
public:
  MPI_Comm get_comm() const { return comm_ ; }
  void set_comm(MPI_Comm c) ;
```

The constructor initializes `comm_` to `MPI_COMM_WORLD`.
`set_comm` duplicates the communicator (via `MPI_Comm_dup`) so that
`fact_db` owns its own copy, and caches the corresponding rank and
size in two new private members (`comm_rank_`, `comm_size_`).

### 0.2  Local rank/size accessors

```cpp
int get_comm_rank() const { return comm_rank_ ; }
int get_comm_size() const { return comm_size_ ; }
```

These replace direct reads of the globals `Loci::MPI_rank` and
`Loci::MPI_processes` for any code that has access to a `fact_db`.

### 0.3  Propagate the communicator inside `fact_db`

Replace the two `MPI_Alltoall` calls in `fact_db::get_dist_alloc()`
(currently hardcoded to `MPI_COMM_WORLD`) with calls that use
`comm_`.  Similarly update `fact_db::put_distribute_info()` and any
allocation helpers that reference the globals.

### 0.4  Propagate to `distribute_info`

Add an `MPI_Comm comm` member to the `distribute_info` struct so that
distributed metadata can carry the communicator to redistribution and
I/O routines downstream.

---

## Phase 1 — Thread the communicator through internal APIs

In this phase the communicator stored in `fact_db` is passed down to
internal subsystems.  Because these APIs already accept an optional
`MPI_Comm` parameter with a default of `MPI_COMM_WORLD`, the changes
are mechanical: pass `facts.get_comm()` instead of relying on the
default.  Note to accomplish this when no fact_db is present create a
get_exec_comm() function that returns the communicator from the
exec_current_fact_db unless this pointer is null, then return
MPI_COMM_WORLD.

### 1.1  Scheduler (`scheduler.cc`)

The scheduler holds a pointer to the active `fact_db`.  All
`MPI_Allreduce` calls used for profiling and barrier-style operations
(approximately 30 call sites) should be updated to use the
communicator retrieved from the fact database.

### 1.2  Grid reader and partitioning (`FVMGridReader.cc`, `distribute.cc`)

The grid reader passes `MPI_COMM_WORLD` to collective operations,
partition routines (`ORBPartition`, `parSampleSort`,
`balanceDistribution`, ParMETIS wrappers), and I/O helpers.  At each
call site, replace the literal `MPI_COMM_WORLD` with the communicator
obtained from the `fact_db` that is already available in context.

For the partitioning utilities themselves, the signatures already
accept `MPI_Comm`; no signature changes are needed.

### 1.3  Data redistribution (`distribute_container.cc`, `Map.cc`, etc.)

Container redistribution helpers call `MPI_Alltoall` /
`MPI_Alltoallv` with `MPI_COMM_WORLD`.  These functions should accept
an `MPI_Comm` parameter (defaulting to `MPI_COMM_WORLD`) and the call
sites in the framework should forward the `fact_db` communicator.

### 1.4  Computation support (`comp_tools.cc`, `comp_reduce.cc`)

Point-to-point communication (`MPI_Irecv`, `MPI_Send`) and packing
routines use `MPI_COMM_WORLD`.  Add an `MPI_Comm` parameter to
affected routines and forward it from the execution context, which has
access to the `fact_db`.

### 1.5  I/O subsystem (`distribute_io.cc`, `multiStoreIO.cc/.h`,
`mpi_multiStoreIO.h`)

The core I/O functions (`hdf5CreateFile`, `hdf5OpenFile`,
`writeContainer`, `readContainer`, `readContainerRAW`, etc.) already
accept `MPI_Comm`.  The overloads without a communicator argument
forward to `MPI_COMM_WORLD`.  No signature changes are needed here;
the callers should simply start passing `facts.get_comm()`.

### 1.6  Dynamic execution (`dynamic_exec.cc`)

Dynamic rule execution uses `MPI_Send`/`MPI_Recv`/`MPI_Probe` with
`MPI_COMM_WORLD`.  The executing context carries a `fact_db` pointer;
extract the communicator from it.

### 1.7  Pack/Unpack in store representations

The templated `MPI_Pack`/`MPI_Unpack` calls in `gstore_impl.h`,
`gparameter_impl.h`, `gmultistore_impl.h`, `gmapvec_impl.h`, and GPU
variants use `MPI_COMM_WORLD`.  Since packing is communicator-
independent in MPI (the communicator is only used to look up
datatypes), these can remain as `MPI_COMM_WORLD` or be parameterized
for consistency.

### 1.8  PETSc integration (`petsc.loci`)

PETSc objects (`VecCreate`, `MatCreate`, `KSPCreate`,
`PetscViewerASCIIOpen`) are created with `MPI_COMM_WORLD`.  These
should use the communicator from the `fact_db` available through the
Loci rule context.  This change enables PETSc solvers to operate on
processor subsets transparently.

### 1.9  Overset and adaptation modules

`FVMOverset/interpolation.loci`, `FVMOverset/motion.loci`,
`FVMAdapt/library/write_vog.cc`, and `FVMAdapt/library/color_matrix.cc`
contain many `MPI_COMM_WORLD` references.  These modules operate on
rule-provided stores and have access to the `fact_db`.  Forward the
communicator at each call site.

### 1.10  Tool programs (`FVMtools/`)

Standalone conversion tools (`ugrid2vog.cc`, `vogtools.cc`,
`fluent2vog.cc`, `cgns2vog.cc`, etc.) use `MPI_COMM_WORLD` directly.
Because these are not library code consumed by solver developers, they
can be updated independently and do not affect the public API.

---

## Phase 2 — Update global state management

### 2.1  Context-local rank and size

Today the globals `Loci::MPI_processes` and `Loci::MPI_rank` are set
once in `Loci::Init()` from `MPI_COMM_WORLD`.  After Phase 0, code
inside the framework that has access to a `fact_db` should prefer
`facts.get_comm_rank()` and `facts.get_comm_size()`.

To avoid breaking existing solver code that reads these globals, they
will continue to exist and remain set to the `MPI_COMM_WORLD` values.
A deprecation warning should be emitted (or documented) advising
solvers to migrate to the `fact_db` accessors.

### 2.2  `Loci::SetDefaultComm()` post-initialization

Since `Loci::Init()` calls `MPI_Init()` internally, a communicator
cannot be passed to it (MPI communicators are not valid before
`MPI_Init`).  Instead, a separate function is provided:

```cpp
void SetDefaultComm(MPI_Comm comm) ;
```

This must be called **after** `Loci::Init()`.  It updates the global
`MPI_processes` and `MPI_rank` to reflect the given communicator.
Solvers that wish to run on a sub-communicator call this before
setting up the `fact_db`.

### 2.3  `exec_current_fact_db` communicator forwarding

The global `exec_current_fact_db` pointer is used by convenience
functions like `collect_entitySet`.  These should be updated to pull
the communicator from `exec_current_fact_db->get_comm()` rather than
assuming `MPI_COMM_WORLD`.

---

## Phase 3 — Feature flag for strict communicator mode

### 3.1  Define `LOCI_STRICT_COMM`

Introduce a compile-time macro following the pattern established by
`LOCI_COMPAT_MODE1`.  When `LOCI_STRICT_COMM` is defined, the
`= MPI_COMM_WORLD` default is removed from every public API function
signature.  This forces downstream solver codes to supply a
communicator explicitly, ensuring correctness when running on a
sub-communicator.

The affected public APIs include:

| Header | Functions / Methods |
|--------|---------------------|
| `store_rep.h` | `storeRep::redistribute()` (both overloads) |
| `distribute.h` | `GLOBAL_OR`, `GLOBAL_AND`, `GLOBAL_MAX`, `GLOBAL_MIN` |
| `distribute_io.h` | `hdf5CreateFile`, `hdf5OpenFile`, overloads of `writeContainer`, `readContainer`, etc. |
| `distribute_long.h` | `g_GLOBAL_MAX`, `g_GLOBAL_MIN`, `g_all_collect_entitySet`, etc. |
| `partition.h` | `transposePtn` |
| All concrete store/map headers | `redistribute`, `recompose`, `split_redistribute`, `expand`, `inplace_compose` |

Implementation strategy using the preprocessor:

```cpp
#ifdef LOCI_STRICT_COMM
#define LOCI_DEFAULT_COMM
#else
#define LOCI_DEFAULT_COMM = MPI_COMM_WORLD
#endif
```

Then each declaration becomes:

```cpp
int GLOBAL_OR(int b, MPI_Comm comm LOCI_DEFAULT_COMM) ;
```

### 3.2  Convenience overloads under compat mode

The existing convenience overloads that omit `MPI_Comm` (e.g. the
one-argument `hdf5OpenFile`, `all_collect_vectors(e)`) should be
wrapped in `#ifndef LOCI_STRICT_COMM` guards.  Under strict mode,
only the explicit-communicator versions are available.

### 3.3  Solver migration path

Solver codes that wish to support sub-communicators can adopt
`LOCI_STRICT_COMM` at their own pace:

1. Add `-DLOCI_STRICT_COMM` to their build flags.
2. Fix all resulting compilation errors by passing
   `facts.get_comm()` where a communicator is now required.
3. Verify that the solver works correctly with a sub-communicator by
   running on a subset of MPI ranks.

Solver codes that do not define `LOCI_STRICT_COMM` continue to
compile and run exactly as before.

---

## Phase 4 — Validation and testing

### 4.1  Unit tests

Write tests that create a `fact_db` on a sub-communicator
(`MPI_Comm_split`) and exercise:

- Entity allocation (`get_distributed_alloc`).
- Store redistribution.
- Collective reductions (`GLOBAL_OR`, etc.).
- HDF5 parallel I/O.
- PETSc object creation.

### 4.2  Integration test with an existing solver

Run an existing Loci-based solver (e.g. CHEM) in two modes:

- **Full communicator**: all ranks participate.  Verify bit-identical
  results with the pre-change code.
- **Sub-communicator**: split ranks into two groups, run independent
  solver instances on each.  Verify correctness.

### 4.3  Performance regression

Ensure that passing the communicator through the extra parameter does
not introduce measurable overhead.  The communicator is a lightweight
handle; no performance impact is expected.

---

## Summary of API changes

| Change | Backward compatible? | Requires solver update? |
|--------|:--------------------:|:-----------------------:|
| `fact_db::get_comm()` / `set_comm()` | Yes | No |
| `fact_db::get_comm_rank()` / `get_comm_size()` | Yes | No |
| `Loci::SetDefaultComm(comm)` | Yes (new function) | No |
| Internal call sites forwarding communicator | Yes | No |
| `LOCI_STRICT_COMM` flag removing defaults | Opt-in | Only if flag is set |
| Deprecation of `Loci::MPI_rank` / `MPI_processes` | Yes (still available) | Recommended migration |

---

## File impact estimate

| Area | Files | Approximate call sites |
|------|------:|-----------------------:|
| `fact_db` (header + impl) | 2 | ~5 |
| Scheduler | 1 | ~30 |
| Grid reader / partitioner | 3 | ~80 |
| Data redistribution | 4 | ~30 |
| Computation / reduction | 3 | ~25 |
| I/O subsystem | 5 | ~40 |
| Dynamic execution | 1 | ~50 |
| Store/map representations (headers) | ~20 | ~60 (default parameter macro) |
| PETSc integration | 1 | ~10 |
| Overset / adaptation | 6 | ~50 |
| Tool programs | 8 | ~100 |
| Initialization | 1 | ~5 |
| **Total** | **~55** | **~485** |

---

## Recommended implementation order

1. **Phase 0** — smallest change, highest value.  Enables all
   downstream work and is independently useful for framework
   developers.
2. **Phase 2** — global state management must align with Phase 0
   before internal subsystems are updated.
3. **Phase 1** — bulk of the work; can be done subsystem by subsystem
   in any order.  Each subsystem can be validated independently.
4. **Phase 3** — introduce the feature flag only after Phase 1 is
   complete and tested.
5. **Phase 4** — continuous throughout, with a final comprehensive
   pass after Phase 3.
