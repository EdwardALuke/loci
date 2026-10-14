# Module interface examples

These are dummy applications loaded into one Loci rule/fact database. Their
small formulas let the tests check interface behavior without a fluid or solid
solver. The host supplies exchange rules; each module supplies its own response.

From a configured source tree:

```sh
make -C OBJ -j18 library lpp SPRNG
make -C OBJ/quickTest/SystemTests/module_bindings \
  LOCI_BASE="$PWD/OBJ" check
make -C OBJ/quickTest/SystemTests/module_bindings \
  LOCI_BASE="$PWD/OBJ" check-mpi
```

The serial check is part of `quickTest/SystemTests`. Failures return nonzero;
crashes never satisfy negative tests. Default memory management and chomping
remain enabled. Each library and driver must use the matching rebuilt SDK.

## Two programs

`test_bindings.cc` retains the five original low-level binding checks. Its
`consumer`, `provider` and `writer` libraries exercise direct inputs, local
`Copy(X)` expansion, changing seeds, shared/distinct partners, database copies,
and rejected table/default/reduction/in-place writes.

`test_composition.loci` uses the new declaration/setup API. Its nine test cases
cover independent pairs, nested instance names, reversed declarations, repeated
initialization, private decks/defaults, shared imports, qualified references,
callback discovery, invalid declarations, forbidden writes, failed setup,
known type mismatches, canonical writes exposed by substitution, and mapped
fields on different ranks. Failure injection lives in the dummy solid only.

The application files are deliberately small:

| File | Responsibility |
| --- | --- |
| `dummy_fluid.loci` | Local `Copy(X)` family and heat-flux response. |
| `dummy_solid.loci` | Temperature response and private compliance setup. |
| `dummy_field.loci` | Gather through a bound map and constraint. |
| `canonical_writer.loci` | Forbidden provider write after substitution. |

## Exchange example

The decks select fluid temperature 400, conductance 2 and solid capacity 10.
Both modules import ambient temperature 300. Fluid setup reads that import;
solid setup reads its private deck control. Both run exactly once per instance.

For host wall temperature `Tw`, the fluid computes `q = 2*(400-Tw)` and the
solid returns `Ts = 300 + q/10`. The host's `{k}` loop feeds `Ts` into the next
`Tw`. Starting at 300 gives `300 -> 320 -> 316 -> 316.8`. Repeating the query
with initial value 310 must give 316.72. This is a fixed-point dataflow example,
not a transient heat-conduction discretization.

Native schedules under `work/loop-serial/debug` and `work/loop-mpi2/debug` show:

```text
fluid@Copy(fluid@wallTemperature){k} <- exchange@wallTemperature{k}
fluid@heatFlux{k} <- fluid@Copy(fluid@wallTemperature){k}, ...
solid@temperature{k} <- fluid@heatFlux{k}, ...
exchange@wallTemperature{k+1} <- solid@temperature{k}
```

The serial target checks the loop and bound-family schedule entries. The
fresh-process `work/load-rule-only` and `work/load-full` runs check that loading
through an older overload first does not hide the versioned setup callback.

## Mapped field example

Rank zero owns samples `[2,5]`; the last rank owns two destination entities.
The supplied map reverses their order, so the result must be `[5,2]`. The field,
map and active destination set are three distinct imported facts. At two ranks,
each rank has an empty local domain on one side of the transfer.

The check repeats with static and dynamic maps. C++ exchanges only fixture
entity identities; the Loci rule moves numerical data. `work/map-mpi2/debug`
retains both schedules. The target requires a send from rank zero and a receive
on rank one for `host@sample`. This does not establish physical conservation or
nonmatching-mesh interpolation.

Full-suite artifacts are under `work/serial` and `work/mpi2` in the build tree.
See the [interface contract][design] for initialization order, error handling,
current limits and the next work to keep within dummy applications.

[design]: ../../../docs/developer/module-bindings.md
