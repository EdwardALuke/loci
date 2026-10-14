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

## Three programs

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

`test_semantics.loci` checks against separately authored explicit rules:
equivalent dependencies/results under renaming and disjoint extension, and
equivalent responses to a changing input inside a host iteration.

The application files are deliberately small:

| File | Responsibility |
| --- | --- |
| `dummy_fluid.loci` | Local `Copy(X)` family and heat-flux response. |
| `dummy_solid.loci` | Temperature response and private compliance setup. |
| `dummy_field.loci` | Gather through a bound map and constraint. |
| `canonical_writer.loci` | Forbidden provider write after substitution. |
| `response.loci` | Imported input plus a bias chosen by a local label. |

## Binding semantics example

The response module computes `Response(X) = X + bias_X`, with private
`bias_input=3`. Its `input` is bound to `source@value`, twice a seed.
The explicit reference writes `source@value` directly and uses its own bias.
The comparison checks both expanded dependency signatures and numerical values.

Seed 2 must give 7 and seed 5 must give 13. Root facts named `input` and
`bias_input`, and a partner fact `source@bias_input`, contain different values
to expose accidental capture. The same checks run after renaming/nesting the
instance and after adding an unconnected instance with the same local names.
This static test loads only stationary host rules. The next test supplies the
initial state and termination control required by the host loop.

The iteration check advances the seed from the module's preceding response.
Starting at 2 gives responses `7,17,37,77`; restarting the query from 5 gives
`13,29,61,125`. The explicit reference must agree at the final iteration. Native
schedules in `work/semantics-serial` and `work/semantics-mpi2` show the promoted
source and both response rules. The targets check the loop, changing input and
inherited bias in the native schedule. This is ordinary host iteration, not
independent solver-loop composition. The [draft specification][spec] states
the assumptions and separates these checks from open language-design decisions.

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
[spec]: ../../../docs/developer/module-binding-semantics.md
