# Explicit module interfaces and bindings

This branch provides an opt-in composition API, tested with dummy applications.
It combines rules, decks and setup in one named interface. Stream and FEM are
not migrated. The existing loader and its namespace conventions remain intact.

## A module and its instances

A library exports a versioned description of its local input and output names:

```cpp
extern "C" const Loci::module_interface_v1 *loci_module_interface_v1() {
  static const Loci::module_interface_v1 interface = {
    {"wallTemperature", "ambient"}, {"heatFlux"}
  } ;
  return &interface ;
}
```

Its rules use those names normally. The host declares instances and bindings:

```cpp
using Loci::variable ;
Loci::module_composition applications ;
applications.add("binding_dummy_fluid", {"fluid"}, {
  {variable("wallTemperature"), variable("exchange@wallTemperature")},
  {variable("ambient"), variable("environment@temperature")}
}, "fluid.vars") ;
applications.add("binding_dummy_solid", {"solid"}, {
  {variable("heatFlux"), variable("fluid@heatFlux")},
  {variable("ambient"), variable("environment@temperature")}
}, "solid.vars") ;
applications.prepare(host_rules) ;
applications.initialize(facts) ;
Loci::makeQuery(applications.rules(), facts, "finalTemperature") ;
```

The host supplies ordinary facts and rules for its exchange loop; see the
[complete example][example]. The scope is a vector of components, so
`{"outer", "solid"}` means `outer@solid` without underscore conventions.

Every declared input requires exactly one binding. Unknown/missing inputs,
repeated interface names, input/output overlap and conflicting instance
specifications are errors. Repeating an identical instance declaration before
preparation is harmless. Output declarations document the public names; they do
not hide private facts from Loci queries or guarantee a satisfiable provider.

Bindings use concrete names without arguments, priorities or time annotations.
Destinations are absolute; an unqualified destination deliberately names a root
fact. Several imports may share one destination. Destinations cannot be other
binding keys, excluding chains and cycles. Bindings become fixed at preparation.
A provider need not have storage until the query computes it.

## Preparation and setup

Preparation loads cached library definitions, clones rules per instance,
qualifies each variable, and validates the complete binding table and rules.
One explicitly qualified source does not exempt other variables from relocation.
The new callback is discovered from the library handle regardless of which
legacy loader overload opened the library first. No legacy initializer is
invoked by the new composition API.

Initialization uses an unscoped host fact database. It performs these phases:

1. Check existing imported data and known input types; install defaults.
2. Read all instance decks using their prepared declarations.
3. Invoke each `loci_init_module_v1(module_context&)` callback in declaration
   order, then check the resulting imported data types.

All decks precede all callbacks. Reversing declarations does not change name
resolution, but arbitrary callback dependencies still need explicit ordering.
A callback cannot query a partner merely to obtain a value: computed inputs
become available through the host's later query.

`module_context` provides `resolve`, `read`, `create`, `update` and `erase`.
Relative names are local to the instance; imported reads resolve to the chosen
provider. Writes are checked before storage access and must stay within the
instance and outside its imported inputs, including explicit provider spellings.
Callbacks should use const Loci containers when reading inputs. The context is
valid during the callback and is not a sandbox against arbitrary C++ pointer
mutation or host-created shared storage.

The deck reader reuses `fact_db::read_vars` through an explicit write resolver.
It checks a name before `Input` can modify existing data. Ordinary private
controls keep default/optional behavior. Qualified names within the instance,
including nested scopes, are supported; namespace-switch blocks are rejected.
Raw legacy deck/setup calls do not acquire this policy automatically.

Repeated initialization on the same live fact database does not rerun callbacks.
Initialization failure makes that composition unusable. Earlier private writes
or callback side effects are not rolled back. Prepare/initialize failures throw
`module_binding_error` where they arise from the interface contract; arbitrary
callback exceptions propagate after marking the composition failed.

## Rule resolution and ownership

Existing parametric expansion runs first. For a local family:

```text
Authored:           Copy(X) <- X
Expanded:          fluid@Copy(fluid@wallTemperature) <- fluid@wallTemperature
Bound:             fluid@Copy(fluid@wallTemperature) <- exchange@wallTemperature
```

The argument remains a substitution label; the standalone input identifies
numerical data. Binding preserves time, offset, assignment and priority.
It does not rewrite argument labels, namespace symbols or derived names such as
`wallTemperature_ic`. Existing substitution and promotion decisions are intact.

The resolution pass uses structured variables and existing rule cloning. It
runs after `parametric_rdb` and before map processing, relation generation and
scheduling. Original templates remain reusable. Concrete imported query targets
resolve before extensional-fact checks; results use canonical provider names.
There is no duplicate numerical input store.

The rule database retains every originating instance, including when identical
rules are interned. Expansion carries this metadata to instantiated rules.
Checks before defaults and after expansion reject imported-input writes through
ordinary rules, defaults, unit/apply targets and either in-place write name.
Instance ownership also detects writes spelled as an imported destination;
legitimate provider instances retain their own writes. Query-time binding errors
propagate before schedule execution. The host must handle a failed query.

Known concrete input declarations are compared with each other and supplied
facts using empty frozen representations. This accepts the tested static/dynamic
map forms without copying numerical fields. Types revealed only by substitution
retain normal Loci checks; units and physical compatibility remain application
responsibilities.

The lookup pass adds `O(N log(B+1))` lookup work for `N` occurrences and `B`
bindings, plus existing copying/interning costs and finite origin metadata.
It emits at most one rule per input rule and no new parametric requests. This is
not a complexity claim about the existing expansion or scheduler.

## Evidence and remaining scope

The [tests][tests] include two dummy applications with decks, callbacks, a
host-owned exchange loop, changed initial data and nested instance names.
Separate mapped-field checks bind a store, map and constraint across two owning
ranks, including empty local domains and static/dynamic maps. Native schedules
retain the loop and the send/receive of the bound field.

On 2026-10-14 the core library, the original five binding cases and nine new
composition cases passed in serial and on two ranks with normal memory
management and chomping. Fresh processes verified callback discovery after each
legacy loading overload. The existing namespace catalogue retained 28 passes
and its two documented failures, `nested_scope_restoration` and
`hierarchical_deck`; the parametric-rule and time-level groups also passed.

This is a framework interface demonstration, not a physical CHT/FSI simulation.
The next useful work remains in these dummy applications: time-history and
restart contracts, interface convergence controls, representative transfers,
and setup failures that occur on only one MPI rank.
Do not migrate solvers until that contract is reviewed.

Current limits are explicit:

- Instances and dependencies are declared by the host. Recursive dependency
  loading and legacy `register_module` declarations are not supported here.
- Dynamic rule constraints are rejected in the new mode because the existing
  renaming path does not cover them. GPU execution is not qualified.
- The examples run in one Loci composition, not separate processes coordinating
  independent solver jobs. The host owns the iteration structure.
- Hosts must use the same composition on all participating ranks and coordinate
  setup failure before entering collective queries. The negative tests inject
  the same failure on every rank; asymmetric failure recovery is not qualified.
- Known native container types are checked; this is not a complete type system
  for generated families, arbitrary third-party containers or physical units.
- Header/object layouts changed. Build libraries and hosts against the same
  SDK; no binary compatibility with an older Loci installation is promised.

[example]: ../../quickTest/SystemTests/module_bindings/test_composition.loci
[tests]: ../../quickTest/SystemTests/module_bindings/README.md
