# Explicit module input bindings: initial implementation

This branch implements the name-resolution step of a proposed module interface.
It proves direct inputs and local parametric families using real loaded modules.
It does not yet provide interface declarations, a new loader, or binding-aware
deck parsing and setup callbacks. Do not migrate solver applications to it yet.

## Contract

The host attaches a table to a `rule_db` before querying it:

```cpp
Loci::module_bindings inputs({
  {Loci::variable("plate@pressure"),
   Loci::variable("exchange@pressureOnSolid")}
}) ;
rules.set_module_bindings(inputs) ;
```

Both sides are concrete names without arguments, time levels or priorities.
Keys are qualified local names. Destinations are absolute identities; an
unqualified destination deliberately names a root fact. The destination can be
computed by rules and need not have storage when the table is constructed.

Several imports may share one destination. A destination cannot itself be a
binding key, which excludes chains and cycles. Installing the same table again
is harmless; replacing a nonempty table is an error. Copies of `rule_db` retain
the table. Hosts remain responsible for constructing it before setup or queries.

Bindings resolve whole fact occurrences. They preserve time, offset, assignment
and priority annotations. They do not rewrite arguments, substitute namespace
symbols, or change derived names such as `pressure_ic`. For example:

```text
Authored family:    Copy(X) <- X
Authored request:   copied <- Copy(pressure)
Loaded instance:   plate@copied <- plate@Copy(plate@pressure)
After expansion:   plate@Copy(plate@pressure) <- plate@pressure
After binding:     plate@Copy(plate@pressure) <- exchange@pressureOnSolid
```

The local argument remains a substitution label. The separate input occurrence
identifies the data to read. Neither namespace-valued parameters nor whole-rule
iteration promotion acquire new semantics.

## Integration

`makeQuery` checks explicit imported-input writers before `setupDefaults`.
After the existing `parametric_rdb` expansion, it checks again for writers
exposed by substitution and resolves a copy of the expanded rules. This precedes
map-constraint replacement, GPU rewriting, relation generation and scheduling.
Original rules remain reusable. The pass uses structured `variable` identities
and existing rule cloning/renaming; it does not rewrite rule strings.

A concrete imported query target is resolved before checking extensional facts.
Results are stored under canonical names. Code using raw `fact_db` reads the
provider name, or calls `inputs.resolve(...)` to obtain it; no duplicate input
store is created. Ordinary unbound databases skip rule resolution.

Explicit rule/default/unit/apply targets and either side of an in-place target
cannot name an imported local input. Optional declarations are not writers.
Invalid bindings throw `module_binding_error`, derived from
`std::invalid_argument`; `makeQuery` propagates this error to the caller before
executing a schedule. Hosts must treat a failed composition as unusable.

The lookup work is `O(N log(B+1))` for `N` visited occurrences and `B` bindings,
in addition to existing interning, cloning and database-insertion costs.
The pass emits at most one rule per input rule and no new parametric requests.
This is not a complexity claim about the existing expansion or scheduler.

## Deliberate limits

- The table is a low-level primitive, not a declaration of required inputs or
  public outputs. It cannot diagnose an omitted or misspelled interface input.
- No module-origin metadata is retained yet. Writes through the imported local
  spelling are rejected; this does not establish ownership of explicit provider
  names. Provider rules must still be allowed to write their own facts.
- Raw deck reads and setup callbacks do not consult the table. Do not create
  facts under imported local names. The next stage must reject such writes and
  resolve reads through a shared setup context, before changing stored data.
- Databases using dynamic rule constraints are rejected when bindings are
  enabled: existing descriptor renaming does not cover those constraints.
- Maps and ordinary constraints use existing descriptor renaming. Distributed
  bound fields, empty owning domains, time-loop execution and type compatibility
  still need dedicated acceptance cases. The current two-rank test uses params.
- The legacy module loader is unchanged, including its existing namespace and
  callback behavior. The tests use its full application-loading overload.
- `rule_db` has new data members. Rebuild the host and modules against the same
  SDK; this work does not preserve binary compatibility with an older library.

## Evidence and next stage

See [the focused tests](../../quickTest/SystemTests/module_bindings/README.md).
On 2026-10-14, a fresh build based on `origin/dev` at `e2da2fe4e` passed five
test cases (86 assertions) in serial and on each of two MPI ranks, with normal
dynamic memory management and chomping. The generated schedule contains:

```text
left@Copy(left@pressure)<-a@pressure
right@Copy(right@pressure)<-b@pressure
```

The existing `stream_cases` namespace catalogue was run unchanged against this
build: 28 passed; `nested_scope_restoration` and `hierarchical_deck` failed as
previously documented. Those issues are not repaired or waived by this branch.
This is a framework proof, not a Stream/FEM coupling qualification.

The next stage should introduce declared interfaces and per-instance provenance,
then a common resolver/write policy for loading, decks and setup callbacks.
Qualify descriptor coverage and distributed stores before migrating a solver
adapter. Mesh transfer and iteration coordination remain numerical operations.
