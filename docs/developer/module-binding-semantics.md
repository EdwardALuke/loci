# Draft specification: explicit module bindings

## Goal and scope

Loci modules should compose with little adapter code and clear dependencies.
A module names its local inputs and supplies ordinary rules. The host chooses
the input providers. Loci constructs the computation from those relationships.
The common case needs no copy rule, readiness flag, or module-specific loop.

This is a proposal for review on `namespace_feature`, motivated by
[issue 468](https://github.com/EdwardALuke/loci/issues/468). It specifies the
limited binding feature demonstrated here, not a replacement namespace language
or an approved upstream contract. The [implementation guide][guide] describes
the C++ API, setup phases, diagnostics and unsupported constructs.

| Concept | Responsibility |
| --- | --- |
| Namespace | Identify facts and separate instance-local names. |
| Interface | Declare local inputs and advertised outputs. |
| Binding | Select the concrete fact supplying one local input. |
| Ordinary rules | Compute, transfer, reduce and express iteration structure. |

## 1. A binding selects identity

For an instance named `plate`, the host may declare:

```text
plate@input -> source@value
```

The arrow is specification notation, not new Loci rule syntax. Every declared
input has exactly one destination. Several inputs may share a destination.
Binding endpoints are concrete names without arguments, priorities or time
annotations. An unqualified destination names a root fact.

Selection is fixed before scheduling. Availability of another fact or provider
cannot redirect the binding. The selected fact may be supplied as data or
computed by rules; its producers retain normal Loci selection semantics.
An unavailable selected fact leaves an unsatisfied dependency, without fallback
to a similarly named root, local or enclosing-scope fact.

A binding adds no numerical store and no numerical operation. The host must
connect compatible entity domains, representations, units and reference frames.
Native container-type checks establish only part of that compatibility.
Different meshes require explicit maps/transfer rules whose result is bound.

## 2. Substitution labels remain local

Consider this module family and consumer:

```text
Response(X) <- X, bias_X
result      <- Response(input)
```

Substitution first selects the local label `input`. Binding then redirects its
standalone data reference to `source@value`. The private `bias_input` and the
local family result keep their identities. In the current branch's spelling:

```text
plate@Response(plate@input) <- source@value, plate@bias_input
plate@result               <- plate@Response(plate@input)
```

Binding `input` does not import a partner's `bias_value`, gradient, initial
condition, or family of rules. Those are separate facts requiring explicit
relationships or bindings. Arguments are substitution labels, not a mechanism
for asking the host to search a partner's namespace.

The loader currently qualifies argument representations along with local fact
names. The examples retain that spelling. Whether local labels should instead
be represented without qualification remains a review decision; success of
these examples does not settle general qualified-argument semantics.

## 3. Time and rule meaning are preserved

Let `B` be the binding table, `x` a concrete input name, and `a` its time,
offset, assignment and priority annotations. Resolution is:

```text
resolve(x with a) = B(x) with a       when x is an input
resolve(x with a) = x with a          otherwise
```

Argument-bearing expressions and their labels are not rewritten by this pass.
Maps, constraints and conditionals that reference a bound concrete fact must
refer to the same destination as ordinary rule inputs. No loop index is renamed.

The order is instance qualification, existing parametric expansion, concrete
binding resolution, then existing relation generation/promotion and scheduling.
Imported query names also resolve before checking for supplied data. These
operations need not commute; their order is part of this proposal.

Binding `plate@input{n,it}` selects `source@value{n,it}`. It does not choose a
previous iteration, create a snapshot, synchronize independent loops or infer
physical time. The host expresses those relationships through ordinary Loci
rules. A stationary relationship can then be promoted using existing semantics.

## 4. Imports are read-only through the interface

An importing instance must not produce or modify its imported fact, through
either its local spelling or the destination spelling. Check writes before
defaults/setup and again after substitution reveals concrete rule targets.
This covers ordinary producers, defaults, unit/apply and both in-place names.
The provider's own writers remain eligible under normal Loci rules.

Shared reduction contributions are a separate design question. For this
profile, modules expose their contributions as separate outputs and host rules
combine them. Invalid bindings/writes are errors, not requests to select a
different provider or to cancel another module's contribution.

Rules, deck writes and setup operations must use the same identity table.
Decks cannot assign imports. Setup may read available imported data through
const containers; computations requiring query results belong in rules.
Callbacks remain imperative and ordered. Reordering callbacks is safe only
when their effects and data dependencies are independent.

This is an interface policy, not C++ memory isolation. Output declarations do
not hide other facts or prove that a provider exists. Explicit foreign rule
references are still possible, so declarations do not yet enumerate every
dependency. Do not infer complete encapsulation from this API.

## 5. Bounded resolution and its limits

For this profile, a destination cannot itself be a binding key. Thus chains and
cycles are rejected, resolution is idempotent, and one lookup is sufficient per
concrete occurrence. Bindings are frozen for the composition's lifetime.

For `N` occurrences in already expanded rules and `B` bindings, table lookup
cost is `O(N log(B+1))`, plus rule copying/interning and origin metadata costs.
The pass emits at most one rule per input rule and requests no new parametric
instantiations. This bound assumes expansion has completed; it is not a bound
on Loci's expansion, promotion or scheduling as a whole.

Sharing destinations can introduce ordinary dependency cycles or incompatible
uses. Fast name resolution does not prove that a composition is schedulable.
Future re-exports could use finite alias-graph traversal with cycle detection;
the current chain ban is a scope restriction, not a claim that nesting is
fundamentally unsafe.

## 6. Review criteria and executable evidence

For deterministic rules with declared dependencies and no hidden cross-instance
effects, the intended criteria are:

| Criterion | Expected observation |
| --- | --- |
| Explicit-wiring equivalence | Same concrete dependencies and values. |
| Consistent instance renaming | Same behavior after renaming all references. |
| Disjoint extension | An unconnected instance does not change the query. |
| Temporal fidelity | Inputs vary or are inherited as in explicit rules. |
| Ownership | Substitution cannot conceal a forbidden imported write. |

[test_semantics.loci][semantics] checks the first four with one small response
module and a separately authored explicit reference. It compares concrete
dependency signatures as well as results, uses misleading root/partner names,
renames and nests an instance, adds a disjoint instance, and repeats a query
inside a changing host iteration. These examples support the criteria; they
are not a proof for every Loci construct or a physical coupling qualification.
In particular, disjoint extension is tested for a stationary module. Unused
applications containing independently configured loops require separate study.

The existing [binding tests][tests] cover preserved annotations, table errors,
ownership after substitution, decks/setup and mapped MPI data. Missing-binding
declarations are checked before setup. General unsatisfied-query diagnostics
remain Loci's existing behavior.

Before extending the feature, review three decisions: representation of local
argument labels, whether declarations should enforce dependency/privacy limits,
and how reusable compositions should re-export inputs. Defaults, optional
imports and shared-contribution interfaces must not introduce silent fallback.
Real solver migration remains outside this stage.

[guide]: module-bindings.md
[semantics]: ../../quickTest/SystemTests/module_bindings/test_semantics.loci
[tests]: ../../quickTest/SystemTests/module_bindings/README.md
