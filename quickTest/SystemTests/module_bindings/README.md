# Module input binding proof

One driver loads three small shared libraries through Loci's application loader.
`consumer.loci` contains direct pressure use and a local `Copy(X)` family;
`provider.loci` computes pressure from a seed. Two consumer instances must read
their selected providers while retaining private normals and family names.
`writer.loci` supplies forbidden writes for rejection tests.

From a configured Loci source tree:

```sh
make -C OBJ -j18 library lpp SPRNG
make -C OBJ/quickTest/SystemTests/module_bindings \
  LOCI_BASE="$PWD/OBJ" check
make -C OBJ/quickTest/SystemTests/module_bindings \
  LOCI_BASE="$PWD/OBJ" check-mpi
```

The serial check is also part of `quickTest/SystemTests`. Logs and native
schedules remain in this test's build directory under `work/serial` and
`work/mpi2`. An assertion failure or nonzero process exit fails the target.
The normal Loci memory-management and chomping settings remain enabled.

The five test cases check:

1. Concrete identity resolution, preserved annotations, untouched argument
   labels and underscore-derived names, and absolute root destinations.
2. Rejected invalid tables/chains/rebinding, and table retention in copies.
3. Actual parametric expansion and query results for distinct/shared providers,
   reversed loading, changed seeds and copied databases. A no-binding control
   uses ordinary private inputs. Queries of imported spellings publish only the
   canonical provider fact. Original expanded rules retain their local inputs.
4. A parametric writer that is legal as a template but targets an imported
   input after substitution. It must throw the binding error without executing.
5. Rejected defaults, unit/apply contributions and both in-place write names.

For seeds 2 and 5 and private normals 3 and 7, the distinct-provider case gives
tractions 6 and 35, and copies 2 and 5. Changing the first seed to 4 produces
tractions 12 and 35. Its final native schedule must retain these relationships:

```text
left@traction<-a@pressure,left@normal
left@Copy(left@pressure)<-a@pressure
right@Copy(right@pressure)<-b@pressure
```

The MPI test repeats the param contract on two ranks; it is not a distributed
mesh-transfer test. Deck/setup bindings and the new module declaration API are
future work; see the [implementation note][design].

[design]: ../../../docs/developer/module-bindings.md
