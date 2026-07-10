#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/variable.h>

/// A nested time level must follow its parent and remain navigable when another
/// time level is prepended.
TEST_CASE("variable time levels preserve their hierarchy") {
  const Loci::variable base("rho{n}");
  const Loci::variable nested("rho{n,it}");

  CHECK(base.time().before(nested.time()));
  CHECK_FALSE(nested.time().before(base.time()));
  CHECK(nested.parent() == base);

  const Loci::time_ident prepended_level("outer", Loci::time_ident());
  const Loci::variable prepended(prepended_level, nested);

  CHECK(prepended == Loci::variable("rho{outer,n,it}"));
  CHECK(prepended.parent() == Loci::variable("rho{outer,n}"));
}

/// Adding a namespace must update a parametric variable and its arguments while
/// preserving its time level.
TEST_CASE("adding a namespace updates parametric variable arguments") {
  const Loci::variable original("flux(species,temperature){n}");
  const Loci::variable namespaced = original.add_namespace("chem");

  CHECK(namespaced ==
        Loci::variable("chem@flux(chem@species,chem@temperature){n}"));
  CHECK(namespaced.time() == original.time());
}
