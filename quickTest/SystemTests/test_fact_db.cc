#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

/// Clearing derived (intensional) facts must preserve supplied (extensional)
/// facts and their values.
TEST_CASE("erasing intensional facts preserves extensional facts") {
  fact_db facts;

  param<int> supplied;
  *supplied = 11;
  facts.create_fact("supplied", supplied);

  param<int> derived;
  *derived = 22;
  facts.create_intensional_fact("derived", derived);

  CHECK(facts.get_extensional_facts().inSet(variable("supplied")));
  CHECK(facts.get_intensional_facts().inSet(variable("derived")));

  facts.erase_intensional_facts();

  storeRepP supplied_rep = facts.get_variable("supplied");
  REQUIRE(supplied_rep != static_cast<storeRep *>(0));
  param<int> supplied_after_erase(supplied_rep);
  CHECK(*supplied_after_erase == 11);
  CHECK(facts.get_variable("derived") == static_cast<storeRep *>(0));
}

/// A namespace must support local lookup while active and qualified lookup
/// after the namespace is removed.
TEST_CASE("fact namespaces support local and qualified lookup") {
  fact_db facts;
  facts.set_namespace("fluid");

  param<int> density;
  *density = 17;
  facts.create_fact("rho", density);

  storeRepP local_rep = facts.get_variable("rho");
  REQUIRE(local_rep != static_cast<storeRep *>(0));
  param<int> local_density(local_rep);
  CHECK(*local_density == 17);

  facts.remove_namespace();
  CHECK(facts.get_variable("rho") == static_cast<storeRep *>(0));

  storeRepP qualified_rep = facts.get_variable("fluid@rho");
  REQUIRE(qualified_rep != static_cast<storeRep *>(0));
  param<int> qualified_density(qualified_rep);
  CHECK(*qualified_density == 17);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
