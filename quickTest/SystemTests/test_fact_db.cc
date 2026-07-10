#include <list>

#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

/// Clearing intensional facts must preserve extensional facts and their values.
TEST_CASE("erasing intensional facts preserves extensional facts") {
  fact_db facts;

  param<int> extensional;
  *extensional = 11;
  facts.create_fact("extensional", extensional);

  param<int> intensional;
  *intensional = 22;
  facts.create_intensional_fact("intensional", intensional);

  CHECK(facts.get_extensional_facts().inSet(variable("extensional")));
  CHECK(facts.get_intensional_facts().inSet(variable("intensional")));

  facts.erase_intensional_facts();

  storeRepP extensional_rep = facts.get_variable("extensional");
  REQUIRE(extensional_rep != static_cast<storeRep *>(0));
  param<int> extensional_after_erase(extensional_rep);
  CHECK(*extensional_after_erase == 11);
  CHECK(facts.get_variable("intensional") == static_cast<storeRep *>(0));
}

/// A variable must be found by its local name in the active namespace and by
/// its namespace-qualified name after leaving that namespace.
TEST_CASE("namespaced variables support local and qualified lookup") {
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

/// At the end of an iteration, rotating variables must cycle their values
/// through the requested variable order.
TEST_CASE("rotating variables cycles their values in the requested order") {
  fact_db facts;

  param<int> first;
  *first = 1;
  facts.create_fact("first", first);

  param<int> second;
  *second = 2;
  facts.create_fact("second", second);

  param<int> third;
  *third = 3;
  facts.create_fact("third", third);

  std::list<variable> rotation_order;
  rotation_order.push_back(variable("first"));
  rotation_order.push_back(variable("second"));
  rotation_order.push_back(variable("third"));

  facts.rotate_vars(rotation_order);

  param<int> rotated_first(facts.get_variable("first"));
  param<int> rotated_second(facts.get_variable("second"));
  param<int> rotated_third(facts.get_variable("third"));

  CHECK(*rotated_first == 3);
  CHECK(*rotated_second == 1);
  CHECK(*rotated_third == 2);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
