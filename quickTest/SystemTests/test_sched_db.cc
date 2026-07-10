#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

/// Initializing sched_db from fact_db must make each fact domain available as
/// the corresponding variable's existence.
TEST_CASE("sched_db initializes variable existence from fact domains") {
  fact_db facts;

  const entitySet cell_domain = interval(4, 6);
  store<int> cells;
  cells.allocate(cell_domain);
  facts.create_fact("cells", cells);

  sched_db scheds(facts);

  CHECK(scheds.get_typed_variables().inSet(variable("cells")));
  CHECK(scheds.variable_existence(variable("cells")) == cell_domain);
}

/// A rule must be scheduled only where its target variable is both requested
/// and produced by that rule.
TEST_CASE("rule requests are limited to entities produced by the rule") {
  fact_db facts;

  const entitySet produced_cells = interval(10, 12);
  store<int> cells;
  cells.allocate(produced_cells);
  facts.create_fact("cells", cells);

  sched_db scheds(facts);
  const rule producer(
    "source(upstream),target(cells),qualifier(test_producer)");
  scheds.set_existential_info(variable("cells"), producer, produced_cells);
  scheds.variable_request(variable("cells"), interval(11, 13));

  const entitySet expected_rule_request = interval(11, 12);
  CHECK(scheds.get_variable_request(producer, variable("cells")) ==
        expected_rule_request);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
