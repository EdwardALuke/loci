#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

/// Initializing sched_db from a fact_db records each fact's domain as the
/// corresponding variable's existence.
TEST_CASE("sched_db records fact domains as variable existence") {
  fact_db facts;

  const entitySet cell_domain = interval(4, 6);
  store<int> cells;
  cells.allocate(cell_domain);
  facts.create_fact("cells", cells);

  sched_db scheds(facts);

  CHECK(scheds.get_typed_variables().inSet(variable("cells")));
  CHECK(scheds.variable_existence(variable("cells")) == cell_domain);
}

/// A rule's request for a target variable is the overlap between the requested
/// entities and the entities where that rule can compute the target.
TEST_CASE("rule request contains requested entities where the rule can "
          "compute its target") {
  fact_db facts;

  const entitySet target_cells = interval(10, 12);
  const entitySet requested_cells = interval(11, 13);
  store<int> cells;
  cells.allocate(target_cells);
  facts.create_fact("cells", cells);

  sched_db scheds(facts);
  const rule cell_rule(
    "source(upstream),target(cells),qualifier(cell_rule)");
  scheds.set_existential_info(variable("cells"), cell_rule, target_cells);
  scheds.variable_request(variable("cells"), requested_cells);

  const entitySet expected_rule_request = interval(11, 12);
  CHECK(scheds.get_variable_request(cell_rule, variable("cells")) ==
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
