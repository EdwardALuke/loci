#include <Loci.h>
#include <depend_graph.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

namespace Loci {
  void set_var_types(fact_db &facts, const digraph &graph, sched_db &scheds);
}

namespace {

  class target_type_rule : public pointwise_rule {
    const_store<int> source_cells;
    store<int> target_cells;

  public:
    target_type_rule() {
      name_store("source_cells", source_cells);
      name_store("target_cells", target_cells);
      input("source_cells");
      output("target_cells");
    }

    void compute(const sequence &) {}
  };

} // namespace

/// Variable type setup must create an intensional fact for an untyped target
/// variable using the variable type provided by a rule that generates it.
TEST_CASE("set_var_types creates an intensional fact for an untyped target "
          "variable") {
  fact_db facts;
  store<int> source_cells;
  source_cells.allocate(interval(1, 3));
  facts.create_fact("source_cells", source_cells);

  sched_db scheds(facts);
  const rule generating_rule(
    rule_implP(new copy_rule_impl<target_type_rule>));
  digraph graph;
  graph.add_edges(generating_rule.sources(), generating_rule.ident());
  graph.add_edges(generating_rule.ident(), generating_rule.targets());

  set_var_types(facts, graph, scheds);

  const storeRepP target_rep = facts.get_variable("target_cells");
  REQUIRE(target_rep != static_cast<storeRep *>(0));
  CHECK(isSTORE(target_rep));
  CHECK(facts.get_intensional_facts().inSet(variable("target_cells")));
  CHECK(scheds.get_typed_variables().inSet(variable("target_cells")));
}

/// Variable type setup must preserve an existing extensional fact for a target
/// variable, including its domain and values.
TEST_CASE("set_var_types preserves an extensional fact for a target variable") {
  fact_db facts;
  store<int> source_cells;
  source_cells.allocate(interval(1, 3));
  facts.create_fact("source_cells", source_cells);

  const entitySet target_domain = interval(7, 8);
  store<int> target_cells;
  target_cells.allocate(target_domain);
  target_cells[7] = 41;
  target_cells[8] = 42;
  facts.create_fact("target_cells", target_cells);

  sched_db scheds(facts);
  const rule generating_rule(
    rule_implP(new copy_rule_impl<target_type_rule>));
  digraph graph;
  graph.add_edges(generating_rule.sources(), generating_rule.ident());
  graph.add_edges(generating_rule.ident(), generating_rule.targets());

  set_var_types(facts, graph, scheds);

  const storeRepP target_rep = facts.get_variable("target_cells");
  REQUIRE(target_rep != static_cast<storeRep *>(0));
  store<int> preserved_target(target_rep);
  CHECK(preserved_target.domain() == target_domain);
  CHECK(preserved_target[7] == 41);
  CHECK(preserved_target[8] == 42);
  CHECK(facts.get_extensional_facts().inSet(variable("target_cells")));
  CHECK_FALSE(facts.get_intensional_facts().inSet(variable("target_cells")));
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
