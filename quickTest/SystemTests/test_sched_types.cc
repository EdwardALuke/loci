#include <Loci.h>
#include <depend_graph.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

namespace Loci {
  void set_var_types(fact_db &facts, const digraph &graph, sched_db &scheds);
}

namespace {

  class store_target_rule : public pointwise_rule {
    const_store<int> source_cells;
    store<int> target_cells;

  public:
    store_target_rule() {
      name_store("source_cells", source_cells);
      name_store("target_cells", target_cells);
      input("source_cells");
      output("target_cells");
    }

    void compute(const sequence &) {}
  };

} // namespace

/// Variable type setup must add an untyped rule target to fact_db as an
/// intensional fact using the container type declared by the rule.
TEST_CASE("set_var_types installs an untyped target as an intensional fact") {
  fact_db facts;
  store<int> source_cells;
  source_cells.allocate(interval(1, 3));
  facts.create_fact("source_cells", source_cells);

  sched_db scheds(facts);
  const rule typing_rule(
    rule_implP(new copy_rule_impl<store_target_rule>));
  digraph graph;
  graph.add_edges(typing_rule.sources(), typing_rule.ident());
  graph.add_edges(typing_rule.ident(), typing_rule.targets());

  set_var_types(facts, graph, scheds);

  const storeRepP target = facts.get_variable("target_cells");
  REQUIRE(target != static_cast<storeRep *>(0));
  CHECK(isSTORE(target));
  CHECK(facts.get_intensional_facts().inSet(variable("target_cells")));
  CHECK(scheds.get_typed_variables().inSet(variable("target_cells")));
}

/// Variable type setup must preserve an existing extensional target fact,
/// including its domain and values.
TEST_CASE("set_var_types preserves an extensional target fact") {
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
  const rule typing_rule(
    rule_implP(new copy_rule_impl<store_target_rule>));
  digraph graph;
  graph.add_edges(typing_rule.sources(), typing_rule.ident());
  graph.add_edges(typing_rule.ident(), typing_rule.targets());

  set_var_types(facts, graph, scheds);

  const storeRepP target = facts.get_variable("target_cells");
  REQUIRE(target != static_cast<storeRep *>(0));
  store<int> target_after_typing(target);
  CHECK(target_after_typing.domain() == target_domain);
  CHECK(target_after_typing[7] == 41);
  CHECK(target_after_typing[8] == 42);
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
