#include <Loci.h>
#include <depend_graph.h>
#include "../../src/System/sched_tools.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <sstream>

using namespace Loci;

namespace {

  template <class TRule>
  rule make_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  variableSet vars(const char *name) {
    variableSet result;
    result += variable(name);
    return result;
  }

  entitySet span(Entity first, Entity last) {
    entitySet result;
    result += interval(first, last);
    return result;
  }

  void create_source_fact(fact_db &facts, const char *name) {
    store<int> source;
    source.allocate(interval(4, 6));
    source[4] = 40;
    source[5] = 50;
    source[6] = 60;
    facts.create_fact(name, source);
  }

  digraph graph_for_rule(const rule &r) {
    digraph graph;
    graph.add_edges(r.sources(), r.ident());
    graph.add_edges(r.ident(), r.targets());
    return graph;
  }

  class sched_comp_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    sched_comp_rule() {
      disable_threading();
      name_store("sc_source", source);
      name_store("sc_target", target);
      input("sc_source");
      output("sc_target");
    }

    void compute(const sequence &seq) {
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si)
        target[*si] = source[*si] + 5;
    }

    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

} // namespace

TEST_CASE("graph_compiler registers concrete rules and top-level supernodes") {
  const rule r = make_rule<sched_comp_rule>();
  const digraph graph = graph_for_rule(r);

  decomposed_graph decomp(graph, vars("sc_source"), vars("sc_target"));
  graph_compiler compiler(decomp, variableSet());

  CHECK(compiler.rule_process.find(r) != compiler.rule_process.end());
  CHECK(compiler.rule_process.find(compiler.baserule) != compiler.rule_process.end());
  CHECK_FALSE(compiler.super_rules.empty());
  CHECK(compiler.baserule.targets().inSet(variable("sc_target")));
}

TEST_CASE("graph_compiler drives existence request and execution schedule for concrete rules") {
  const rule r = make_rule<sched_comp_rule>();
  const digraph graph = graph_for_rule(r);

  fact_db facts;
  create_source_fact(facts, "sc_source");

  sched_db scheds;
  scheds.init(facts);
  set_var_types(facts, graph, scheds);

  decomposed_graph decomp(graph, vars("sc_source"), vars("sc_target"));
  graph_compiler compiler(decomp, variableSet());

  compiler.compile(facts, scheds, vars("sc_source"), vars("sc_target"));
  compiler.existential_analysis(facts, scheds);

  CHECK(scheds.get_existential_info(variable("sc_target"), r) == span(4, 6));
  CHECK(scheds.get_variable_requests(variable("sc_target")) == span(4, 6));
  CHECK(scheds.get_variable_requests(variable("sc_source")) == span(4, 6));
  CHECK(scheds.get_exec_seq(r) == span(4, 6));

  executeP schedule = compiler.execution_schedule(facts, scheds, variableSet());
  executeP empty_schedule;
  REQUIRE(schedule != empty_schedule);

  std::ostringstream printed;
  schedule->Print(printed);
  CHECK(printed.str().find("sc_target") != std::string::npos);
  CHECK(printed.str().find("over sequence") != std::string::npos);
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
