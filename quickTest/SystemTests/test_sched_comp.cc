#include <Loci.h>
#include <depend_graph.h>
#include "../../src/System/loci_globs.h"
#include "../../src/System/sched_tools.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <initializer_list>
#include <sstream>

using namespace Loci;

namespace Loci {
  extern bool memory_greedy_schedule;
  extern bool randomized_memory_greedy_schedule;
}

namespace {

  struct SchedulerFlagGuard {
    const bool saved_dynamic_memory;
    const bool saved_chomp;
    const bool saved_memory_greedy;
    const bool saved_randomized_memory_greedy;

    SchedulerFlagGuard()
      : saved_dynamic_memory(use_dynamic_memory),
        saved_chomp(use_chomp),
        saved_memory_greedy(memory_greedy_schedule),
        saved_randomized_memory_greedy(randomized_memory_greedy_schedule) {
      use_dynamic_memory = false;
      use_chomp = false;
      memory_greedy_schedule = false;
      randomized_memory_greedy_schedule = false;
    }

    ~SchedulerFlagGuard() {
      use_dynamic_memory = saved_dynamic_memory;
      use_chomp = saved_chomp;
      memory_greedy_schedule = saved_memory_greedy;
      randomized_memory_greedy_schedule = saved_randomized_memory_greedy;
    }
  };

  template <class TRule>
  rule make_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  variableSet vars(const char *name) {
    variableSet result;
    result += variable(name);
    return result;
  }

  variableSet vars(std::initializer_list<const char *> names) {
    variableSet result;
    for(std::initializer_list<const char *>::const_iterator i = names.begin();
        i != names.end(); ++i)
      result += variable(*i);
    return result;
  }

  entitySet span(Entity first, Entity last) {
    entitySet result;
    result += interval(first, last);
    return result;
  }

  void create_source_fact(fact_db &facts, const char *name,
                          const entitySet &domain = span(4, 6)) {
    store<int> source;
    source.allocate(domain);
    FORALL(domain, entity) {
      source[entity] = entity * 10;
    } ENDFORALL;
    facts.create_fact(name, source);
  }

  void create_constraint_fact(fact_db &facts, const char *name,
                              const entitySet &domain) {
    Loci::constraint mask;
    mask = domain;
    facts.create_fact(name, mask);
  }

  digraph graph_for_rules(std::initializer_list<rule> rules) {
    digraph graph;
    for(std::initializer_list<rule>::const_iterator i = rules.begin();
        i != rules.end(); ++i) {
      graph.add_edges(i->sources() + i->get_info().constraints(), i->ident());
      graph.add_edges(i->ident(), i->targets());
    }
    return graph;
  }

  digraph graph_for_rule(const rule &r) {
    return graph_for_rules({r});
  }

  void compile_and_analyze(graph_compiler &compiler, fact_db &facts,
                           sched_db &scheds, const variableSet &given,
                           const variableSet &target) {
    compiler.compile(facts, scheds, given, target);
    compiler.existential_analysis(facts, scheds);
    CHECK_FALSE(scheds.errors_found());
  }

  // These small pointwise rules exercise real graph_compiler/impl_compiler paths.
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

  class sched_mid_rule : public pointwise_rule {
    const_store<int> source;
    store<int> mid;
  public:
    sched_mid_rule() {
      disable_threading();
      name_store("sc_source", source);
      name_store("sc_mid", mid);
      input("sc_source");
      output("sc_mid");
    }

    void compute(const sequence &seq) {
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si)
        mid[*si] = source[*si] + 1;
    }

    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  class sched_target_from_mid_rule : public pointwise_rule {
    const_store<int> mid;
    store<int> target;
  public:
    sched_target_from_mid_rule() {
      disable_threading();
      name_store("sc_mid", mid);
      name_store("sc_target", target);
      input("sc_mid");
      output("sc_target");
    }

    void compute(const sequence &seq) {
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si)
        target[*si] = mid[*si] * 2;
    }

    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  class sched_masked_rule : public pointwise_rule {
    const_store<int> source;
    Loci::constraint mask;
    store<int> target;
  public:
    sched_masked_rule() {
      disable_threading();
      name_store("sc_source", source);
      name_store("sc_mask", mask);
      name_store("sc_masked_target", target);
      input("sc_source");
      constraint("sc_mask");
      output("sc_masked_target");
    }

    void compute(const sequence &seq) {
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si)
        target[*si] = source[*si] + 5;
    }

    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

} // namespace

TEST_CASE("graph_compiler registers concrete rules and top-level supernodes") {
  SchedulerFlagGuard flags;
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
  SchedulerFlagGuard flags;
  const rule r = make_rule<sched_comp_rule>();
  const digraph graph = graph_for_rule(r);

  fact_db facts;
  create_source_fact(facts, "sc_source");

  sched_db scheds;
  scheds.init(facts);
  set_var_types(facts, graph, scheds);

  decomposed_graph decomp(graph, vars("sc_source"), vars("sc_target"));
  graph_compiler compiler(decomp, variableSet());

  compile_and_analyze(compiler, facts, scheds, vars("sc_source"), vars("sc_target"));

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

  schedule->execute(facts, scheds);
  store<int> target(facts.get_variable("sc_target"));
  CHECK(target.domain() == span(4, 6));
  CHECK(target[4] == 45);
  CHECK(target[5] == 55);
  CHECK(target[6] == 65);
}

TEST_CASE("graph_compiler schedules dependent concrete rules in order") {
  SchedulerFlagGuard flags;
  const rule mid_rule = make_rule<sched_mid_rule>();
  const rule target_rule = make_rule<sched_target_from_mid_rule>();
  const digraph graph = graph_for_rules({mid_rule, target_rule});

  fact_db facts;
  create_source_fact(facts, "sc_source");

  sched_db scheds;
  scheds.init(facts);
  set_var_types(facts, graph, scheds);

  decomposed_graph decomp(graph, vars("sc_source"), vars("sc_target"));
  graph_compiler compiler(decomp, variableSet());

  compile_and_analyze(compiler, facts, scheds, vars("sc_source"), vars("sc_target"));

  CHECK(scheds.get_existential_info(variable("sc_mid"), mid_rule) == span(4, 6));
  CHECK(scheds.get_existential_info(variable("sc_target"), target_rule) == span(4, 6));
  CHECK(scheds.get_variable_requests(variable("sc_source")) == span(4, 6));
  CHECK(scheds.get_variable_requests(variable("sc_mid")) == span(4, 6));
  CHECK(scheds.get_variable_requests(variable("sc_target")) == span(4, 6));
  CHECK(scheds.get_exec_seq(mid_rule) == span(4, 6));
  CHECK(scheds.get_exec_seq(target_rule) == span(4, 6));

  executeP schedule = compiler.execution_schedule(facts, scheds, variableSet());
  executeP empty_schedule;
  REQUIRE(schedule != empty_schedule);

  schedule->execute(facts, scheds);
  store<int> mid(facts.get_variable("sc_mid"));
  store<int> target(facts.get_variable("sc_target"));

  CHECK(mid.domain() == span(4, 6));
  CHECK(target.domain() == span(4, 6));
  CHECK(mid[4] == 41);
  CHECK(mid[5] == 51);
  CHECK(mid[6] == 61);
  CHECK(target[4] == 82);
  CHECK(target[5] == 102);
  CHECK(target[6] == 122);
}

TEST_CASE("graph_compiler respects constraint facts when choosing execution context") {
  SchedulerFlagGuard flags;
  const rule r = make_rule<sched_masked_rule>();
  const digraph graph = graph_for_rule(r);

  fact_db facts;
  create_source_fact(facts, "sc_source");
  create_constraint_fact(facts, "sc_mask", span(5, 6));

  sched_db scheds;
  scheds.init(facts);
  set_var_types(facts, graph, scheds);

  decomposed_graph decomp(graph, vars({"sc_source", "sc_mask"}),
                          vars("sc_masked_target"));
  graph_compiler compiler(decomp, variableSet());

  compile_and_analyze(compiler, facts, scheds, vars({"sc_source", "sc_mask"}),
                      vars("sc_masked_target"));

  CHECK(scheds.get_existential_info(variable("sc_masked_target"), r) == span(5, 6));
  CHECK(scheds.variable_existence(variable("sc_mask")) == span(5, 6));
  CHECK(scheds.get_variable_requests(variable("sc_source")) == span(5, 6));
  CHECK(scheds.get_variable_requests(variable("sc_masked_target")) == span(5, 6));
  CHECK(scheds.get_exec_seq(r) == span(5, 6));

  executeP schedule = compiler.execution_schedule(facts, scheds, variableSet());
  executeP empty_schedule;
  REQUIRE(schedule != empty_schedule);

  schedule->execute(facts, scheds);
  store<int> target(facts.get_variable("sc_masked_target"));
  CHECK(target.domain() == span(5, 6));
  CHECK(target[5] == 55);
  CHECK(target[6] == 65);
}

TEST_CASE("graph_compiler handles empty concrete rule domains") {
  SchedulerFlagGuard flags;
  const rule r = make_rule<sched_comp_rule>();
  const digraph graph = graph_for_rule(r);

  fact_db facts;
  create_source_fact(facts, "sc_source", EMPTY);

  sched_db scheds;
  scheds.init(facts);
  set_var_types(facts, graph, scheds);

  decomposed_graph decomp(graph, vars("sc_source"), vars("sc_target"));
  graph_compiler compiler(decomp, variableSet());

  compile_and_analyze(compiler, facts, scheds, vars("sc_source"), vars("sc_target"));

  CHECK(scheds.get_existential_info(variable("sc_target"), r) == EMPTY);
  CHECK(scheds.get_variable_requests(variable("sc_target")) == EMPTY);
  CHECK(scheds.get_variable_requests(variable("sc_source")) == EMPTY);
  CHECK(scheds.get_exec_seq(r) == EMPTY);

  executeP schedule = compiler.execution_schedule(facts, scheds, variableSet());
  executeP empty_schedule;
  REQUIRE(schedule != empty_schedule);

  schedule->execute(facts, scheds);
  store<int> target(facts.get_variable("sc_target"));
  CHECK(target.domain() == EMPTY);
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
