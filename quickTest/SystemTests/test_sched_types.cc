#include <Loci.h>
#include <depend_graph.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <sstream>
#include <string>

using namespace Loci;

namespace Loci {
  void set_var_types(fact_db &facts, const digraph &dg, sched_db &scheds);
}

namespace {

  class CErrCapture {
    std::streambuf *old_;
    std::ostringstream captured_;
  public:
    CErrCapture() : old_(std::cerr.rdbuf(captured_.rdbuf())) {}
    ~CErrCapture() { std::cerr.rdbuf(old_); }

    std::string str() const { return captured_.str(); }
  };

  template <class TRule>
  rule make_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  template <class TRule>
  rule add_rule(rule_db &rdb) {
    const rule r = make_rule<TRule>();
    rdb.add_rule(r);
    return r;
  }

  void add_public_rule_edges(const rule &r, digraph &graph) {
    graph.add_edges(r.sources(), r.ident());
    graph.add_edges(r.ident(), r.targets());
  }

  variableSet vars(const char *name) {
    variableSet result;
    result += variable(name);
    return result;
  }

  void create_source_fact(fact_db &facts, const char *name) {
    store<int> source;
    source.allocate(interval(1, 3));
    facts.create_fact(name, source);
  }

  void type_graph(fact_db &facts, const digraph &graph, sched_db &scheds) {
    scheds.init(facts);
    set_var_types(facts, graph, scheds);
  }

  class sched_source_to_target_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    sched_source_to_target_rule() {
      name_store("st_source", source);
      name_store("st_target", target);
      input("st_source");
      output("st_target");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  class sched_priority_target_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    sched_priority_target_rule() {
      name_store("st_priority_source", source);
      name_store("p::st_priority_target", target);
      input("st_priority_source");
      output("p::st_priority_target");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  class sched_rename_target_rule : public pointwise_rule {
    const_store<int> source;
    store<int> old_value;
    store<int> new_value;
  public:
    sched_rename_target_rule() {
      name_store("st_rename_source", source);
      name_store("st_old_value", old_value);
      name_store("st_new_value", new_value);
      input("st_rename_source");
      output("st_new_value=st_old_value");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

} // namespace

//----------------------------------------------------------------------------
// Regression Tests For Known Bugs
//----------------------------------------------------------------------------

// Current regression:
// If an assignment target aliases an already-typed fact that is not present in
// the dependency graph, `set_var_types()` emits an "undefined var_to_cluster"
// diagnostic before it falls back to aliasing the variable.
//
// Suggested fix in `src/System/sched_types.cc`:
// Seed `var_to_cluster`/`cluster_map` with typed variables from `fact_db`, or
// handle rename sources outside the graph without reporting a scheduler
// diagnostic when their type is already known.
TEST_CASE("rename assignment to an existing fact does not emit missing-cluster diagnostic"
          * doctest::may_fail()) {
  rule_db rdb;
  add_rule<sched_rename_target_rule>(rdb);

  fact_db facts;
  create_source_fact(facts, "st_rename_source");
  create_source_fact(facts, "st_old_value");

  const digraph graph =
    dependency_graph2(rdb, vars("st_rename_source"), vars("st_new_value")).get_graph();

  sched_db scheds;
  CErrCapture capture;
  type_graph(facts, graph, scheds);

  CHECK(capture.str().find("undefined var_to_cluster") == std::string::npos);
}

//----------------------------------------------------------------------------
// Additional behavior/invariant tests
//----------------------------------------------------------------------------

TEST_CASE("set_var_types creates an intensional store type for a graph target") {
  rule_db rdb;
  add_rule<sched_source_to_target_rule>(rdb);

  fact_db facts;
  create_source_fact(facts, "st_source");

  const digraph graph =
    dependency_graph2(rdb, vars("st_source"), vars("st_target")).get_graph();

  sched_db scheds;
  type_graph(facts, graph, scheds);

  storeRepP target = facts.get_variable("st_target");
  REQUIRE(target != static_cast<storeRep *>(0));
  CHECK(isSTORE(target));
  CHECK(target->domain() == EMPTY);
  CHECK(facts.get_typed_variables().inSet(variable("st_target")));
  CHECK(facts.get_intensional_facts().inSet(variable("st_target")));
  CHECK_FALSE(facts.get_extensional_facts().inSet(variable("st_target")));
  CHECK(scheds.try_get_aliases(variable("st_target")).inSet(variable("st_target")));
  CHECK(scheds.get_typed_variables().inSet(variable("st_target")));
}

TEST_CASE("set_var_types links priority-qualified targets as synonyms of the base target") {
  rule_db rdb;
  add_rule<sched_priority_target_rule>(rdb);

  fact_db facts;
  create_source_fact(facts, "st_priority_source");

  const digraph graph =
    dependency_graph2(rdb, vars("st_priority_source"), vars("st_priority_target")).get_graph();

  sched_db scheds;
  type_graph(facts, graph, scheds);

  CHECK(facts.get_typed_variables().inSet(variable("st_priority_target")));
  CHECK(facts.get_typed_variables().inSet(variable("p::st_priority_target")));
  CHECK(scheds.try_get_synonyms(variable("st_priority_target")).inSet(variable("p::st_priority_target")));
  CHECK(scheds.try_get_synonyms(variable("p::st_priority_target")).inSet(variable("st_priority_target")));
}

TEST_CASE("set_var_types turns target assignment expressions into scheduler aliases") {
  fact_db facts;
  create_source_fact(facts, "st_rename_source");
  create_source_fact(facts, "st_old_value");

  rule r = make_rule<sched_rename_target_rule>();
  digraph graph;
  add_public_rule_edges(r, graph);
  graph.add_edge(variable("st_old_value").ident(), variable("st_old_value").ident());

  sched_db scheds;
  type_graph(facts, graph, scheds);

  CHECK(facts.get_typed_variables().inSet(variable("st_new_value")));
  CHECK(scheds.try_get_aliases(variable("st_old_value")).inSet(variable("st_new_value")));
  CHECK(scheds.try_get_antialiases(variable("st_new_value")).inSet(variable("st_old_value")));
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
