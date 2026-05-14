#include <Loci.h>
#include "../../src/System/comp_tools.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <initializer_list>
#include <list>
#include <utility>
#include <vector>

using namespace Loci;

namespace {

  /*
    comp_tools.cc contains the small pieces that convert rule descriptors and
    scheduler state into executable entity sets.  These tests keep the fact
    database serial and tiny so failures point at vmap/request logic instead of
    distributed communication setup.
  */

  template <class TRule>
  rule make_loci_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  entitySet set_of(std::initializer_list<Entity> values) {
    entitySet result;
    for(std::initializer_list<Entity>::const_iterator i = values.begin();
        i != values.end(); ++i)
      result += *i;
    return result;
  }

  variableSet single_variable(const char *name) {
    variableSet result;
    result += variable(name);
    return result;
  }

  /*
    ct_src intentionally includes both rule-context entities and mapped-image
    entities.  That lets one small fixture cover source maps and target maps.
  */
  entitySet mapped_source_domain() {
    return set_of({1, 2, 3, 10, 20});
  }

  void add_store(fact_db &facts, const char *name, const entitySet &domain) {
    store<int> values;
    values.allocate(domain);
    facts.create_fact(name, values);
  }

  void add_map(fact_db &facts, const char *name) {
    // 2 and 3 both map to 20 so preimage paths exercise fan-in.
    Map map;
    map.allocate(interval(1, 3));
    map[1] = 10;
    map[2] = 20;
    map[3] = 20;
    facts.create_fact(name, map);
  }

  /*
    These rules are metadata fixtures: comp_tools only reads their variable,
    mapping, conditional, constraint, and output declarations.
  */

  // ct_tgt[i] depends on ct_src[ct_m[i]] and the local condition ct_cond[i].
  class mapped_source_rule : public pointwise_rule {
    const_Map map;
    const_store<int> source;
    const_store<int> condition;
    store<int> target;
  public:
    mapped_source_rule() {
      name_store("ct_m", map);
      name_store("ct_src", source);
      name_store("ct_cond", condition);
      name_store("ct_tgt", target);
      input("ct_m->ct_src");
      conditional("ct_cond");
      output("ct_tgt");
    }

    void compute(const sequence &) {}
  };

  // ct_tgt[ct_m[i]] depends on ct_src[i].
  class mapped_target_rule : public pointwise_rule {
    const_Map map;
    const_store<int> source;
    store<int> target;
  public:
    mapped_target_rule() {
      name_store("ct_m", map);
      name_store("ct_src", source);
      name_store("ct_tgt", target);
      input("ct_src");
      output("ct_m->ct_tgt");
    }

    void compute(const sequence &) {}
  };

  // ct_constraint_tgt[i] is only valid when ct_constraint[ct_m[i]] exists.
  class mapped_constraint_rule : public pointwise_rule {
    const_Map map;
    const_store<int> source;
    const_store<int> constraint_store;
    store<int> target;
  public:
    mapped_constraint_rule() {
      name_store("ct_m", map);
      name_store("ct_src", source);
      name_store("ct_constraint", constraint_store);
      name_store("ct_constraint_tgt", target);
      input("ct_src");
      constraint("ct_m->ct_constraint");
      output("ct_constraint_tgt");
    }

    void compute(const sequence &) {}
  };

  // Blackbox rules request their whole target extent before computing context.
  class mapped_blackbox_rule : public blackbox_rule {
    const_Map map;
    const_store<int> source;
    const_store<int> condition;
    store<int> target;
  public:
    mapped_blackbox_rule() {
      name_store("ct_m", map);
      name_store("ct_src", source);
      name_store("ct_cond", condition);
      name_store("ct_tgt", target);
      input("ct_m->ct_src");
      conditional("ct_cond");
      output("ct_tgt");
    }

    void compute(const sequence &) {}
  };

  struct CompToolsFixture {
    fact_db facts;
    sched_db scheds;

    CompToolsFixture() {
      add_map(facts, "ct_m");
      add_store(facts, "ct_src", mapped_source_domain());
      add_store(facts, "ct_cond", set_of({1, 2, 3}));
      add_store(facts, "ct_constraint", set_of({20}));
      add_store(facts, "ct_tgt", EMPTY);
      add_store(facts, "ct_constraint_tgt", EMPTY);
      add_store(facts, "ct_other", EMPTY);
      scheds.init(facts);
    }
  };

} // namespace

TEST_CASE("schedule_dag emits dependency waves after available vertices") {
  digraph graph;
  graph.add_edge(1, 3);
  graph.add_edge(2, 3);
  graph.add_edge(3, 4);
  graph.add_edge(5, 4);

  std::vector<digraph::vertexSet> schedule =
    schedule_dag(graph, EMPTY, interval(UNIVERSE_MIN, UNIVERSE_MAX));

  REQUIRE(schedule.size() == 3);
  CHECK(schedule[0] == set_of({1, 2, 5}));
  CHECK(schedule[1] == set_of({3}));
  CHECK(schedule[2] == set_of({4}));

  std::vector<digraph::vertexSet> with_available_roots =
    schedule_dag(graph, set_of({1, 5}), interval(UNIVERSE_MIN, UNIVERSE_MAX));

  REQUIRE(with_available_roots.size() == 3);
  CHECK(with_available_roots[0] == set_of({2}));
  CHECK(with_available_roots[1] == set_of({3}));
  CHECK(with_available_roots[2] == set_of({4}));
}

TEST_CASE("existential_rule_analysis and process_rule_requests follow source mappings") {
  CompToolsFixture fixture;
  rule rule_under_test = make_loci_rule<mapped_source_rule>();

  existential_rule_analysis(rule_under_test, fixture.facts, fixture.scheds);
  CHECK(fixture.scheds.get_existential_info(variable("ct_tgt"),
                                            rule_under_test) ==
        set_of({1, 2, 3}));

  fixture.scheds.variable_request(variable("ct_tgt"), set_of({2, 3}));
  CHECK(process_rule_requests(rule_under_test, fixture.facts, fixture.scheds) ==
        set_of({2, 3}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_src")) ==
        set_of({20}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_cond")) ==
        set_of({2, 3}));
}

TEST_CASE("constraint mappings restrict existential context") {
  CompToolsFixture fixture;
  rule rule_under_test = make_loci_rule<mapped_constraint_rule>();

  CHECK(input_variables(rule_under_test).inSet(variable("ct_src")));
  CHECK(input_variables(rule_under_test).inSet(variable("ct_constraint")));
  CHECK(input_variables_with_mapping(rule_under_test).inSet(
          variable("ct_constraint")));
  CHECK_FALSE(input_variables_with_mapping(rule_under_test).inSet(
                variable("ct_src")));

  existential_rule_analysis(rule_under_test, fixture.facts, fixture.scheds);
  CHECK(fixture.scheds.get_existential_info(variable("ct_constraint_tgt"),
                                            rule_under_test) ==
        set_of({2, 3}));
}

TEST_CASE("vmap existence helpers map source and target domains") {
  CompToolsFixture fixture;

  vmap_info mapped_source_vmap;
  mapped_source_vmap.var += variable("ct_src");
  mapped_source_vmap.mapping.push_back(single_variable("ct_m"));

  CHECK(vmap_source_exist(mapped_source_vmap, fixture.facts, fixture.scheds) ==
        set_of({1, 2, 3}));

  vmap_info mapped_target_vmap;
  mapped_target_vmap.var += variable("ct_tgt");
  mapped_target_vmap.mapping.push_back(single_variable("ct_m"));

  CHECK(vmap_target_exist(mapped_target_vmap, fixture.facts, set_of({1, 3}),
                          fixture.scheds) ==
        set_of({10, 20}));
}

TEST_CASE("output mappings expose metadata and pull requests back to context") {
  CompToolsFixture fixture;
  rule rule_under_test = make_loci_rule<mapped_target_rule>();

  CHECK(rule_has_mapping_in_output(rule_under_test));
  CHECK(input_variables(rule_under_test).inSet(variable("ct_src")));
  CHECK_FALSE(input_variables_with_mapping(rule_under_test).inSet(
                variable("ct_src")));

  existential_rule_analysis(rule_under_test, fixture.facts, fixture.scheds);
  CHECK(fixture.scheds.get_existential_info(variable("ct_tgt"),
                                            rule_under_test) ==
        set_of({10, 20}));

  fixture.scheds.variable_request(variable("ct_tgt"), set_of({20}));
  CHECK(process_rule_requests(rule_under_test, fixture.facts, fixture.scheds) ==
        set_of({2, 3}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_src")) ==
        set_of({2, 3}));
}

TEST_CASE("intensive output map detection compares target mapping structure") {
  CompToolsFixture fixture;
  rule mapped_target = make_loci_rule<mapped_target_rule>();
  rule mapped_source = make_loci_rule<mapped_source_rule>();

  CHECK_FALSE(is_intensive_rule_output_mapping(mapped_target, fixture.facts));

  std::vector<variableSet> intensive_mapping;
  intensive_mapping.push_back(single_variable("ct_m"));
  fixture.facts.intensive_output_maps.insert(intensive_mapping);

  CHECK(is_intensive_rule_output_mapping(mapped_target, fixture.facts));
  CHECK_FALSE(is_intensive_rule_output_mapping(mapped_source, fixture.facts));
}

TEST_CASE("vmap request helpers map target requests backward and source requests forward") {
  CompToolsFixture fixture;

  vmap_info target_vmap;
  target_vmap.var += variable("ct_tgt");
  target_vmap.var += variable("ct_other");

  vdefmap unmapped_requests;
  unmapped_requests[variable("ct_tgt")] = set_of({1});
  unmapped_requests[variable("ct_other")] = set_of({3});

  CHECK(vmap_target_requests(target_vmap, unmapped_requests, fixture.facts,
                             fixture.scheds, false) ==
        set_of({1, 3}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_tgt")) == EMPTY);
  CHECK(fixture.scheds.get_variable_requests(variable("ct_other")) == EMPTY);

  CHECK(vmap_target_requests(target_vmap, unmapped_requests, fixture.facts,
                             fixture.scheds) ==
        set_of({1, 3}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_tgt")) ==
        set_of({1, 3}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_other")) ==
        set_of({1, 3}));

  vmap_info mapped_target_vmap;
  mapped_target_vmap.var += variable("ct_tgt");
  mapped_target_vmap.mapping.push_back(single_variable("ct_m"));

  vdefmap mapped_requests;
  mapped_requests[variable("ct_tgt")] = set_of({10});

  CHECK(vmap_target_requests(mapped_target_vmap, mapped_requests,
                             fixture.facts, fixture.scheds, false) ==
        set_of({1}));

  mapped_requests[variable("ct_tgt")] = set_of({20});

  CHECK(vmap_target_requests(mapped_target_vmap, mapped_requests,
                             fixture.facts, fixture.scheds, false) ==
        set_of({2, 3}));

  mapped_requests[variable("ct_tgt")] = set_of({10, 20});

  CHECK(vmap_target_requests(mapped_target_vmap, mapped_requests,
                             fixture.facts, fixture.scheds, false) ==
        set_of({1, 2, 3}));

  vmap_info mapped_source_vmap;
  mapped_source_vmap.var += variable("ct_src");
  mapped_source_vmap.mapping.push_back(single_variable("ct_m"));

  CHECK(vmap_source_requests(mapped_source_vmap, fixture.facts, set_of({1, 2}),
                             fixture.scheds) ==
        set_of({10, 20}));
}

TEST_CASE("blackbox request analysis uses full target extent and mapped sources") {
  CompToolsFixture fixture;
  rule rule_under_test = make_loci_rule<mapped_blackbox_rule>();

  existential_blackboxrule_analysis(rule_under_test, fixture.facts,
                                    fixture.scheds);
  CHECK(fixture.scheds.get_existential_info(variable("ct_tgt"),
                                            rule_under_test) ==
        set_of({1, 2, 3}));

  CHECK(process_blackboxrule_requests(rule_under_test, fixture.facts,
                                      fixture.scheds) ==
        set_of({1, 2, 3}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_tgt")) ==
        set_of({1, 2, 3}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_src")) ==
        set_of({10, 20}));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_cond")) ==
        set_of({1, 2, 3}));
}

TEST_CASE("communication helpers are no-ops for non-distributed facts") {
  CompToolsFixture fixture;
  std::list<comm_info> clist;

  CHECK(send_requests(set_of({1, 2}), variable("ct_src"), fixture.facts,
                      clist) == EMPTY);
  CHECK(clist.empty());

  std::vector<std::pair<variable, entitySet> > sends;
  sends.push_back(std::make_pair(variable("ct_src"), set_of({1, 2})));
  CHECK(put_precomm_info(sends, fixture.facts).empty());

  variableSet request_vars = single_variable("ct_src");
  fixture.scheds.variable_request(variable("ct_src"), set_of({1}));
  CHECK(barrier_process_rule_requests(request_vars, fixture.facts,
                                      fixture.scheds).empty());
  CHECK(fixture.scheds.get_variable_requests(variable("ct_src")) ==
        set_of({1}));
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
