#include <Loci.h>
#include "../../src/System/comp_tools.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <list>
#include <utility>
#include <vector>

using namespace Loci;

namespace {

  template <class TRule>
  rule make_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  entitySet ids(Entity a) {
    entitySet result;
    result += a;
    return result;
  }

  entitySet ids(Entity a, Entity b) {
    entitySet result;
    result += a;
    result += b;
    return result;
  }

  entitySet ids(Entity a, Entity b, Entity c) {
    entitySet result;
    result += a;
    result += b;
    result += c;
    return result;
  }

  entitySet source_domain() {
    entitySet result = ids(1, 2, 3);
    result += ids(10, 20);
    return result;
  }

  variableSet vars(const char *name) {
    variableSet result;
    result += variable(name);
    return result;
  }

  void add_store(fact_db &facts, const char *name, const entitySet &domain) {
    store<int> s;
    s.allocate(domain);
    facts.create_fact(name, s);
  }

  void add_map(fact_db &facts, const char *name) {
    Map m;
    m.allocate(interval(1, 3));
    m[1] = 10;
    m[2] = 20;
    m[3] = 20;
    facts.create_fact(name, m);
  }

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
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

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
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  struct CompToolsFixture {
    fact_db facts;
    sched_db scheds;

    CompToolsFixture() {
      add_map(facts, "ct_m");
      add_store(facts, "ct_src", source_domain());
      add_store(facts, "ct_cond", ids(1, 2, 3));
      add_store(facts, "ct_tgt", EMPTY);
      add_store(facts, "ct_other", EMPTY);
      scheds.init(facts);
    }
  };

} // namespace

TEST_CASE("existential_rule_analysis and process_rule_requests follow source mappings") {
  CompToolsFixture fixture;
  rule r = make_rule<mapped_source_rule>();

  existential_rule_analysis(r, fixture.facts, fixture.scheds);
  CHECK(fixture.scheds.get_existential_info(variable("ct_tgt"), r) == ids(1, 2, 3));

  fixture.scheds.variable_request(variable("ct_tgt"), ids(2, 3));
  CHECK(process_rule_requests(r, fixture.facts, fixture.scheds) == ids(2, 3));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_src")) == ids(20));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_cond")) == ids(2, 3));
}

TEST_CASE("output mappings expose metadata and pull requests back to context") {
  CompToolsFixture fixture;
  rule r = make_rule<mapped_target_rule>();

  CHECK(rule_has_mapping_in_output(r));
  CHECK(input_variables(r).inSet(variable("ct_src")));
  CHECK_FALSE(input_variables_with_mapping(r).inSet(variable("ct_src")));

  existential_rule_analysis(r, fixture.facts, fixture.scheds);
  CHECK(fixture.scheds.get_existential_info(variable("ct_tgt"), r) == ids(10, 20));

  fixture.scheds.variable_request(variable("ct_tgt"), ids(20));
  CHECK(process_rule_requests(r, fixture.facts, fixture.scheds) == ids(2, 3));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_src")) == ids(2, 3));
}

TEST_CASE("vmap_target_requests can compute context without mutating requests") {
  CompToolsFixture fixture;
  vmap_info vmi;
  vmi.var += variable("ct_tgt");
  vmi.var += variable("ct_other");

  vdefmap requests;
  requests[variable("ct_tgt")] = ids(1);
  requests[variable("ct_other")] = ids(3);

  CHECK(vmap_target_requests(vmi, requests, fixture.facts, fixture.scheds, false) == ids(1, 3));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_tgt")) == EMPTY);
  CHECK(fixture.scheds.get_variable_requests(variable("ct_other")) == EMPTY);

  CHECK(vmap_target_requests(vmi, requests, fixture.facts, fixture.scheds) == ids(1, 3));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_tgt")) == ids(1, 3));
  CHECK(fixture.scheds.get_variable_requests(variable("ct_other")) == ids(1, 3));
}

TEST_CASE("communication helpers are no-ops for non-distributed facts") {
  CompToolsFixture fixture;
  std::list<comm_info> clist;

  CHECK(send_requests(ids(1, 2), variable("ct_src"), fixture.facts, clist) == EMPTY);
  CHECK(clist.empty());

  std::vector<std::pair<variable, entitySet> > sends;
  sends.push_back(std::make_pair(variable("ct_src"), ids(1, 2)));
  CHECK(put_precomm_info(sends, fixture.facts).empty());

  variableSet request_vars = vars("ct_src");
  fixture.scheds.variable_request(variable("ct_src"), ids(1));
  CHECK(barrier_process_rule_requests(request_vars, fixture.facts, fixture.scheds).empty());
  CHECK(fixture.scheds.get_variable_requests(variable("ct_src")) == ids(1));
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
