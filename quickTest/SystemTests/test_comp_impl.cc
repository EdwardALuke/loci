#include <Loci.h>

#include "../../src/System/comp_tools.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <string>
#include <vector>

using namespace Loci;

namespace Loci {
  extern int current_rule_id;
}

namespace {

  /*
    FVM regressions already cover makeQuery and full rule-graph scheduling with
    real .loci modules. These fixtures stay below that path and pin only
    comp_impl contracts that are hard to observe directly from those tests.
  */

  template <class TRule>
  rule make_loci_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  entitySet entities(Entity first, Entity last) {
    entitySet result;
    result += interval(first, last);
    return result;
  }

  void fill_store(store<int> &values, int scale) {
    entitySet dom = values.domain();
    for(entitySet::const_iterator i = dom.begin(); i != dom.end(); ++i)
      values[*i] = scale * *i;
  }

  void create_int_fact_on_2_to_4(fact_db &facts, const char *name,
                                  int scale = 10) {
    store<int> source;
    source.allocate(interval(2, 4));
    fill_store(source, scale);
    facts.create_fact(name, source);
  }

  class lifecycle_execute_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    static std::vector<std::string> events;
    static std::vector<int> sequence_seen_in_compute;
    static int rule_id_seen_in_prelude;
    static int rule_id_seen_in_compute;
    static int rule_id_seen_in_postlude;

    lifecycle_execute_rule() {
      disable_threading();
      name_store("ci_exec_src", source);
      name_store("ci_exec_tgt", target);
      input("ci_exec_src");
      output("ci_exec_tgt");
    }

    void prelude(const sequence &seq) {
      rule_id_seen_in_prelude = current_rule_id;
      events.push_back(std::string("prelude:") + std::to_string(seq.size()));
    }

    void compute(const sequence &seq) {
      rule_id_seen_in_compute = current_rule_id;
      events.push_back(std::string("compute:") + std::to_string(seq.size()));
      sequence_seen_in_compute.clear();
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
        sequence_seen_in_compute.push_back(*si);
        target[*si] = source[*si] + 5;
      }
    }

    void postlude(const sequence &seq) {
      rule_id_seen_in_postlude = current_rule_id;
      events.push_back(std::string("postlude:") + std::to_string(seq.size()));
    }
  };

  std::vector<std::string> lifecycle_execute_rule::events;
  std::vector<int> lifecycle_execute_rule::sequence_seen_in_compute;
  int lifecycle_execute_rule::rule_id_seen_in_prelude = 0;
  int lifecycle_execute_rule::rule_id_seen_in_compute = 0;
  int lifecycle_execute_rule::rule_id_seen_in_postlude = 0;

  class request_tracking_blackbox_rule : public blackbox_rule {
    const_store<int> source;
    store<int> target;
  public:
    request_tracking_blackbox_rule() {
      name_store("ci_blackbox_src", source);
      name_store("ci_blackbox_tgt", target);
      input("ci_blackbox_src");
      output("ci_blackbox_tgt");
    }

    void compute(const sequence &) {}
  };

  class observing_super_rule : public super_rule {
    store<int> target;
  public:
    static int existential_calls;
    static int request_calls;
    static int last_existential_rule_id;
    static int last_request_rule_id;

    observing_super_rule() {
      name_store("ci_super_tgt", target);
      output("ci_super_tgt");
    }

    void process_existential(rule r, fact_db &, sched_db &scheds) {
      ++existential_calls;
      last_existential_rule_id = r.ident();
      scheds.set_existential_info(variable("ci_super_tgt"), r, entities(3, 5));
    }

    void process_requests(rule r, fact_db &, sched_db &scheds) {
      ++request_calls;
      last_request_rule_id = r.ident();
      scheds.variable_request(variable("ci_super_tgt"), entities(4, 4));
    }

    void compute(const sequence &) {}
  };

  int observing_super_rule::existential_calls = 0;
  int observing_super_rule::request_calls = 0;
  int observing_super_rule::last_existential_rule_id = 0;
  int observing_super_rule::last_request_rule_id = 0;

} // namespace

TEST_CASE("execute_rule preserves lifecycle order and rule context") {
  const rule r = make_loci_rule<lifecycle_execute_rule>();

  fact_db facts;
  create_int_fact_on_2_to_4(facts, "ci_exec_src", 10);

  store<int> target_fact;
  target_fact.allocate(interval(2, 4));
  for(Entity i = 2; i <= 4; ++i)
    target_fact[i] = -1;
  facts.create_fact("ci_exec_tgt", target_fact);

  sequence reverse_exec;
  reverse_exec += interval(4, 2);

  sched_db scheds;
  lifecycle_execute_rule::events.clear();
  lifecycle_execute_rule::sequence_seen_in_compute.clear();
  lifecycle_execute_rule::rule_id_seen_in_prelude = 0;
  lifecycle_execute_rule::rule_id_seen_in_compute = 0;
  lifecycle_execute_rule::rule_id_seen_in_postlude = 0;
  current_rule_id = -1;

  execute_rule execute(r, reverse_exec, facts, scheds);
  execute.execute(facts, scheds);

  REQUIRE(lifecycle_execute_rule::events.size() == 3);
  CHECK(lifecycle_execute_rule::events[0] == "prelude:3");
  CHECK(lifecycle_execute_rule::events[1] == "compute:3");
  CHECK(lifecycle_execute_rule::events[2] == "postlude:3");

  CHECK(lifecycle_execute_rule::rule_id_seen_in_prelude == r.ident());
  CHECK(lifecycle_execute_rule::rule_id_seen_in_compute == r.ident());
  CHECK(lifecycle_execute_rule::rule_id_seen_in_postlude == r.ident());
  CHECK(current_rule_id == 0);

  REQUIRE(lifecycle_execute_rule::sequence_seen_in_compute.size() == 3);
  CHECK(lifecycle_execute_rule::sequence_seen_in_compute[0] == 4);
  CHECK(lifecycle_execute_rule::sequence_seen_in_compute[1] == 3);
  CHECK(lifecycle_execute_rule::sequence_seen_in_compute[2] == 2);

  store<int> target(facts.get_variable("ci_exec_tgt"));
  CHECK(target[2] == 25);
  CHECK(target[3] == 35);
  CHECK(target[4] == 45);
}

TEST_CASE("blackbox_compiler requests source and target extents") {
  const rule r = make_loci_rule<request_tracking_blackbox_rule>();

  fact_db facts;
  create_int_fact_on_2_to_4(facts, "ci_blackbox_src");

  store<int> target_fact;
  target_fact.allocate(interval(7, 8));
  facts.create_fact("ci_blackbox_tgt", target_fact);

  sched_db scheds;
  scheds.init(facts);

  blackbox_compiler compiler(r);
  compiler.set_var_existence(facts, scheds);
  compiler.process_var_requests(facts, scheds);

  CHECK(scheds.get_existential_info(variable("ci_blackbox_tgt"), r) == ~EMPTY);
  CHECK(scheds.get_variable_requests(variable("ci_blackbox_src")) ==
        entities(2, 4));
  CHECK(scheds.get_variable_requests(variable("ci_blackbox_tgt")) == ~EMPTY);
  CHECK(scheds.get_exec_seq(r) == ~EMPTY);

  executeP schedule = compiler.create_execution_schedule(facts, scheds);
  executeP empty_schedule;
  REQUIRE(schedule != empty_schedule);
}

TEST_CASE("superRule_compiler delegates existence and request phases") {
  const rule r = make_loci_rule<observing_super_rule>();

  fact_db facts;
  store<int> target_fact;
  target_fact.allocate(interval(3, 5));
  facts.create_fact("ci_super_tgt", target_fact);

  sched_db scheds;
  scheds.init(facts);

  observing_super_rule::existential_calls = 0;
  observing_super_rule::request_calls = 0;
  observing_super_rule::last_existential_rule_id = 0;
  observing_super_rule::last_request_rule_id = 0;

  superRule_compiler compiler(r);
  compiler.set_var_existence(facts, scheds);
  compiler.process_var_requests(facts, scheds);

  CHECK(observing_super_rule::existential_calls == 1);
  CHECK(observing_super_rule::request_calls == 1);
  CHECK(observing_super_rule::last_existential_rule_id == r.ident());
  CHECK(observing_super_rule::last_request_rule_id == r.ident());
  CHECK(scheds.get_existential_info(variable("ci_super_tgt"), r) ==
        entities(3, 5));
  CHECK(scheds.get_variable_requests(variable("ci_super_tgt")) ==
        entities(4, 4));

  executeP schedule = compiler.create_execution_schedule(facts, scheds);
  executeP empty_schedule;
  REQUIRE(schedule != empty_schedule);
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
