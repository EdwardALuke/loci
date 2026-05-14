#include <Loci.h>
#include <scheduler.h>

#include "../../src/System/comp_tools.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <sstream>
#include <string>
#include <vector>

using namespace Loci;

namespace Loci {
  extern int current_rule_id;
}

namespace {

  /*
    These tests exercise comp_impl.cc with deliberately small rules and fact
    databases. Each fixture uses unique variable names so a test can ask the
    scheduler for one target without depending on any registered application
    rules.
  */

  template <class TRule>
  rule make_loci_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  template <class TRule>
  rule add_loci_rule(rule_db &rdb) {
    const rule r = make_loci_rule<TRule>();
    rdb.add_rule(r);
    return r;
  }

  variableSet single_variable(const char *name) {
    variableSet result;
    result += variable(name);
    return result;
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

  void allocate_store_on_2_to_4(store<int> &values, int value) {
    values.allocate(interval(2, 4));
    for(Entity i = 2; i <= 4; ++i)
      values[i] = value;
  }

  /*
    The concrete pointwise rules below model the part of a user rule that the
    compiler needs: named input/output stores and a deterministic compute body.
    Threading is disabled so ordering-sensitive assertions stay about
    comp_impl behavior rather than worker scheduling.
  */

  class schedule_pointwise_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    schedule_pointwise_rule() {
      disable_threading();
      name_store("ci_sched_src", source);
      name_store("ci_sched_tgt", target);
      input("ci_sched_src");
      output("ci_sched_tgt");
    }

    void compute(const sequence &seq) {
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si)
        target[*si] = source[*si] + 1;
    }
  };

  class query_pointwise_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    query_pointwise_rule() {
      disable_threading();
      name_store("ci_query_src", source);
      name_store("ci_query_tgt", target);
      input("ci_query_src");
      output("ci_query_tgt");
    }

    void compute(const sequence &seq) {
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si)
        target[*si] = source[*si] + 10;
    }
  };

  class lifecycle_pointwise_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    static std::vector<std::string> events;

    lifecycle_pointwise_rule() {
      disable_threading();
      name_store("ci_life_src", source);
      name_store("ci_life_tgt", target);
      input("ci_life_src");
      output("ci_life_tgt");
    }

    void prelude(const sequence &seq) {
      events.push_back(std::string("prelude:") + std::to_string(seq.size()));
    }

    void compute(const sequence &seq) {
      events.push_back(std::string("compute:") + std::to_string(seq.size()));
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si)
        target[*si] = source[*si] * 2;
    }

    void postlude(const sequence &seq) {
      events.push_back(std::string("postlude:") + std::to_string(seq.size()));
    }
  };

  std::vector<std::string> lifecycle_pointwise_rule::events;

  class direct_execute_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    static int rule_id_seen_in_compute;
    static std::vector<int> sequence_seen_in_compute;

    direct_execute_rule() {
      disable_threading();
      name_store("ci_direct_src", source);
      name_store("ci_direct_tgt", target);
      input("ci_direct_src");
      output("ci_direct_tgt");
    }

    void compute(const sequence &seq) {
      rule_id_seen_in_compute = current_rule_id;
      sequence_seen_in_compute.clear();
      for(sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
        sequence_seen_in_compute.push_back(*si);
        target[*si] = source[*si] + 100;
      }
    }
  };

  int direct_execute_rule::rule_id_seen_in_compute = 0;
  std::vector<int> direct_execute_rule::sequence_seen_in_compute;

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

TEST_CASE("create_execution_schedule records the pointwise execution domain") {
  rule_db rdb;
  const rule r = add_loci_rule<schedule_pointwise_rule>(rdb);

  fact_db facts;
  create_int_fact_on_2_to_4(facts, "ci_sched_src");
  facts.setupDefaults(rdb);

  sched_db scheds;
  executeP schedule =
    create_execution_schedule(rdb, facts, scheds, single_variable("ci_sched_tgt"));

  executeP empty_schedule;
  REQUIRE(schedule != empty_schedule);
  CHECK(scheds.get_existential_info(variable("ci_sched_tgt"), r) ==
        entities(2, 4));
  CHECK(scheds.get_variable_requests(variable("ci_sched_src")) ==
        entities(2, 4));
  CHECK(scheds.get_variable_requests(variable("ci_sched_tgt")) ==
        entities(2, 4));
  CHECK(scheds.get_exec_seq(r) == entities(2, 4));

  std::ostringstream printed;
  schedule->Print(printed);
  CHECK(printed.str().find("ci_sched_tgt") != std::string::npos);
  CHECK(printed.str().find("over sequence") != std::string::npos);
}

TEST_CASE("makeQuery executes a pointwise rule and installs the computed fact") {
  rule_db rdb;
  add_loci_rule<query_pointwise_rule>(rdb);

  fact_db facts;
  create_int_fact_on_2_to_4(facts, "ci_query_src");

  REQUIRE(makeQuery(rdb, facts, "ci_query_tgt"));

  store<int> target(facts.get_variable("ci_query_tgt"));
  CHECK(target.domain() == entities(2, 4));
  CHECK(target[2] == 30);
  CHECK(target[3] == 40);
  CHECK(target[4] == 50);
}

TEST_CASE("execute_rule calls prelude compute and postlude in order") {
  rule_db rdb;
  add_loci_rule<lifecycle_pointwise_rule>(rdb);

  fact_db facts;
  create_int_fact_on_2_to_4(facts, "ci_life_src", 1);

  lifecycle_pointwise_rule::events.clear();
  REQUIRE(makeQuery(rdb, facts, "ci_life_tgt"));

  REQUIRE(lifecycle_pointwise_rule::events.size() == 3);
  CHECK(lifecycle_pointwise_rule::events[0] == "prelude:3");
  CHECK(lifecycle_pointwise_rule::events[1] == "compute:3");
  CHECK(lifecycle_pointwise_rule::events[2] == "postlude:3");

  store<int> target(facts.get_variable("ci_life_tgt"));
  CHECK(target[2] == 4);
  CHECK(target[3] == 6);
  CHECK(target[4] == 8);
}

TEST_CASE("execute_rule keeps current_rule_id live and can override a target store") {
  const rule r = make_loci_rule<direct_execute_rule>();

  fact_db facts;
  create_int_fact_on_2_to_4(facts, "ci_direct_src", 3);

  store<int> original_target;
  allocate_store_on_2_to_4(original_target, -1);
  facts.create_fact("ci_direct_tgt", original_target);

  store<int> replacement_target;
  allocate_store_on_2_to_4(replacement_target, -9);
  storeRepP replacement_rep = replacement_target.Rep();

  sequence reverse_exec;
  reverse_exec += interval(4, 2);

  sched_db scheds;
  direct_execute_rule::rule_id_seen_in_compute = 0;
  direct_execute_rule::sequence_seen_in_compute.clear();
  current_rule_id = -1;

  execute_rule execute(r, reverse_exec, facts, variable("ci_direct_tgt"),
                       replacement_rep, scheds);
  execute.execute(facts, scheds);

  CHECK(direct_execute_rule::rule_id_seen_in_compute == r.ident());
  CHECK(current_rule_id == 0);
  REQUIRE(direct_execute_rule::sequence_seen_in_compute.size() == 3);
  CHECK(direct_execute_rule::sequence_seen_in_compute[0] == 4);
  CHECK(direct_execute_rule::sequence_seen_in_compute[1] == 3);
  CHECK(direct_execute_rule::sequence_seen_in_compute[2] == 2);

  CHECK(replacement_target[2] == 106);
  CHECK(replacement_target[3] == 109);
  CHECK(replacement_target[4] == 112);
  CHECK(original_target[2] == -1);
  CHECK(original_target[3] == -1);
  CHECK(original_target[4] == -1);
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
  CHECK(schedule->getName() == "execute_rule");
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
  CHECK(schedule->getName() == "execute_rule");
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
