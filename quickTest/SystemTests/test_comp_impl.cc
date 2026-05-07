#include <Loci.h>
#include <scheduler.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <sstream>
#include <string>
#include <vector>

using namespace Loci;

namespace {

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

  void create_input_fact(fact_db &facts, const char *name, int scale = 10) {
    store<int> source;
    source.allocate(interval(2, 4));
    source[2] = scale * 2;
    source[3] = scale * 3;
    source[4] = scale * 4;
    facts.create_fact(name, source);
  }

  class concrete_schedule_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    concrete_schedule_rule() {
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

    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  class concrete_query_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    concrete_query_rule() {
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

    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  class concrete_lifecycle_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    static std::vector<std::string> events;

    concrete_lifecycle_rule() {
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

    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  std::vector<std::string> concrete_lifecycle_rule::events;

} // namespace

TEST_CASE("create_execution_schedule records execution sequence for a concrete rule") {
  rule_db rdb;
  const rule r = add_rule<concrete_schedule_rule>(rdb);

  fact_db facts;
  create_input_fact(facts, "ci_sched_src");
  facts.setupDefaults(rdb);

  sched_db scheds;
  executeP schedule = create_execution_schedule(rdb, facts, scheds, vars("ci_sched_tgt"));

  executeP empty_schedule;
  REQUIRE(schedule != empty_schedule);
  CHECK(scheds.get_existential_info(variable("ci_sched_tgt"), r) == span(2, 4));
  CHECK(scheds.get_variable_requests(variable("ci_sched_src")) == span(2, 4));
  CHECK(scheds.get_variable_requests(variable("ci_sched_tgt")) == span(2, 4));
  CHECK(scheds.get_exec_seq(r) == span(2, 4));

  std::ostringstream printed;
  schedule->Print(printed);
  CHECK(printed.str().find("ci_sched_tgt") != std::string::npos);
  CHECK(printed.str().find("over sequence") != std::string::npos);
}

TEST_CASE("makeQuery executes a concrete rule and installs the computed fact") {
  rule_db rdb;
  add_rule<concrete_query_rule>(rdb);

  fact_db facts;
  create_input_fact(facts, "ci_query_src");

  REQUIRE(makeQuery(rdb, facts, "ci_query_tgt"));

  store<int> target(facts.get_variable("ci_query_tgt"));
  CHECK(target.domain() == span(2, 4));
  CHECK(target[2] == 30);
  CHECK(target[3] == 40);
  CHECK(target[4] == 50);
}

TEST_CASE("execute_rule calls prelude compute and postlude in order") {
  rule_db rdb;
  add_rule<concrete_lifecycle_rule>(rdb);

  fact_db facts;
  create_input_fact(facts, "ci_life_src", 1);

  concrete_lifecycle_rule::events.clear();
  REQUIRE(makeQuery(rdb, facts, "ci_life_tgt"));

  REQUIRE(concrete_lifecycle_rule::events.size() == 3);
  CHECK(concrete_lifecycle_rule::events[0] == "prelude:3");
  CHECK(concrete_lifecycle_rule::events[1] == "compute:3");
  CHECK(concrete_lifecycle_rule::events[2] == "postlude:3");

  store<int> target(facts.get_variable("ci_life_tgt"));
  CHECK(target[2] == 4);
  CHECK(target[3] == 6);
  CHECK(target[4] == 8);
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
