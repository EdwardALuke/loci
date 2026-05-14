#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

namespace {

  /*
    These tests intentionally reach dependency_graph2 only through makeQuery.
    That keeps the assertions on scheduler-visible behavior instead of on the
    exact graph representation, which is free to change as scheduling evolves.
  */

  template <class TRule>
  rule add_rule(rule_db &rdb) {
    const rule r(rule_implP(new copy_rule_impl<TRule>));
    rdb.add_rule(r);
    return r;
  }

  void fill_store(store<int> &values, const entitySet &domain, int offset) {
    values.allocate(domain);
    for(entitySet::const_iterator i = domain.begin(); i != domain.end(); ++i)
      values[*i] = offset + *i;
  }

  void check_store_values(const store<int> &values, const entitySet &domain,
                          int offset) {
    CHECK(values.domain() == domain);
    for(entitySet::const_iterator i = domain.begin(); i != domain.end(); ++i) {
      CAPTURE(*i);
      CHECK(values[*i] == offset + *i);
    }
  }

  int source_to_mid_computes = 0;
  int mid_to_target_computes = 0;
  int priority_target_computes = 0;

  void reset_counts() {
    source_to_mid_computes = 0;
    mid_to_target_computes = 0;
    priority_target_computes = 0;
  }

  class source_to_mid_rule : public pointwise_rule {
    const_store<int> source;
    store<int> mid;
  public:
    source_to_mid_rule() {
      name_store("dg_source", source);
      name_store("dg_mid", mid);
      input("dg_source");
      output("dg_mid");
    }

    void calculate(Entity e) {
      mid[e] = source[e] + 1;
    }

    void compute(const sequence &seq) {
      ++source_to_mid_computes;
      do_loop(seq, this);
    }
  };

  class mid_to_target_rule : public pointwise_rule {
    const_store<int> mid;
    store<int> target;
  public:
    mid_to_target_rule() {
      name_store("dg_mid", mid);
      name_store("dg_target", target);
      input("dg_mid");
      output("dg_target");
    }

    void calculate(Entity e) {
      target[e] = mid[e] + 10;
    }

    void compute(const sequence &seq) {
      ++mid_to_target_computes;
      do_loop(seq, this);
    }
  };

  class priority_target_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;
  public:
    priority_target_rule() {
      name_store("dg_priority_source", source);
      name_store("p::dg_priority_target", target);
      input("dg_priority_source");
      output("p::dg_priority_target");
    }

    void calculate(Entity e) {
      target[e] = source[e] + 30;
    }

    void compute(const sequence &seq) {
      ++priority_target_computes;
      do_loop(seq, this);
    }
  };

} // namespace

TEST_CASE("makeQuery uses an extensional intermediate without recomputing it") {
  /*
    Value over the FVM regression: this isolates a scheduler bug where a rule
    that can produce an already-given fact is still scheduled and overwrites
    user-provided data. FVM usually does not seed intentionally conflicting
    extensional values, so it can miss this failure mode or make it hard to
    identify as dependency-graph pruning.
  */
  reset_counts();

  rule_db rdb;
  add_rule<source_to_mid_rule>(rdb);
  add_rule<mid_to_target_rule>(rdb);

  fact_db facts;
  const entitySet domain(interval(1, 3));
  store<int> source;
  store<int> mid;
  fill_store(source, domain, 1000);
  fill_store(mid, domain, 2000);
  facts.create_fact("dg_source", source);
  facts.create_fact("dg_mid", mid);

  CHECK(makeQuery(rdb, facts, "dg_target"));

  store<int> target(facts.get_variable("dg_target"));
  check_store_values(target, domain, 2010);
  CHECK(source_to_mid_computes == 0);
  CHECK(mid_to_target_computes == 1);
}

TEST_CASE("makeQuery satisfies an unqualified target from a priority output") {
  /*
    Value over the FVM regression: this is a minimal public-query check for the
    priority bridge from p::x to x. A larger FVM case may exercise priority
    indirectly, but this catches the specific regression where the scheduler
    builds a producer for the qualified fact yet fails to satisfy the user's
    unqualified query.
  */
  reset_counts();

  rule_db rdb;
  add_rule<priority_target_rule>(rdb);

  fact_db facts;
  const entitySet domain(interval(5, 7));
  store<int> source;
  fill_store(source, domain, 3000);
  facts.create_fact("dg_priority_source", source);

  CHECK(makeQuery(rdb, facts, "dg_priority_target"));

  store<int> target(facts.get_variable("dg_priority_target"));
  check_store_values(target, domain, 3030);
  CHECK(priority_target_computes == 1);
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
