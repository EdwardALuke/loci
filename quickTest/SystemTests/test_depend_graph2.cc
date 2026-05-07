#include <Loci.h>
#include <depend_graph.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <sstream>

using namespace Loci;

namespace {

  /// @brief Captures `std::cerr` so negative-path tests can assert diagnostics.
  class CErrCapture {
    std::streambuf *old_;
    std::ostringstream captured_;
  public:
    CErrCapture() : old_(std::cerr.rdbuf(captured_.rdbuf())) {}
    ~CErrCapture() { std::cerr.rdbuf(old_); }

    std::string str() const { return captured_.str(); }
  };

  /// @brief Wraps a concrete rule fixture in the public `rule` handle.
  template <class TRule>
  rule make_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  /// @brief Adds a rule fixture to `rdb` and returns its canonical handle.
  template <class TRule>
  rule add_rule(rule_db &rdb) {
    const rule r = make_rule<TRule>();
    rdb.add_rule(r);
    return r;
  }

  /// @brief Pointwise rule used as the first edge in dependency-chain tests.
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

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Pointwise rule used as the final edge in dependency-chain tests.
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

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Independent pointwise rule whose target exercises priority edges.
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

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Returns a single-variable set for compact graph setup.
  variableSet vars(const char *name) {
    variableSet result;
    result += variable(name);
    return result;
  }

} // namespace

TEST_CASE("dependency_graph2 builds source-rule-target edges for a rule chain") {
  rule_db rdb;
  const rule first = add_rule<source_to_mid_rule>(rdb);
  const rule second = add_rule<mid_to_target_rule>(rdb);

  const digraph graph =
    dependency_graph2(rdb, vars("dg_source"), vars("dg_target")).get_graph();

  CHECK(graph.is_edge(variable("dg_source").ident(), first.ident()));
  CHECK(graph.is_edge(first.ident(), variable("dg_mid").ident()));
  CHECK(graph.is_edge(variable("dg_mid").ident(), second.ident()));
  CHECK(graph.is_edge(second.ident(), variable("dg_target").ident()));
}

TEST_CASE("dependency_graph2 returns an empty graph when no rule can produce the target") {
  rule_db rdb;
  add_rule<source_to_mid_rule>(rdb);

  CErrCapture err;
  const digraph graph =
    dependency_graph2(rdb, vars("dg_source"), vars("dg_missing")).get_graph();

  CHECK(graph.get_all_vertices() == EMPTY);
  CHECK(graph.get_target_vertices() == EMPTY);
  CHECK(err.str().find("dg_missing cannot be inferred") != std::string::npos);
}

TEST_CASE("dependency_graph2 prunes rules whose targets are already given") {
  rule_db rdb;
  const rule pruned = add_rule<source_to_mid_rule>(rdb);
  const rule needed = add_rule<mid_to_target_rule>(rdb);

  const digraph graph =
    dependency_graph2(rdb, vars("dg_mid"), vars("dg_target")).get_graph();

  CHECK_FALSE(graph.get_all_vertices().inSet(pruned.ident()));
  CHECK_FALSE(graph.is_edge(variable("dg_source").ident(), pruned.ident()));
  CHECK(graph.is_edge(variable("dg_mid").ident(), needed.ident()));
  CHECK(graph.is_edge(needed.ident(), variable("dg_target").ident()));
}

TEST_CASE("dependency_graph2 inserts priority edges from qualified target to base target") {
  rule_db rdb;
  const rule producer = add_rule<priority_target_rule>(rdb);
  const rule priority("source(p::dg_priority_target),target(dg_priority_target),qualifier(priority)");

  const digraph graph =
    dependency_graph2(rdb, vars("dg_priority_source"), vars("dg_priority_target")).get_graph();

  CHECK(graph.is_edge(variable("dg_priority_source").ident(), producer.ident()));
  CHECK(graph.is_edge(producer.ident(), variable("p::dg_priority_target").ident()));
  CHECK(graph.is_edge(variable("p::dg_priority_target").ident(), priority.ident()));
  CHECK(graph.is_edge(priority.ident(), variable("dg_priority_target").ident()));
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
