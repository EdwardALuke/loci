#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <map>
#include <sstream>
#include <string>

using namespace Loci;

namespace {

  /// @brief Captures `std::cerr` so tests can assert on emitted diagnostics.
  class CErrCapture {
    std::streambuf *old_;
    std::ostringstream captured_;
  public:
    CErrCapture() : old_(std::cerr.rdbuf(captured_.rdbuf())) {}
    ~CErrCapture() { std::cerr.rdbuf(old_); }

    /// @brief Returns the text captured from `std::cerr`.
    std::string str() const { return captured_.str(); }
  };

  bool check_perm_bits_quietly(const rule_implP &impl) {
    CErrCapture capture;
    return impl->check_perm_bits();
  }

  std::string generated_map_constraint_name(const std::string &map_name) {
    return "__" + map_name + "_MAP_constraint__";
  }

  /// @brief Wraps a concrete rule type in the external `rule` handle used by
  ///   the public APIs under test.
  template <class TRule>
  rule make_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  /// @brief Returns true when low-level rule metadata contains `v`.
  bool desc_has_var(const std::set<vmap_info> &desc, const variable &v) {
    for(std::set<vmap_info>::const_iterator i = desc.begin(); i != desc.end(); ++i)
      if(i->var.inSet(v))
        return true;
    return false;
  }

  /// @brief Returns true when low-level rule metadata contains `v` as a map.
  bool desc_has_map_var(const std::set<vmap_info> &desc, const variable &v) {
    for(std::set<vmap_info>::const_iterator i = desc.begin(); i != desc.end(); ++i)
      for(size_t j = 0; j < i->mapping.size(); ++j)
        if(i->mapping[j].inSet(v))
          return true;
    return false;
  }

  /// @brief Pointwise rule whose priority-qualified target exercises
  ///   add/remove lookup across both qualified and base target names.
  class priority_target_rule_for_removal : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    priority_target_rule_for_removal() {
      name_store("src", src);
      name_store("p::tgt", tgt);
      input("src");
      output("p::tgt");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Baseline pointwise rule used by add/remove lookup tests.
  class simple_rule_for_indexes : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    simple_rule_for_indexes() {
      name_store("src_idx", src);
      name_store("tgt_idx", tgt);
      input("src_idx");
      output("tgt_idx");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Second independent pointwise rule used to test bulk removal paths.
  class simple_rule_for_indexes_b : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    simple_rule_for_indexes_b() {
      name_store("src_idx_b", src);
      name_store("tgt_idx_b", tgt);
      input("src_idx_b");
      output("tgt_idx_b");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Pointwise rule whose target carries a two-level priority chain.
  class deep_priority_target_rule : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    deep_priority_target_rule() {
      name_store("src_pri", src);
      name_store("p1::p2::tgt_pri", tgt);
      input("src_pri");
      output("p1::p2::tgt_pri");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Pointwise rule assigned to a non-default keyspace.
  class custom_keyspace_rule : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    custom_keyspace_rule() {
      keyspace_tag("ks_test");
      keyspace_dist_hint();
      name_store("ks_src", src);
      name_store("ks_tgt", tgt);
      input("ks_src");
      output("ks_tgt");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Valid `default_rule` fixture for permission-bit checks.
  class default_category_rule : public default_rule {
    param<int> out;
  public:
    default_category_rule() {
      name_store("default_out", out);
      output("default_out");
    }

    void compute(const sequence &) {}
  };

  /// @brief Invalid `default_rule` fixture with a store target for negative
  ///   permission checks.
  class default_bad_store_target_rule : public default_rule {
    store<int> out;
  public:
    default_bad_store_target_rule() {
      name_store("default_bad_store_out", out);
      output("default_bad_store_out");
    }

    void compute(const sequence &) {}
  };

  /// @brief Valid `optional_rule` fixture for rule-set routing checks.
  class optional_category_rule : public optional_rule {
    param<int> out;
  public:
    optional_category_rule() {
      name_store("optional_out", out);
      output("optional_out");
    }

    void compute(const sequence &) {}
  };

  /// @brief Invalid pointwise fixture with a parameter target for negative
  ///   permission checks.
  class pointwise_bad_param_target_rule : public pointwise_rule {
    const_store<int> src;
    param<int> tgt;
  public:
    pointwise_bad_param_target_rule() {
      name_store("src_bad_pw", src);
      name_store("tgt_bad_pw", tgt);
      input("src_bad_pw");
      output("tgt_bad_pw");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Valid `map_rule` fixture for permission-bit checks.
  class map_good_rule : public map_rule {
    const_store<int> src;
    Map out;
  public:
    map_good_rule() {
      name_store("src_map_good", src);
      name_store("map_good_out", out);
      input("src_map_good");
      output("map_good_out");
    }

    void compute(const sequence &) {}
  };

  /// @brief Invalid `map_rule` fixture with a store target for negative
  ///   permission checks.
  class map_bad_store_target_rule : public map_rule {
    const_store<int> src;
    store<int> out;
  public:
    map_bad_store_target_rule() {
      name_store("src_map_bad", src);
      name_store("map_bad_out", out);
      input("src_map_bad");
      output("map_bad_out");
    }

    void compute(const sequence &) {}
  };

  /// @brief Valid `constraint_rule` fixture for permission-bit checks.
  class constraint_good_rule : public constraint_rule {
    const_store<int> src;
    Constraint out;
  public:
    constraint_good_rule() {
      name_store("src_con_good", src);
      name_store("con_good_out", out);
      input("src_con_good");
      output("con_good_out");
    }

    void compute(const sequence &) {}
  };

  /// @brief Invalid `constraint_rule` fixture with a store target for negative
  ///   permission checks.
  class constraint_bad_store_target_rule : public constraint_rule {
    const_store<int> src;
    store<int> out;
  public:
    constraint_bad_store_target_rule() {
      name_store("src_con_bad", src);
      name_store("con_bad_out", out);
      input("src_con_bad");
      output("con_bad_out");
    }

    void compute(const sequence &) {}
  };

  /// @brief Pointwise rule that mixes map and non-map constraints for
  ///   `replace_map_constraints()` coverage.
  class map_constraint_replace_rule : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    map_constraint_replace_rule() {
      name_store("src_map_con", src);
      name_store("tgt_map_con", tgt);
      input("src_map_con");
      output("tgt_map_con");
      constraint("map_con,keep_con");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Pointwise rule with both static and dynamic-candidate constraints.
  class split_constraint_rule : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    split_constraint_rule() {
      name_store("src_split_con", src);
      name_store("tgt_split_con", tgt);
      input("src_split_con");
      output("tgt_split_con");
      constraint("static_con,dyn_con");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Singleton rule whose priority-qualified target exercises
  ///   priority-override validation.
  class singleton_priority_override_rule : public singleton_rule {
    const_param<int> in;
    param<int> out;
  public:
    singleton_priority_override_rule() {
      name_store("in", in);
      name_store("p::out", out);
      input("in");
      output("p::out");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Default rule whose target name is assigned to an unnamed alias.
  class unnamed_alias_target_rule : public default_rule {
    param<int> real_tgt;
  public:
    unnamed_alias_target_rule() {
      name_store("real_tgt", real_tgt);
      output("alias_tgt=real_tgt");
    }

    void compute(const sequence &) {}
  };

  /// @brief Pointwise rule whose alias target expression exercises target-type
  ///   mixing diagnostics in `rule::info` parsing.
  class assign_param_alias_rule : public pointwise_rule {
    const_param<int> in;
    param<int> x;
    param<int> y;
  public:
    assign_param_alias_rule() {
      name_store("in", in);
      name_store("x", x);
      name_store("y", y);
      input("in");
      output("x=y");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// @brief Pointwise rule that carries mappings, a constraint, and a
  ///   conditional for dependency, rename, and time-shift coverage.
  class mapping_conditional_rule : public pointwise_rule {
    const_Map in_map;
    const_Map out_map;
    const_store<int> src;
    const_store<int> c;
    const_param<bool> cond;
    store<int> tgt;
  public:
    mapping_conditional_rule() {
      name_store("m_in", in_map);
      name_store("m_out", out_map);
      name_store("src_mc", src);
      name_store("c_mc", c);
      name_store("cond_mc", cond);
      name_store("tgt_mc", tgt);
      input("m_in->src_mc");
      output("m_out->tgt_mc");
      constraint("c_mc");
      conditional("cond_mc");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

}

//----------------------------------------------------------------------------
// Regression tests for bugs that were visible through normal rule construction
// or scheduler lookup.  Prefer the public `rule` and `rule_db` views here; direct
// `rule_impl` probes are used only for validation helpers that have no higher
// level API.
//----------------------------------------------------------------------------

TEST_CASE("remove_rule removes both priority and base target indexes") {
  rule_db rdb;
  rule r = make_rule<priority_target_rule_for_removal>();

  const variable prio_target("p::tgt");
  const variable base_target("tgt");

  // The rule database answers the scheduler question "which rule can derive
  // this fact?" for both the priority-qualified override and the base fact.
  rdb.add_rule(r);
  CHECK(rdb.rules_by_target(prio_target).inSet(r));
  CHECK(rdb.rules_by_target(base_target).inSet(r));

  // Removing the rule should clear both scheduler-visible target forms.
  rdb.remove_rule(r);
  CHECK_FALSE(rdb.rules_by_target(prio_target).inSet(r));
  CHECK_FALSE(rdb.rules_by_target(base_target).inSet(r));
}

TEST_CASE("internal prepend preserves mapping variables") {
  rule internal_rule("qualifier(internal_map_test),source(m->src),target(tgt)");

  time_ident n_level("n", time_ident());
  rule prepended(n_level, internal_rule);

  const rule::info &info = prepended.get_info();
  const variable expected_mapping_var(n_level, variable("m"));
  const variable wrong_mapping_var(n_level, variable("src"));

  // For `m->src`, the map remains its own scheduler dependency when a time
  // label is prepended; the payload variable must not be mistaken for the map.
  CHECK(info.maps().inSet(expected_mapping_var));
  CHECK_FALSE(info.maps().inSet(wrong_mapping_var));
}

TEST_CASE("singleton priority override makes check_perm_bits fail") {
  rule_implP impl = new copy_rule_impl<singleton_priority_override_rule>;

  // Priority-qualified singleton targets would look like override facts, but
  // singletons produce shared parameters rather than entity-local facts.  The
  // validation path should reject that authoring mistake.
  CHECK_FALSE(check_perm_bits_quietly(impl));
}

TEST_CASE("parameter alias targets do not imply mixed parameter/store outputs") {
  CErrCapture capture;
  rule r = make_rule<assign_param_alias_rule>();

  // A target alias such as `x=y` should be typed by its backing variable, so a
  // pure parameter alias should not look like a mixed parameter/store target.
  std::string err = capture.str();
  CHECK(err.find("can't mix parameters and stores in target") == std::string::npos);

  (void)r;
}

TEST_CASE("get_store resolves unnamed alias targets through their backing variable") {
  rule_implP impl = new copy_rule_impl<unnamed_alias_target_rule>;

  storeRepP alias = impl->get_store(variable("alias_tgt"));
  storeRepP backing = impl->get_store(variable("real_tgt"));

  REQUIRE(alias != static_cast<storeRep *>(0));
  REQUIRE(backing != static_cast<storeRep *>(0));
  CHECK(&*alias == &*backing);
}

//----------------------------------------------------------------------------
// Additional behavior tests.  These favor API-observable rule behavior over
// exact implementation layout.
//----------------------------------------------------------------------------

TEST_CASE("internal rule string parser captures qualifier source target constraint conditional") {
  rule r("qualifier(q_test),source(a_test),target(b_test),constraint(c_test),conditional(cond_test)");
  const rule::info &ri = r.get_info();

  // INTERNAL rules model scheduler-created links.  The parsed rule should expose
  // the same sources, targets, constraints, and guards through `rule::info`.
  CHECK(r.type() == rule::INTERNAL);
  CHECK(ri.qualifier() == "q_test");
  CHECK(ri.sources().inSet(variable("a_test")));
  CHECK(ri.sources().inSet(variable("cond_test")));
  CHECK(ri.targets().inSet(variable("b_test")));
  CHECK(ri.constraints().inSet(variable("c_test")));
}

TEST_CASE("internal rule rename_vars renames dependencies and targets") {
  rule r("qualifier(q_rename),source(a_rn),target(b_rn),conditional(cond_rn)");
  std::map<variable, variable> rm;
  rm[variable("a_rn")] = variable("a_rn_new");
  rm[variable("b_rn")] = variable("b_rn_new");
  rm[variable("cond_rn")] = variable("cond_rn_new");

  rule renamed = r.rename_vars(rm);
  const rule::info &ri = renamed.get_info();

  // Renaming is the public mechanism behind promotion and internal scheduler
  // rewrites, so every visible dependency and target should follow the map.
  CHECK(ri.sources().inSet(variable("a_rn_new")));
  CHECK_FALSE(ri.sources().inSet(variable("a_rn")));
  CHECK(ri.sources().inSet(variable("cond_rn_new")));
  CHECK_FALSE(ri.sources().inSet(variable("cond_rn")));
  CHECK(ri.targets().inSet(variable("b_rn_new")));
  CHECK_FALSE(ri.targets().inSet(variable("b_rn")));
}

TEST_CASE("promote_rule appends time level to variables") {
  rule r = make_rule<simple_rule_for_indexes>();
  time_ident n_level("n_prom", time_ident());
  rule promoted = promote_rule(r, n_level);

  // Promotion makes a stationary rule usable at a time level: conceptually
  // `tgt <- src` becomes `tgt{n} <- src{n}`.
  CHECK(promoted.sources().inSet(variable(variable("src_idx"), n_level)));
  CHECK(promoted.targets().inSet(variable(variable("tgt_idx"), n_level)));
  CHECK_FALSE(promoted.sources().inSet(variable("src_idx")));
  CHECK_FALSE(promoted.targets().inSet(variable("tgt_idx")));
}

TEST_CASE("prepend_rule prepends time level to variables") {
  rule r = make_rule<simple_rule_for_indexes>();
  time_ident n_level("n_pre", time_ident());
  rule prepended = prepend_rule(r, n_level);

  // Internal scheduler rewrites use prefix-style labels; both the dependency
  // and derived fact should move together.
  CHECK(prepended.sources().inSet(variable(n_level, variable("src_idx"))));
  CHECK(prepended.targets().inSet(variable(n_level, variable("tgt_idx"))));
  CHECK_FALSE(prepended.sources().inSet(variable("src_idx")));
  CHECK_FALSE(prepended.targets().inSet(variable("tgt_idx")));
}

TEST_CASE("input_vars and output_vars contain expected variables for mapping and conditionals") {
  rule_implP impl = new copy_rule_impl<mapping_conditional_rule>;
  const rule_impl::info &info = impl->get_info();
  variableSet in = info.input_vars();
  variableSet out = info.output_vars();

  // The dependency set used by scheduler construction must include ordinary
  // inputs, constraints, conditionals, and maps crossed by mapped signatures.
  CHECK(in.inSet(variable("src_mc")));
  CHECK(in.inSet(variable("c_mc")));
  CHECK(in.inSet(variable("cond_mc")));
  CHECK(in.inSet(variable("m_in")));
  CHECK(in.inSet(variable("m_out")));
  CHECK(desc_has_map_var(info.sources, variable("m_in")));
  CHECK(desc_has_map_var(info.targets, variable("m_out")));

  // Mapped target variables are dependencies; the derived fact is still just
  // the target payload.
  CHECK(out.inSet(variable("tgt_mc")));
  CHECK_FALSE(out.inSet(variable("src_mc")));
}

TEST_CASE("add_rule indexes deep priority chain targets at each priority level") {
  rule_db rdb;
  rule r = make_rule<deep_priority_target_rule>();

  // A deeper priority-qualified fact should be discoverable as itself and
  // through each less-specific override layer.
  rdb.add_rule(r);
  CHECK(rdb.rules_by_target(variable("p1::p2::tgt_pri")).inSet(r));
  CHECK(rdb.rules_by_target(variable("p2::tgt_pri")).inSet(r));
  CHECK(rdb.rules_by_target(variable("tgt_pri")).inSet(r));
}

TEST_CASE("remove_rules clears target indexes for every rule in the input ruleSet") {
  rule_db rdb;
  rule r1 = make_rule<simple_rule_for_indexes>();
  rule r2 = make_rule<simple_rule_for_indexes_b>();
  ruleSet rs;

  // Seed the database and removal set with two independent derived facts.
  rdb.add_rule(r1);
  rdb.add_rule(r2);
  rs += r1;
  rs += r2;

  // Bulk removal should make both facts unavailable to scheduler lookup.
  rdb.remove_rules(rs);

  CHECK_FALSE(rdb.rules_by_target(variable("tgt_idx")).inSet(r1));
  CHECK_FALSE(rdb.rules_by_target(variable("tgt_idx_b")).inSet(r2));
}

TEST_CASE("add_rule and remove_rule maintain source indexes") {
  rule_db rdb;
  rule r = make_rule<mapping_conditional_rule>();

  // Source lookup is how the scheduler finds rules affected by a newly available
  // dependency, including mapped dependencies and conditionals.
  rdb.add_rule(r);
  CHECK(rdb.rules_by_source(variable("src_mc")).inSet(r));
  CHECK(rdb.rules_by_source(variable("c_mc")).inSet(r));
  CHECK(rdb.rules_by_source(variable("cond_mc")).inSet(r));
  CHECK(rdb.rules_by_source(variable("m_in")).inSet(r));
  CHECK(rdb.rules_by_source(variable("m_out")).inSet(r));

  rdb.remove_rule(r);
  CHECK_FALSE(rdb.rules_by_source(variable("src_mc")).inSet(r));
  CHECK_FALSE(rdb.rules_by_source(variable("c_mc")).inSet(r));
  CHECK_FALSE(rdb.rules_by_source(variable("cond_mc")).inSet(r));
  CHECK_FALSE(rdb.rules_by_source(variable("m_in")).inSet(r));
  CHECK_FALSE(rdb.rules_by_source(variable("m_out")).inSet(r));
}

TEST_CASE("default and optional rules are routed to their dedicated rule sets") {
  rule_db rdb;
  rule dr = make_rule<default_category_rule>();
  rule orule = make_rule<optional_category_rule>();

  // Add one rule from each special category.
  rdb.add_rule(dr);
  rdb.add_rule(orule);

  // Defaults and optional inputs seed user-facing facts; they should not be
  // scheduled as ordinary executable rules.
  CHECK(rdb.get_default_rules().inSet(dr));
  CHECK(rdb.get_optional_rules().inSet(orule));
  CHECK_FALSE(rdb.all_rules().inSet(dr));
  CHECK_FALSE(rdb.all_rules().inSet(orule));

  rdb.remove_rule(dr);
  rdb.remove_rule(orule);

  CHECK_FALSE(rdb.get_default_rules().inSet(dr));
  CHECK_FALSE(rdb.get_optional_rules().inSet(orule));
}

TEST_CASE("rule_db partitions non-default rules by keyspace") {
  rule_db rdb;
  rule r = make_rule<custom_keyspace_rule>();

  // Keyspace tagging lets the scheduler ask for rules in one partition without
  // also seeing rules from the main keyspace.
  rdb.add_rule(r);
  CHECK(rdb.all_rules().inSet(r));
  CHECK(rdb.all_rules("ks_test").inSet(r));
  CHECK_FALSE(rdb.all_rules("main").inSet(r));

  rdb.remove_rule(r);
  CHECK_FALSE(rdb.all_rules("ks_test").inSet(r));
}

TEST_CASE("mapping and conditional variables appear in sources() for external rules") {
  rule r = make_rule<mapping_conditional_rule>();
  const variableSet &srcs = r.sources();

  // External source discovery should include every fact that must exist before
  // the mapped rule can participate in a schedule.
  CHECK(srcs.inSet(variable("m_in")));
  CHECK(srcs.inSet(variable("m_out")));
  CHECK(srcs.inSet(variable("src_mc")));
  CHECK(srcs.inSet(variable("c_mc")));
  CHECK(srcs.inSet(variable("cond_mc")));
}

TEST_CASE("rename_vars on external rules renames mappings targets constraints and conditionals") {
  rule r = make_rule<mapping_conditional_rule>();
  std::map<variable, variable> rm;
  rm[variable("m_in")] = variable("m_in_new");
  rm[variable("m_out")] = variable("m_out_new");
  rm[variable("src_mc")] = variable("src_mc_new");
  rm[variable("tgt_mc")] = variable("tgt_mc_new");
  rm[variable("c_mc")] = variable("c_mc_new");
  rm[variable("cond_mc")] = variable("cond_mc_new");

  rule renamed = r.rename_vars(rm);
  const rule::info &ri = renamed.get_info();

  // Renaming should preserve the logical signature: maps, payload facts,
  // constraints, and conditionals all move together.
  CHECK(ri.sources().inSet(variable("m_in_new")));
  CHECK(ri.sources().inSet(variable("m_out_new")));
  CHECK(ri.sources().inSet(variable("src_mc_new")));
  CHECK(ri.sources().inSet(variable("c_mc_new")));
  CHECK(ri.sources().inSet(variable("cond_mc_new")));
  CHECK(ri.targets().inSet(variable("tgt_mc_new")));
  CHECK_FALSE(ri.targets().inSet(variable("tgt_mc")));
}

TEST_CASE("add_namespace namespaces mappings payload variables constraints and conditionals") {
  rule namespaced = make_rule<mapping_conditional_rule>().add_namespace("ns");
  const rule::info &ri = namespaced.get_info();

  CHECK(ri.sources().inSet(variable("src_mc").add_namespace("ns")));
  CHECK(ri.sources().inSet(variable("c_mc").add_namespace("ns")));
  CHECK(ri.sources().inSet(variable("cond_mc").add_namespace("ns")));
  CHECK(ri.sources().inSet(variable("m_in").add_namespace("ns")));
  CHECK(ri.sources().inSet(variable("m_out").add_namespace("ns")));
  CHECK(ri.targets().inSet(variable("tgt_mc").add_namespace("ns")));
  CHECK(ri.maps().inSet(variable("m_in").add_namespace("ns")));
  CHECK(ri.maps().inSet(variable("m_out").add_namespace("ns")));

  CHECK_FALSE(ri.sources().inSet(variable("src_mc")));
  CHECK_FALSE(ri.targets().inSet(variable("tgt_mc")));
}

TEST_CASE("remove_rule followed by add_rule restores target index") {
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();

  // Round-trip the same rule through add, remove, and add again.
  rdb.add_rule(r);
  CHECK(rdb.rules_by_target(variable("tgt_idx")).inSet(r));
  rdb.remove_rule(r);
  CHECK_FALSE(rdb.rules_by_target(variable("tgt_idx")).inSet(r));
  rdb.add_rule(r);

  // Re-adding should fully restore scheduler lookup for the derived fact.
  CHECK(rdb.rules_by_target(variable("tgt_idx")).inSet(r));
}

TEST_CASE("replace_map_constraints converts map constraints to generated constraints") {
  rule_implP impl = new copy_rule_impl<map_constraint_replace_rule>;
  fact_db facts;
  const std::string generated_name = generated_map_constraint_name("map_con");

  // Seed one map-backed legality set and one ordinary constraint.
  Map map_con;
  map_con.allocate(interval(1, 3));
  facts.create_fact("map_con", map_con);

  Constraint keep_con;
  *keep_con = interval(2, 4);
  facts.create_fact("keep_con", keep_con);

  // The scheduler needs constraints to be set-valued facts.  A map listed as a
  // constraint should be materialized as a generated constraint over its domain.
  impl->replace_map_constraints(facts);
  rule r(impl);
  const rule::info &info = r.get_info();

  CHECK(info.constraints().inSet(variable(generated_name)));
  CHECK(info.constraints().inSet(variable("keep_con")));
  CHECK_FALSE(info.constraints().inSet(variable("map_con")));

  // The generated fact should be materialized as a constraint in `fact_db`.
  storeRepP generated = facts.get_variable(generated_name);
  REQUIRE(generated != static_cast<storeRep *>(0));
  CHECK(isCONSTRAINT(generated));
}

TEST_CASE("replace_map_constraints is idempotent for generated constraint facts") {
  rule_implP impl = new copy_rule_impl<map_constraint_replace_rule>;
  fact_db facts;
  const std::string generated_name = generated_map_constraint_name("map_con");

  // Seed fresh map-backed and ordinary constraints for the idempotence check.
  Map map_con;
  map_con.allocate(interval(7, 9));
  facts.create_fact("map_con", map_con);

  Constraint keep_con;
  *keep_con = interval(3, 5);
  facts.create_fact("keep_con", keep_con);

  // The first rewrite establishes the generated constraint fact and its domain.
  impl->replace_map_constraints(facts);
  storeRepP first = facts.get_variable(generated_name);
  REQUIRE(first != static_cast<storeRep *>(0));
  const entitySet first_dom = first->domain();

  // A second rewrite should leave that generated fact stable, so repeated
  // scheduler preparation does not keep changing the fact database.
  impl->replace_map_constraints(facts);
  storeRepP second = facts.get_variable(generated_name);
  REQUIRE(second != static_cast<storeRep *>(0));
  CHECK(second->domain() == first_dom);
}

TEST_CASE("split_constraints moves selected constraints to the dynamic set") {
  rule_implP impl = new copy_rule_impl<split_constraint_rule>;
  variableSet dynamic_constraints;
  dynamic_constraints += variable("dyn_con");

  impl->split_constraints(dynamic_constraints);
  const rule_impl::info &info = impl->get_info();

  // Dynamic constraints are still legality sets, but they are held separately
  // so later schedule construction can re-evaluate them.
  CHECK(desc_has_var(info.constraints, variable("static_con")));
  CHECK_FALSE(desc_has_var(info.constraints, variable("dyn_con")));
  CHECK(desc_has_var(info.dynamic_constraints, variable("dyn_con")));
}

TEST_CASE("check_perm_bits type matrix catches bad and accepts good rule targets") {
  rule_implP good_default = new copy_rule_impl<default_category_rule>;
  rule_implP bad_default = new copy_rule_impl<default_bad_store_target_rule>;
  rule_implP bad_pointwise = new copy_rule_impl<pointwise_bad_param_target_rule>;
  rule_implP good_map = new copy_rule_impl<map_good_rule>;
  rule_implP bad_map = new copy_rule_impl<map_bad_store_target_rule>;
  rule_implP good_constraint = new copy_rule_impl<constraint_good_rule>;
  rule_implP bad_constraint = new copy_rule_impl<constraint_bad_store_target_rule>;

  // These are rule-authoring contracts: pointwise rules derive stores, defaults
  // derive parameters, map rules derive maps, and constraint rules derive
  // activation sets.
  CHECK(check_perm_bits_quietly(good_default));
  CHECK_FALSE(check_perm_bits_quietly(bad_default));
  CHECK_FALSE(check_perm_bits_quietly(bad_pointwise));
  CHECK(check_perm_bits_quietly(good_map));
  CHECK_FALSE(check_perm_bits_quietly(bad_map));
  CHECK(check_perm_bits_quietly(good_constraint));
  CHECK_FALSE(check_perm_bits_quietly(bad_constraint));
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
