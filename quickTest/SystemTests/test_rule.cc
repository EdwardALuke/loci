#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <map>
#include <string>

using namespace Loci;

namespace {

  /// @brief Wraps a concrete rule type in the external `rule` handle used by
  ///   the public APIs under test.
  template <class TRule>
  rule make_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  /// @brief Returns true when `desc` contains `v` in one of its variable sets.
  bool desc_has_var(const std::set<vmap_info> &desc, const variable &v) {
    for(std::set<vmap_info>::const_iterator i = desc.begin(); i != desc.end(); ++i)
      if(i->var.inSet(v))
        return true;
    return false;
  }

  /// @brief Returns true when `desc` contains `v` in one of its mapping sets.
  bool desc_has_map_var(const std::set<vmap_info> &desc, const variable &v) {
    for(std::set<vmap_info>::const_iterator i = desc.begin(); i != desc.end(); ++i)
      for(size_t j = 0; j < i->mapping.size(); ++j)
        if(i->mapping[j].inSet(v))
          return true;
    return false;
  }

  /// @brief Baseline pointwise rule used by the general add/remove/index tests.
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

  /// @brief Pointwise rule that carries mappings, a constraint, and a
  ///   conditional for descriptor, rename, and time-shift coverage.
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
// Additional behavior/invariant tests
//----------------------------------------------------------------------------

TEST_CASE("add_rule indexes target for simple pointwise rule") {
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();
  const variable tgt("tgt_idx");

  // Register the rule in the database.
  rdb.add_rule(r);

  // The reverse target index should now return that rule.
  CHECK(rdb.rules_by_target(tgt).inSet(r));
}

TEST_CASE("remove_rule clears target index for simple pointwise rule") {
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();
  const variable tgt("tgt_idx");

  // Add then remove the same rule through the public database API.
  rdb.add_rule(r);
  rdb.remove_rule(r);

  // The target index should be empty for that rule after removal.
  CHECK_FALSE(rdb.rules_by_target(tgt).inSet(r));
}

TEST_CASE("internal rule string parser captures qualifier source target constraint conditional") {
  rule r("qualifier(q_test),source(a_test),target(b_test),constraint(c_test),conditional(cond_test)");
  const rule::info &ri = r.get_info();

  // INTERNAL rules built from strings should preserve the parsed descriptor
  // structure.
  CHECK(r.type() == rule::INTERNAL);
  CHECK(ri.qualifier() == "q_test");
  CHECK(desc_has_var(ri.desc.sources, variable("a_test")));
  CHECK(desc_has_var(ri.desc.targets, variable("b_test")));
  CHECK(desc_has_var(ri.desc.constraints, variable("c_test")));
  CHECK(ri.desc.conditionals.inSet(variable("cond_test")));
}

TEST_CASE("internal rule rename_vars renames source target and conditional variables") {
  rule r("qualifier(q_rename),source(a_rn),target(b_rn),conditional(cond_rn)");
  std::map<variable, variable> rm;
  rm[variable("a_rn")] = variable("a_rn_new");
  rm[variable("b_rn")] = variable("b_rn_new");
  rm[variable("cond_rn")] = variable("cond_rn_new");

  rule renamed = r.rename_vars(rm);
  const rule::info &ri = renamed.get_info();

  // Source, target, and conditional names should all follow the rename map.
  CHECK(desc_has_var(ri.desc.sources, variable("a_rn_new")));
  CHECK_FALSE(desc_has_var(ri.desc.sources, variable("a_rn")));
  CHECK(desc_has_var(ri.desc.targets, variable("b_rn_new")));
  CHECK_FALSE(desc_has_var(ri.desc.targets, variable("b_rn")));
  CHECK(ri.desc.conditionals.inSet(variable("cond_rn_new")));
  CHECK_FALSE(ri.desc.conditionals.inSet(variable("cond_rn")));
}

TEST_CASE("promote_rule appends time level to variables") {
  rule r = make_rule<simple_rule_for_indexes>();
  time_ident n_level("n_prom", time_ident());
  rule promoted = promote_rule(r, n_level);

  // Promotion should append the time level to both source and target names.
  CHECK(promoted.sources().inSet(variable(variable("src_idx"), n_level)));
  CHECK(promoted.targets().inSet(variable(variable("tgt_idx"), n_level)));
  CHECK_FALSE(promoted.sources().inSet(variable("src_idx")));
  CHECK_FALSE(promoted.targets().inSet(variable("tgt_idx")));
}

TEST_CASE("prepend_rule prepends time level to variables") {
  rule r = make_rule<simple_rule_for_indexes>();
  time_ident n_level("n_pre", time_ident());
  rule prepended = prepend_rule(r, n_level);

  // Prepending should use the prefix-style time decoration on both ends.
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

  // Inputs should include data dependencies, constraints, conditionals, and
  // mapping variables referenced by the descriptor.
  CHECK(in.inSet(variable("src_mc")));
  CHECK(in.inSet(variable("c_mc")));
  CHECK(in.inSet(variable("cond_mc")));
  CHECK(in.inSet(variable("m_in")));
  CHECK(in.inSet(variable("m_out")));
  CHECK(desc_has_map_var(info.sources, variable("m_in")));
  CHECK(desc_has_map_var(info.targets, variable("m_out")));

  // Outputs should still include only the computed target variable.
  CHECK(out.inSet(variable("tgt_mc")));
  CHECK_FALSE(out.inSet(variable("src_mc")));
}

TEST_CASE("rule_identifier includes source and target variable names") {
  rule_implP impl = new copy_rule_impl<simple_rule_for_indexes>;
  const std::string rid = impl->get_info().rule_identifier();

  // The identifier should remain useful in diagnostics and logs.
  CHECK(rid.find("src_idx") != std::string::npos);
  CHECK(rid.find("tgt_idx") != std::string::npos);
  CHECK(rid.find("<-") != std::string::npos);
}

TEST_CASE("add_rule indexes deep priority chain targets at each priority level") {
  rule_db rdb;
  rule r = make_rule<deep_priority_target_rule>();

  // Adding the rule should populate every progressively de-prioritized target
  // form that lookups may ask for.
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

  // Seed the database and removal set with two independent rules.
  rdb.add_rule(r1);
  rdb.add_rule(r2);
  rs += r1;
  rs += r2;

  // Bulk removal should clear each rule from its target index.
  rdb.remove_rules(rs);

  CHECK_FALSE(rdb.rules_by_target(variable("tgt_idx")).inSet(r1));
  CHECK_FALSE(rdb.rules_by_target(variable("tgt_idx_b")).inSet(r2));
}

TEST_CASE("default and optional rules are routed to their dedicated rule sets") {
  rule_db rdb;
  rule dr = make_rule<default_category_rule>();
  rule orule = make_rule<optional_category_rule>();

  // Add one rule from each special category.
  rdb.add_rule(dr);
  rdb.add_rule(orule);

  // They should be routed into their dedicated collections instead of the
  // regular executable rule set.
  CHECK(rdb.get_default_rules().inSet(dr));
  CHECK(rdb.get_optional_rules().inSet(orule));
  CHECK_FALSE(rdb.all_rules().inSet(dr));
  CHECK_FALSE(rdb.all_rules().inSet(orule));
}

TEST_CASE("mapping and conditional variables appear in sources() for external rules") {
  rule r = make_rule<mapping_conditional_rule>();
  const variableSet &srcs = r.sources();

  // External source discovery should include every source-side dependency the
  // descriptor references.
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

  // The rename map should be applied consistently across mapping variables,
  // payload variables, constraints, and conditionals.
  CHECK(ri.sources().inSet(variable("m_in_new")));
  CHECK(ri.sources().inSet(variable("m_out_new")));
  CHECK(ri.sources().inSet(variable("src_mc_new")));
  CHECK(ri.sources().inSet(variable("c_mc_new")));
  CHECK(ri.sources().inSet(variable("cond_mc_new")));
  CHECK(ri.targets().inSet(variable("tgt_mc_new")));
  CHECK_FALSE(ri.targets().inSet(variable("tgt_mc")));
}

TEST_CASE("remove_rule followed by add_rule restores target index") {
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();

  // Round-trip the same rule through add, remove, and add again.
  rdb.add_rule(r);
  rdb.remove_rule(r);
  rdb.add_rule(r);

  // Re-adding should fully restore the target reverse index.
  CHECK(rdb.rules_by_target(variable("tgt_idx")).inSet(r));
}

TEST_CASE("replace_map_constraints converts map constraints to generated constraints") {
  rule_implP impl = new copy_rule_impl<map_constraint_replace_rule>;
  fact_db facts;

  // Seed one map constraint and one ordinary constraint.
  Map map_con;
  map_con.allocate(interval(1, 3));
  facts.create_fact("map_con", map_con);

  Constraint keep_con;
  *keep_con = interval(2, 4);
  facts.create_fact("keep_con", keep_con);

  // Replacing map constraints should rewrite only the map-backed constraint.
  impl->replace_map_constraints(facts);
  const rule_impl::info &info = impl->get_info();

  CHECK(desc_has_var(info.constraints, variable("__map_con_MAP_constraint__")));
  CHECK(desc_has_var(info.constraints, variable("keep_con")));
  CHECK_FALSE(desc_has_var(info.constraints, variable("map_con")));

  // The generated fact should be materialized as a constraint in `fact_db`.
  storeRepP generated = facts.get_variable("__map_con_MAP_constraint__");
  REQUIRE(generated != static_cast<storeRep *>(0));
  CHECK(isCONSTRAINT(generated));
}

TEST_CASE("replace_map_constraints is idempotent for generated constraint facts") {
  rule_implP impl = new copy_rule_impl<map_constraint_replace_rule>;
  fact_db facts;

  // Seed fresh map and ordinary constraints for the idempotence check.
  Map map_con;
  map_con.allocate(interval(7, 9));
  facts.create_fact("map_con", map_con);

  Constraint keep_con;
  *keep_con = interval(3, 5);
  facts.create_fact("keep_con", keep_con);

  // The first rewrite establishes the generated constraint fact and its domain.
  impl->replace_map_constraints(facts);
  storeRepP first = facts.get_variable("__map_con_MAP_constraint__");
  REQUIRE(first != static_cast<storeRep *>(0));
  const entitySet first_dom = first->domain();

  // A second rewrite should leave that generated fact stable.
  impl->replace_map_constraints(facts);
  storeRepP second = facts.get_variable("__map_con_MAP_constraint__");
  REQUIRE(second != static_cast<storeRep *>(0));
  CHECK(second->domain() == first_dom);
}

TEST_CASE("check_perm_bits type matrix catches bad and accepts good rule targets") {
  rule_implP good_default = new copy_rule_impl<default_category_rule>;
  rule_implP bad_default = new copy_rule_impl<default_bad_store_target_rule>;
  rule_implP bad_pointwise = new copy_rule_impl<pointwise_bad_param_target_rule>;
  rule_implP good_map = new copy_rule_impl<map_good_rule>;
  rule_implP bad_map = new copy_rule_impl<map_bad_store_target_rule>;
  rule_implP good_constraint = new copy_rule_impl<constraint_good_rule>;
  rule_implP bad_constraint = new copy_rule_impl<constraint_bad_store_target_rule>;

  // Each category should accept the container type it expects.
  CHECK(good_default->check_perm_bits());
  CHECK_FALSE(bad_default->check_perm_bits());
  CHECK_FALSE(bad_pointwise->check_perm_bits());
  CHECK(good_map->check_perm_bits());
  CHECK_FALSE(bad_map->check_perm_bits());
  CHECK(good_constraint->check_perm_bits());
  CHECK_FALSE(bad_constraint->check_perm_bits());
}

int main(int argc, char **argv) {
  Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Finalize();
  return result;
}
