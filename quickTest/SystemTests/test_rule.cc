#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <map>
#include <sstream>
#include <string>

using namespace Loci;

namespace {

  // Captures stderr so tests can assert on diagnostics emitted by rule parsing
  // and permission validation paths.
  class CErrCapture {
    std::streambuf *old_;
    std::ostringstream captured_;
  public:
    CErrCapture() : old_(std::cerr.rdbuf(captured_.rdbuf())) {}
    ~CErrCapture() { std::cerr.rdbuf(old_); }

    std::string str() const { return captured_.str(); }
  };

  // Builds an external rule wrapper from a concrete rule type so tests can use
  // rule_db/rule APIs without hand-writing boilerplate constructors.
  template <class TRule>
  rule make_rule() {
    return rule(rule_implP(new copy_rule_impl<TRule>));
  }

  // Checks whether a rule_impl descriptor contains a specific target variable.
  // Used by name/indexing tests that fetch cloned rule_impl objects.
  bool impl_has_target(const rule_implP &impl, const variable &target_var) {
    const std::set<vmap_info> &targets = impl->get_info().targets;
    for(std::set<vmap_info>::const_iterator i = targets.begin(); i != targets.end(); ++i)
      if(i->var.inSet(target_var))
        return true;
    return false;
  }

  // Lightweight ruleSet cardinality helper used by duplicate insertion tests.
  size_t count_rules(const ruleSet &rs) {
    size_t c = 0;
    for(ruleSet::const_iterator i = rs.begin(); i != rs.end(); ++i)
      ++c;
    return c;
  }

  // Descriptor membership helper for source/target/constraint sets.
  bool desc_has_var(const std::set<vmap_info> &desc, const variable &v) {
    for(std::set<vmap_info>::const_iterator i = desc.begin(); i != desc.end(); ++i)
      if(i->var.inSet(v))
        return true;
    return false;
  }

  // Descriptor membership helper for mapping variables inside vmap_info.
  bool desc_has_map_var(const std::set<vmap_info> &desc, const variable &v) {
    for(std::set<vmap_info>::const_iterator i = desc.begin(); i != desc.end(); ++i)
      for(size_t j = 0; j < i->mapping.size(); ++j)
        if(i->mapping[j].inSet(v))
          return true;
    return false;
  }

  // Uses a priority-qualified target so remove_rule target-index cleanup can be
  // checked across both qualified and base variable forms.
  class priority_target_rule_for_removal : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    priority_target_rule_for_removal() {
      rule_name("priority_target_rule_for_removal");
      name_store("src", src);
      name_store("p::tgt", tgt);
      input("src");
      output("p::tgt");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Baseline pointwise rule used by most add/remove/index/rename/timing tests.
  class simple_rule_for_indexes : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    simple_rule_for_indexes() {
      rule_name("simple_rule_for_indexes");
      name_store("src_idx", src);
      name_store("tgt_idx", tgt);
      input("src_idx");
      output("tgt_idx");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Second independent baseline rule to validate multi-rule remove_rules paths.
  class simple_rule_for_indexes_b : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    simple_rule_for_indexes_b() {
      rule_name("simple_rule_for_indexes_b");
      name_store("src_idx_b", src);
      name_store("tgt_idx_b", tgt);
      input("src_idx_b");
      output("tgt_idx_b");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Deep priority chain target used to validate indexing at each priority level.
  class deep_priority_target_rule : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    deep_priority_target_rule() {
      rule_name("deep_priority_target_rule");
      name_store("src_pri", src);
      name_store("p1::p2::tgt_pri", tgt);
      input("src_pri");
      output("p1::p2::tgt_pri");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Valid default_rule shape for check_perm_bits type-matrix coverage.
  class default_category_rule : public default_rule {
    param<int> out;
  public:
    default_category_rule() {
      rule_name("default_category_rule");
      name_store("default_out", out);
      output("default_out");
    }

    void compute(const sequence &) {}
  };

  // Intentionally invalid default_rule target type for negative validation.
  class default_bad_store_target_rule : public default_rule {
    store<int> out;
  public:
    default_bad_store_target_rule() {
      rule_name("default_bad_store_target_rule");
      name_store("default_bad_store_out", out);
      output("default_bad_store_out");
    }

    void compute(const sequence &) {}
  };

  // Valid optional_rule shape for category-routing checks in rule_db.
  class optional_category_rule : public optional_rule {
    param<int> out;
  public:
    optional_category_rule() {
      rule_name("optional_category_rule");
      name_store("optional_out", out);
      output("optional_out");
    }

    void compute(const sequence &) {}
  };

  // Intentionally invalid pointwise target type for negative check_perm_bits.
  class pointwise_bad_param_target_rule : public pointwise_rule {
    const_store<int> src;
    param<int> tgt;
  public:
    pointwise_bad_param_target_rule() {
      rule_name("pointwise_bad_param_target_rule");
      name_store("src_bad_pw", src);
      name_store("tgt_bad_pw", tgt);
      input("src_bad_pw");
      output("tgt_bad_pw");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Valid map_rule shape for map target permission checks.
  class map_good_rule : public map_rule {
    const_store<int> src;
    Map out;
  public:
    map_good_rule() {
      rule_name("map_good_rule");
      name_store("src_map_good", src);
      name_store("map_good_out", out);
      input("src_map_good");
      output("map_good_out");
    }

    void compute(const sequence &) {}
  };

  // Intentionally invalid map_rule target type for negative validation.
  class map_bad_store_target_rule : public map_rule {
    const_store<int> src;
    store<int> out;
  public:
    map_bad_store_target_rule() {
      rule_name("map_bad_store_target_rule");
      name_store("src_map_bad", src);
      name_store("map_bad_out", out);
      input("src_map_bad");
      output("map_bad_out");
    }

    void compute(const sequence &) {}
  };

  // Valid constraint_rule shape for constraint target permission checks.
  class constraint_good_rule : public constraint_rule {
    const_store<int> src;
    Constraint out;
  public:
    constraint_good_rule() {
      rule_name("constraint_good_rule");
      name_store("src_con_good", src);
      name_store("con_good_out", out);
      input("src_con_good");
      output("con_good_out");
    }

    void compute(const sequence &) {}
  };

  // Intentionally invalid constraint_rule target type for negative validation.
  class constraint_bad_store_target_rule : public constraint_rule {
    const_store<int> src;
    store<int> out;
  public:
    constraint_bad_store_target_rule() {
      rule_name("constraint_bad_store_target_rule");
      name_store("src_con_bad", src);
      name_store("con_bad_out", out);
      input("src_con_bad");
      output("con_bad_out");
    }

    void compute(const sequence &) {}
  };

  // Carries map + non-map constraints so replace_map_constraints behavior can
  // be asserted.
  class map_constraint_replace_rule : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    map_constraint_replace_rule() {
      rule_name("map_constraint_replace_rule");
      name_store("src_map_con", src);
      name_store("tgt_map_con", tgt);
      input("src_map_con");
      output("tgt_map_con");
      constraint("map_con,keep_con");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Same public rule name as name_key_rule_v2 to exercise name2rule key reuse
  // and remove/add behavior.
  class name_key_rule_v1 : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    name_key_rule_v1() {
      rule_name("name_key_reuse_rule");
      name_store("src", src);
      name_store("tgt_a", tgt);
      input("src");
      output("tgt_a");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Same public rule name as name_key_rule_v1 with a different target, used to
  // detect stale name-key entries after remove/add.
  class name_key_rule_v2 : public pointwise_rule {
    const_store<int> src;
    store<int> tgt;
  public:
    name_key_rule_v2() {
      rule_name("name_key_reuse_rule");
      name_store("src", src);
      name_store("tgt_b", tgt);
      input("src");
      output("tgt_b");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Singleton rule with priority-qualified target; used to expose priority
  // override validation behavior in check_perm_bits.
  class singleton_priority_override_rule : public singleton_rule {
    const_param<int> in;
    param<int> out;
  public:
    singleton_priority_override_rule() {
      rule_name("singleton_priority_override_rule");
      name_store("in", in);
      name_store("p::out", out);
      input("in");
      output("p::out");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Alias target expression x=y is used to exercise parameter/store mixing
  // diagnostics in rule::info parsing.
  class assign_param_alias_rule : public pointwise_rule {
    const_param<int> in;
    param<int> x;
    param<int> y;
  public:
    assign_param_alias_rule() {
      rule_name("assign_param_alias_rule");
      name_store("in", in);
      name_store("x", x);
      name_store("y", y);
      input("in");
      output("x=y");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  // Includes source/target mappings, constraint, and conditional to cover
  // descriptor input_vars/output_vars, rename, and time-shift behavior.
  class mapping_conditional_rule : public pointwise_rule {
    const_Map in_map;
    const_Map out_map;
    const_store<int> src;
    const_store<int> c;
    const_param<bool> cond;
    store<int> tgt;
  public:
    mapping_conditional_rule() {
      rule_name("mapping_conditional_rule");
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

  // Two constraints so split_constraints can move selected variables to the
  // dynamic constraint descriptor set.
  class constraint_split_rule : public pointwise_rule {
    const_store<int> src;
    const_store<int> c1;
    const_store<int> c2;
    store<int> tgt;
  public:
    constraint_split_rule() {
      rule_name("constraint_split_rule");
      name_store("src_split", src);
      name_store("c_split1", c1);
      name_store("c_split2", c2);
      name_store("tgt_split", tgt);
      input("src_split");
      output("tgt_split");
      constraint("c_split1,c_split2");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

}

//----------------------------------------------------------------------------
// Regression Tests For Known Bugs
//----------------------------------------------------------------------------

// Why this fails:
// add_rule() indexes priority-qualified targets at each priority level, but
// remove_rule() mistakenly does "+=" instead of "-=" in the priority loop.
// That leaves stale entries in trgt2rule for priority forms after removal.
//
// Potential fix:
// In src/System/rule.cc (remove_rule), change trgt2rule[v] += f to
// trgt2rule[v] -= f inside the while(v.get_info().priority.size() != 0) loop.
TEST_CASE("remove_rule removes both priority and base target indexes") {
  rule_db rdb;
  rule r = make_rule<priority_target_rule_for_removal>();

  const variable prio_target("p::tgt");
  const variable base_target("tgt");

  rdb.add_rule(r);
  CHECK(rdb.rules_by_target(prio_target).inSet(r));
  CHECK(rdb.rules_by_target(base_target).inSet(r));

  rdb.remove_rule(r);
  CHECK_FALSE(rdb.rules_by_target(prio_target).inSet(r));
  CHECK_FALSE(rdb.rules_by_target(base_target).inSet(r));
}

// Why this fails:
// add_rule() stores rule by either get_name() or info.name() depending on
// collision, but remove_rule() uses inverted key-selection logic when erasing.
// This can erase the wrong key and leave a stale rule in name2rule.
//
// Potential fix:
// In src/System/rule.cc (remove_rule), mirror the key choice used in add_rule()
// so erase uses the exact key that was inserted for this rule.
TEST_CASE("remove_rule erases primary name key used by add_rule") {
  rule_db rdb;
  rule r1 = make_rule<name_key_rule_v1>();
  rule r2 = make_rule<name_key_rule_v2>();

  std::string shared_name("name_key_reuse_rule");

  rdb.add_rule(r1);
  rdb.remove_rule(r1);
  rdb.add_rule(r2);

  rule_implP by_name = rdb.get_rule(shared_name);

  CHECK(impl_has_target(by_name, variable("tgt_b")));
  CHECK_FALSE(impl_has_target(by_name, variable("tgt_a")));
}

// Why this fails:
// prepend_vmap_info() prepends mapping entries from in.var instead of from each
// actual mapping set entry (*i). This corrupts mapping chains when prepending
// time levels on INTERNAL rules.
//
// Potential fix:
// In src/System/rule.cc, prepend_vmap_info(), replace
//   res.mapping.push_back(prepend_set(in.var, tl))
// with
//   res.mapping.push_back(prepend_set(*i, tl))
TEST_CASE("internal prepend preserves mapping variables") {
  rule internal_rule("qualifier(internal_map_test),source(m->src),target(tgt)");

  time_ident n_level("n", time_ident());
  rule prepended(n_level, internal_rule);

  const std::set<vmap_info> &sources = prepended.get_info().desc.sources;
  REQUIRE(!sources.empty());

  const vmap_info &first_source = *(sources.begin());
  REQUIRE(!first_source.mapping.empty());

  const variable expected_mapping_var(n_level, variable("m"));
  const variable wrong_mapping_var(n_level, variable("src"));

  CHECK(first_source.mapping.front().inSet(expected_mapping_var));
  CHECK_FALSE(first_source.mapping.front().inSet(wrong_mapping_var));
}

// Why this fails:
// check_perm_bits() logs that singleton rules cannot use priority override, but
// does not set retval=false in that branch, so validation still reports success.
//
// Potential fix:
// In src/System/rule.cc singleton-check branch, set retval=false when priority
// override is detected.
TEST_CASE("singleton priority override makes check_perm_bits fail") {
  rule_implP impl = new copy_rule_impl<singleton_priority_override_rule>;
  CHECK_FALSE(impl->check_perm_bits());
}

// Why this fails:
// In rule::info construction, iterator 'i' (from tvar_types) is compared
// against tvars.begin(). Those iterators are from different containers, so
// first-element checks are wrong and may trigger false "mix parameter/store"
// diagnostics for alias targets like x=y.
//
// Potential fix:
// In src/System/rule.cc around the parameter/store mixing checks, compare
// against tvar_types.begin() (or use an explicit first-element flag) instead
// of comparing to tvars.begin().
TEST_CASE("parameter alias targets do not imply mixed parameter/store outputs") {
  CErrCapture capture;
  rule r = make_rule<assign_param_alias_rule>();

  std::string err = capture.str();
  CHECK(err.find("can't mix parameters and stores in target") == std::string::npos);

  (void)r;
}

//----------------------------------------------------------------------------
// Additional behavior/invariant tests
//----------------------------------------------------------------------------

// Verifies that adding a normal pointwise rule populates both reverse indexes:
// source->rules and target->rules. If either index misses the rule, dependency
// discovery by variable lookup becomes incomplete.
TEST_CASE("add_rule indexes source and target for simple pointwise rule") {
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();
  const variable src("src_idx");
  const variable tgt("tgt_idx");

  rdb.add_rule(r);

  CHECK(rdb.rules_by_source(src).inSet(r));
  CHECK(rdb.rules_by_target(tgt).inSet(r));
}

// Verifies remove_rule fully undoes index insertion for a simple external rule.
// After removal, querying by either source or target should not return the rule.
TEST_CASE("remove_rule clears source and target indexes for simple pointwise rule") {
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();
  const variable src("src_idx");
  const variable tgt("tgt_idx");

  rdb.add_rule(r);
  rdb.remove_rule(r);

  CHECK_FALSE(rdb.rules_by_source(src).inSet(r));
  CHECK_FALSE(rdb.rules_by_target(tgt).inSet(r));
}

// Verifies internal rule string parsing populates expected descriptor fields
// (qualifier, sources, targets, constraints, conditionals) and INTERNAL type.
TEST_CASE("internal rule string parser captures qualifier source target constraint conditional") {
  rule r("qualifier(q_test),source(a_test),target(b_test),constraint(c_test),conditional(cond_test)");
  const rule::info &ri = r.get_info();

  CHECK(r.type() == rule::INTERNAL);
  CHECK(ri.qualifier() == "q_test");
  CHECK(desc_has_var(ri.desc.sources, variable("a_test")));
  CHECK(desc_has_var(ri.desc.targets, variable("b_test")));
  CHECK(desc_has_var(ri.desc.constraints, variable("c_test")));
  CHECK(ri.desc.conditionals.inSet(variable("cond_test")));
}

// Verifies rename_vars() on INTERNAL rules rewrites descriptor variables and
// conditionals, and that old names are no longer present.
TEST_CASE("internal rule rename_vars renames source target and conditional variables") {
  rule r("qualifier(q_rename),source(a_rn),target(b_rn),conditional(cond_rn)");
  std::map<variable, variable> rm;
  rm[variable("a_rn")] = variable("a_rn_new");
  rm[variable("b_rn")] = variable("b_rn_new");
  rm[variable("cond_rn")] = variable("cond_rn_new");

  rule renamed = r.rename_vars(rm);
  const rule::info &ri = renamed.get_info();

  CHECK(desc_has_var(ri.desc.sources, variable("a_rn_new")));
  CHECK_FALSE(desc_has_var(ri.desc.sources, variable("a_rn")));
  CHECK(desc_has_var(ri.desc.targets, variable("b_rn_new")));
  CHECK_FALSE(desc_has_var(ri.desc.targets, variable("b_rn")));
  CHECK(ri.desc.conditionals.inSet(variable("cond_rn_new")));
  CHECK_FALSE(ri.desc.conditionals.inSet(variable("cond_rn")));
}

// Verifies promote_rule() appends the given time level to source/target vars
// (postfix time) and removes hits under unpromoted variable names.
TEST_CASE("promote_rule appends time level to variables") {
  rule r = make_rule<simple_rule_for_indexes>();
  time_ident n_level("n_prom", time_ident());
  rule promoted = promote_rule(r, n_level);

  CHECK(promoted.sources().inSet(variable(variable("src_idx"), n_level)));
  CHECK(promoted.targets().inSet(variable(variable("tgt_idx"), n_level)));
  CHECK_FALSE(promoted.sources().inSet(variable("src_idx")));
  CHECK_FALSE(promoted.targets().inSet(variable("tgt_idx")));
}

// Verifies prepend_rule() prepends the given time level to source/target vars
// (prefix time) and removes hits under unqualified variable names.
TEST_CASE("prepend_rule prepends time level to variables") {
  rule r = make_rule<simple_rule_for_indexes>();
  time_ident n_level("n_pre", time_ident());
  rule prepended = prepend_rule(r, n_level);

  CHECK(prepended.sources().inSet(variable(n_level, variable("src_idx"))));
  CHECK(prepended.targets().inSet(variable(n_level, variable("tgt_idx"))));
  CHECK_FALSE(prepended.sources().inSet(variable("src_idx")));
  CHECK_FALSE(prepended.targets().inSet(variable("tgt_idx")));
}

// Verifies add_namespace() consistently rewrites source/target variables to
// namespace-qualified names.
TEST_CASE("add_namespace applies namespace to source and target variables") {
  rule r = make_rule<simple_rule_for_indexes>();
  rule namespaced = r.add_namespace("ns_test");

  CHECK(namespaced.sources().inSet(variable("ns_test@src_idx")));
  CHECK(namespaced.targets().inSet(variable("ns_test@tgt_idx")));
  CHECK_FALSE(namespaced.sources().inSet(variable("src_idx")));
  CHECK_FALSE(namespaced.targets().inSet(variable("tgt_idx")));
}

// Verifies rule_impl::set_variable_times mutates descriptor source/target
// variables directly, matching rule promotion semantics at impl level.
TEST_CASE("set_variable_times updates rule_impl descriptor variables") {
  rule_implP impl = new copy_rule_impl<simple_rule_for_indexes>;
  time_ident n_level("n_set", time_ident());
  impl->set_variable_times(n_level);
  const rule_impl::info &info = impl->get_info();

  CHECK(desc_has_var(info.sources, variable(variable("src_idx"), n_level)));
  CHECK(desc_has_var(info.targets, variable(variable("tgt_idx"), n_level)));
}

// Verifies split_constraints() moves selected constraints into
// dynamic_constraints, leaves non-selected constraints in constraints, and
// updates rule_identifier with the dynamic-constraint marker.
TEST_CASE("split_constraints moves selected constraints to dynamic_constraints") {
  rule_implP impl = new copy_rule_impl<constraint_split_rule>;
  variableSet dynamic_set(expression::create("c_split1"));
  impl->split_constraints(dynamic_set);
  const rule_impl::info &info = impl->get_info();

  CHECK(desc_has_var(info.dynamic_constraints, variable("c_split1")));
  CHECK(desc_has_var(info.constraints, variable("c_split2")));
  CHECK_FALSE(desc_has_var(info.constraints, variable("c_split1")));
  CHECK(info.rule_identifier().find("DYNAMIC_CONSTRAINT(") != std::string::npos);
}

// Verifies input_vars()/output_vars() accounting for a rule with mapping,
// constraints, and conditionals: inputs include all dependencies, outputs
// include only computed targets.
TEST_CASE("input_vars and output_vars contain expected variables for mapping and conditionals") {
  rule_implP impl = new copy_rule_impl<mapping_conditional_rule>;
  const rule_impl::info &info = impl->get_info();
  variableSet in = info.input_vars();
  variableSet out = info.output_vars();

  CHECK(in.inSet(variable("src_mc")));
  CHECK(in.inSet(variable("c_mc")));
  CHECK(in.inSet(variable("cond_mc")));
  CHECK(in.inSet(variable("m_in")));
  CHECK(in.inSet(variable("m_out")));
  CHECK(desc_has_map_var(info.sources, variable("m_in")));
  CHECK(desc_has_map_var(info.targets, variable("m_out")));
  CHECK(out.inSet(variable("tgt_mc")));
  CHECK_FALSE(out.inSet(variable("src_mc")));
}

// Verifies rule_identifier() includes core dependency shape text useful for
// diagnostics: source name, target name, and "<-" relation marker.
TEST_CASE("rule_identifier includes source and target variable names") {
  rule_implP impl = new copy_rule_impl<simple_rule_for_indexes>;
  const std::string rid = impl->get_info().rule_identifier();

  CHECK(rid.find("src_idx") != std::string::npos);
  CHECK(rid.find("tgt_idx") != std::string::npos);
  CHECK(rid.find("<-") != std::string::npos);
}

// Verifies get_rule(name) resolves to a compatible stored rule implementation
// and preserves expected target descriptor membership.
TEST_CASE("rule_db get_rule returns a clone with expected target var") {
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();
  std::string nm("simple_rule_for_indexes");

  rdb.add_rule(r);
  rule_implP got = rdb.get_rule(nm);
  CHECK(impl_has_target(got, variable("tgt_idx")));
}

// Verifies deep priority targets are indexed at each progressive priority
// level, not only at the fully-qualified priority path.
TEST_CASE("add_rule indexes deep priority chain targets at each priority level") {
  rule_db rdb;
  rule r = make_rule<deep_priority_target_rule>();

  rdb.add_rule(r);
  CHECK(rdb.rules_by_target(variable("p1::p2::tgt_pri")).inSet(r));
  CHECK(rdb.rules_by_target(variable("p2::tgt_pri")).inSet(r));
  CHECK(rdb.rules_by_target(variable("tgt_pri")).inSet(r));
}

// Verifies batch remove_rules() clears all source/target reverse indexes for
// each rule in the input ruleSet.
TEST_CASE("remove_rules clears indexes for every rule in the input ruleSet") {
  rule_db rdb;
  rule r1 = make_rule<simple_rule_for_indexes>();
  rule r2 = make_rule<simple_rule_for_indexes_b>();
  ruleSet rs;

  rdb.add_rule(r1);
  rdb.add_rule(r2);
  rs += r1;
  rs += r2;

  rdb.remove_rules(rs);

  CHECK_FALSE(rdb.rules_by_target(variable("tgt_idx")).inSet(r1));
  CHECK_FALSE(rdb.rules_by_source(variable("src_idx")).inSet(r1));
  CHECK_FALSE(rdb.rules_by_target(variable("tgt_idx_b")).inSet(r2));
  CHECK_FALSE(rdb.rules_by_source(variable("src_idx_b")).inSet(r2));
}

// Verifies default/optional rules are routed into dedicated rule_db categories
// and do not appear in the normal all_rules() collection.
TEST_CASE("default and optional rules are routed to their dedicated rule sets") {
  rule_db rdb;
  rule dr = make_rule<default_category_rule>();
  rule orule = make_rule<optional_category_rule>();

  rdb.add_rule(dr);
  rdb.add_rule(orule);

  CHECK(rdb.get_default_rules().inSet(dr));
  CHECK(rdb.get_optional_rules().inSet(orule));
  CHECK_FALSE(rdb.all_rules().inSet(dr));
  CHECK_FALSE(rdb.all_rules().inSet(orule));
}

// Verifies external rule::sources() includes mapping variables, explicit data
// inputs, constraints, and conditionals because all are source-side deps.
TEST_CASE("mapping and conditional variables appear in sources() for external rules") {
  rule r = make_rule<mapping_conditional_rule>();
  const variableSet &srcs = r.sources();

  CHECK(srcs.inSet(variable("m_in")));
  CHECK(srcs.inSet(variable("m_out")));
  CHECK(srcs.inSet(variable("src_mc")));
  CHECK(srcs.inSet(variable("c_mc")));
  CHECK(srcs.inSet(variable("cond_mc")));
}

// Verifies rename_vars() on external rules rewrites names consistently across
// mappings, targets, constraints, and conditionals.
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

  CHECK(ri.sources().inSet(variable("m_in_new")));
  CHECK(ri.sources().inSet(variable("m_out_new")));
  CHECK(ri.sources().inSet(variable("src_mc_new")));
  CHECK(ri.sources().inSet(variable("c_mc_new")));
  CHECK(ri.sources().inSet(variable("cond_mc_new")));
  CHECK(ri.targets().inSet(variable("tgt_mc_new")));
  CHECK_FALSE(ri.targets().inSet(variable("tgt_mc")));
}

// Verifies set_variable_times() also time-shifts conditional variables so
// conditional evaluation stays synchronized with promoted dependencies.
TEST_CASE("set_variable_times also updates conditional variables") {
  rule_implP impl = new copy_rule_impl<mapping_conditional_rule>;
  time_ident n_level("n_cond", time_ident());
  impl->set_variable_times(n_level);
  const rule_impl::info &info = impl->get_info();

  CHECK(info.conditionals.inSet(variable(variable("cond_mc"), n_level)));
  CHECK_FALSE(info.conditionals.inSet(variable("cond_mc")));
}

// Verifies rule::rename() updates the printable rule name used by stream output
// and diagnostics.
TEST_CASE("rule rename overrides printable name") {
  rule r = make_rule<simple_rule_for_indexes>();
  r.rename("simple_rule_alias");
  std::ostringstream ss;
  ss << r;
  CHECK(ss.str() == "simple_rule_alias");
}

// Verifies duplicate add_rule() is de-duplicated in known_rules and emits a
// warning so accidental duplicate insertion is visible.
TEST_CASE("duplicate add_rule emits warning and does not grow known_rules") {
  CErrCapture capture;
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();

  rdb.add_rule(r);
  const size_t n_before = count_rules(rdb.all_rules());
  rdb.add_rule(r);
  const size_t n_after = count_rules(rdb.all_rules());

  CHECK(n_before == n_after);
  CHECK(capture.str().find("Warning, adding duplicate rule to rule database") != std::string::npos);
}

// Verifies remove_rule() performs enough cleanup that the same rule can be
// re-added and re-indexed correctly.
TEST_CASE("remove_rule followed by add_rule allows reinsert of same rule id") {
  rule_db rdb;
  rule r = make_rule<simple_rule_for_indexes>();

  rdb.add_rule(r);
  rdb.remove_rule(r);
  rdb.add_rule(r);

  CHECK(rdb.rules_by_source(variable("src_idx")).inSet(r));
  CHECK(rdb.rules_by_target(variable("tgt_idx")).inSet(r));
}

// Verifies replace_map_constraints() rewrites map-typed constraints to
// generated constraint facts, keeps non-map constraints intact, and creates a
// constraint-typed fact in fact_db.
TEST_CASE("replace_map_constraints converts map constraints to generated constraints") {
  rule_implP impl = new copy_rule_impl<map_constraint_replace_rule>;
  fact_db facts;

  Map map_con;
  map_con.allocate(interval(1, 3));
  facts.create_fact("map_con", map_con);

  Constraint keep_con;
  *keep_con = interval(2, 4);
  facts.create_fact("keep_con", keep_con);

  impl->replace_map_constraints(facts);
  const rule_impl::info &info = impl->get_info();

  CHECK(desc_has_var(info.constraints, variable("__map_con_MAP_constraint__")));
  CHECK(desc_has_var(info.constraints, variable("keep_con")));
  CHECK_FALSE(desc_has_var(info.constraints, variable("map_con")));

  storeRepP generated = facts.get_variable("__map_con_MAP_constraint__");
  REQUIRE(generated != static_cast<storeRep *>(0));
  CHECK(isCONSTRAINT(generated));
}

// Verifies replace_map_constraints() is idempotent: repeated calls should not
// destabilize generated map-constraint facts or their domain.
TEST_CASE("replace_map_constraints is idempotent for generated constraint facts") {
  rule_implP impl = new copy_rule_impl<map_constraint_replace_rule>;
  fact_db facts;

  Map map_con;
  map_con.allocate(interval(7, 9));
  facts.create_fact("map_con", map_con);

  Constraint keep_con;
  *keep_con = interval(3, 5);
  facts.create_fact("keep_con", keep_con);

  impl->replace_map_constraints(facts);
  storeRepP first = facts.get_variable("__map_con_MAP_constraint__");
  REQUIRE(first != static_cast<storeRep *>(0));
  const entitySet first_dom = first->domain();

  impl->replace_map_constraints(facts);
  storeRepP second = facts.get_variable("__map_con_MAP_constraint__");
  REQUIRE(second != static_cast<storeRep *>(0));
  CHECK(second->domain() == first_dom);
}

// Verifies check_perm_bits() with a small positive/negative matrix across rule
// categories: valid target container types pass, intentionally wrong ones fail.
TEST_CASE("check_perm_bits type matrix catches bad and accepts good rule targets") {
  rule_implP good_default = new copy_rule_impl<default_category_rule>;
  rule_implP bad_default = new copy_rule_impl<default_bad_store_target_rule>;
  rule_implP bad_pointwise = new copy_rule_impl<pointwise_bad_param_target_rule>;
  rule_implP good_map = new copy_rule_impl<map_good_rule>;
  rule_implP bad_map = new copy_rule_impl<map_bad_store_target_rule>;
  rule_implP good_constraint = new copy_rule_impl<constraint_good_rule>;
  rule_implP bad_constraint = new copy_rule_impl<constraint_bad_store_target_rule>;

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
