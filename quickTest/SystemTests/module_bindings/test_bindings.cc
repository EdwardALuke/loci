// Copyright 2026, Streamline Numerics, Inc.
// Distributed under LGPL-3.0-or-later; see COPYING.LESSER.
#include <Loci.h>
#include <module_bindings.h>
#include <param_rule.h>
#include <stdexcept>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci ;

namespace {
  void put(fact_db &facts, const char *name, int value) {
    param<int> data ;
    *data = value ;
    facts.create_fact(variable(name), data.Rep()) ;
  }

  int get(fact_db &facts, const char *name) {
    const_param<int> data(facts.get_variable(variable(name))) ;
    return *data ;
  }

  void load(const char *library, const char *instance,
            fact_db &facts, rule_db &rules) {
    std::set<std::string> loaded ;
    // Use the application loader; these fixtures have no setup callback.
    load_module(library, instance, "unused", facts, rules, loaded) ;
  }
}

TEST_CASE("binding identities preserve annotations and argument labels") {
  const module_bindings bindings({
    {variable("plate@pressure"), variable("exchange@pressureOnSolid")}
  }) ;
  CHECK(bindings.resolve(variable("plate@pressure"))
        == variable("exchange@pressureOnSolid")) ;
  CHECK(bindings.resolve(variable("plate@pressure{n,it+1}"))
        == variable("exchange@pressureOnSolid{n,it+1}")) ;
  CHECK(bindings.resolve(variable("plate@pressure{n=0}"))
        == variable("exchange@pressureOnSolid{n=0}")) ;
  CHECK(bindings.resolve(variable("fast::plate@pressure{n}"))
        == variable("fast::exchange@pressureOnSolid{n}")) ;
  CHECK(bindings.resolve(variable("plate@Copy(plate@pressure)"))
        == variable("plate@Copy(plate@pressure)")) ;
  CHECK(bindings.resolve(variable("plate@pressure_ic"))
        == variable("plate@pressure_ic")) ;
  const variable once = bindings.resolve(variable("plate@pressure{n}")) ;
  CHECK(bindings.resolve(once) == once) ;

  // A root destination is absolute, not relative to the importing instance.
  const module_bindings root({{variable("plate@pressure"), variable("p")}}) ;
  CHECK(root.resolve(variable("plate@pressure")) == variable("p")) ;
}

TEST_CASE("invalid tables and rebinding fail before scheduling") {
  CHECK_THROWS_AS(module_bindings({
    {variable("plate@pressure"), variable()}
  }), std::invalid_argument) ;
  CHECK_THROWS_AS(module_bindings({
    {variable("pressure"), variable("host@pressure")}
  }), std::invalid_argument) ;
  CHECK_THROWS_AS(module_bindings({
    {variable("plate@pressure"), variable("host@Copy(pressure)")}
  }), std::invalid_argument) ;
  CHECK_THROWS_AS(module_bindings({
    {variable("plate@pressure{n}"), variable("host@pressure")}
  }), std::invalid_argument) ;
  CHECK_THROWS_AS(module_bindings({
    {variable("plate@pressure"), variable("other@pressure")},
    {variable("other@pressure"), variable("host@pressure")}
  }), std::invalid_argument) ;

  rule_db rules ;
  const module_bindings first({
    {variable("plate@pressure"), variable("host@pressure")}
  }) ;
  rules.set_module_bindings(first) ;
  CHECK_NOTHROW(rules.set_module_bindings(first)) ;
  CHECK_THROWS_AS(rules.set_module_bindings(module_bindings({
    {variable("plate@pressure"), variable("other@pressure")}
  })), std::invalid_argument) ;
  rule_db copied(rules) ;
  CHECK(copied.get_module_bindings().entries() == first.entries()) ;
}

TEST_CASE("two loaded instances consume computed partner facts") {
  bool reverse = false ;
  bool shared = false ;
  bool legacy = false ;
  SUBCASE("legacy private facts without bindings") { legacy = true ; }
  SUBCASE("intentionally shared partner") { shared = true ; }
  SUBCASE("reverse declaration order") { reverse = true ; }
  SUBCASE("distinct partners") {}

  fact_db facts ;
  rule_db rules ;
  load("binding_consumer", reverse ? "right" : "left", facts, rules) ;
  load("binding_consumer", reverse ? "left" : "right", facts, rules) ;
  load("binding_provider", "a", facts, rules) ;
  load("binding_provider", "b", facts, rules) ;
  if(!legacy)
    rules.set_module_bindings(module_bindings({
      {variable("left@pressure"), variable("a@pressure")},
      {variable("right@pressure"),
       variable(shared ? "a@pressure" : "b@pressure")}
    })) ;
  put(facts, "a@seed", 2) ;
  put(facts, "b@seed", 5) ;
  put(facts, "left@normal", 3) ;
  put(facts, "right@normal", 7) ;
  if(legacy) {
    put(facts, "left@pressure", 2) ;
    put(facts, "right@pressure", 5) ;
  }

  // Inspect the expansion itself, not only the final arithmetic.
  const variable family("left@Copy(left@pressure)") ;
  variableSet request ;
  request += family ;
  rule_db expanded = parametric_rdb(rules, request) ;
  rule_db resolved = resolve_module_bindings(
    expanded, rules.get_module_bindings()) ;
  bool found = false ;
  for(rule provider : resolved.all_rules()) {
    rule_implP impl = provider.get_rule_implP() ;
    if(impl->get_info().output_vars().inSet(family)) {
      found = true ;
      CHECK(impl->get_info().input_vars().inSet(variable(
        legacy ? "left@pressure" : "a@pressure"))) ;
    }
  }
  REQUIRE(found) ;
  // The source templates remain available for another composition/query.
  bool original = false ;
  for(rule provider : expanded.all_rules())
    if(provider.get_rule_implP()->get_info().output_vars().inSet(family))
      original = provider.get_rule_implP()->get_info().input_vars()
        .inSet(variable("left@pressure")) ;
  CHECK(original) ;

  const char *query = "left@traction,right@traction,left@copied,right@copied" ;
  REQUIRE(makeQuery(rules, facts, query)) ;
  CHECK(get(facts, "left@traction") == 6) ;
  CHECK(get(facts, "right@traction") == (shared ? 14 : 35)) ;
  CHECK(get(facts, "left@copied") == 2) ;
  CHECK(get(facts, "right@copied") == (shared ? 2 : 5)) ;

  if(!legacy) {
    // Querying an imported spelling publishes only the canonical fact.
    REQUIRE(makeQuery(rules, facts, "left@pressure")) ;
    CHECK(get(facts, "a@pressure") == 2) ;
    CHECK_FALSE(facts.get_typed_variables().inSet(variable("left@pressure"))) ;

    // Copied databases and changed inputs must not reuse stale results.
    fact_db next_facts(facts) ;
    rule_db next_rules(rules) ;
    put(next_facts, "a@seed", 4) ;
    REQUIRE(makeQuery(next_rules, next_facts, query)) ;
    CHECK(get(next_facts, "left@traction") == 12) ;
    CHECK(get(next_facts, "left@copied") == 4) ;
    CHECK(get(next_facts, "right@traction") == (shared ? 28 : 35)) ;
  }
}

TEST_CASE("a target exposed by substitution cannot write an imported input") {
  fact_db facts ;
  rule_db rules ;
  load("binding_writer", "bad", facts, rules) ;
  rules.set_module_bindings(module_bindings({
    {variable("bad@pressure"), variable("host@pressure")}
  })) ;
  put(facts, "bad@seed", 9) ;
  put(facts, "host@pressure", 17) ;
  CHECK_NOTHROW(validate_module_bindings(rules, rules.get_module_bindings())) ;
  CHECK_THROWS_WITH_AS(makeQuery(rules, facts, "bad@result"),
    doctest::Contains("rule writes imported input bad@pressure"),
    std::invalid_argument) ;
  CHECK(get(facts, "host@pressure") == 17) ;
}

TEST_CASE("defaults, reductions and in-place writes cannot supply imports") {
  fact_db facts ;
  rule_db rules ;
  load("binding_writer", "bad", facts, rules) ;
  variable input("bad@prescribed") ;
  SUBCASE("default") {}
  SUBCASE("unit and apply") { input = variable("bad@total") ; }
  SUBCASE("in-place original") { input = variable("bad@original") ; }
  SUBCASE("in-place result") { input = variable("bad@updated") ; }
  const module_bindings bindings({{input, variable("host@value")}}) ;
  // Check each matching writer separately so a unit cannot hide a bad apply.
  ruleSet providers = rules.all_rules() ;
  providers += rules.get_default_rules() ;
  int checked = 0 ;
  for(rule provider : providers) {
    rule_implP impl = provider.get_rule_implP() ;
    if(!impl->get_var_list().inSet(input))
      continue ;
    // The apply merely reads original; only the in-place rule writes it.
    if(input == variable("bad@original")
       && impl->get_rule_class() != rule_impl::POINTWISE)
      continue ;
    rule_db single ;
    single.add_rule(provider) ;
    single.set_module_bindings(bindings) ;
    put(facts, "host@value", 17) ;
    CHECK_THROWS_WITH_AS(makeQuery(single, facts, "host@value"),
      doctest::Contains("rule writes imported input"), module_binding_error) ;
    CHECK(get(facts, "host@value") == 17) ;
    ++checked ;
  }
  CHECK(checked == (input == variable("bad@total") ? 2 : 1)) ;
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv) ;
  doctest::Context context ;
  context.applyCommandLine(argc, argv) ;
  const int result = context.run() ;
  Loci::Finalize() ;
  return result ;
}
