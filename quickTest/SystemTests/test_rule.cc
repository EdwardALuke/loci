#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

namespace {

  /// Minimal pointwise rule with a priority annotation on its target.
  class priority_target_rule : public pointwise_rule {
    const_store<int> source;
    store<int> target;

  public:
    priority_target_rule() {
      name_store("rule_test_source", source);
      name_store("priority::rule_test_target", target);
      input("rule_test_source");
      output("priority::rule_test_target");
    }

    void compute(const sequence &) {}
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0); }
  };

  /// Creates the priority-target rule used by the rule_db test.
  rule make_priority_target_rule() {
    return rule(rule_implP(new copy_rule_impl<priority_target_rule>));
  }

} // namespace

/// A rule for priority::x is found as both priority::x and general x; removal
/// must clear both lookups.
TEST_CASE("rule removal clears priority-annotated and general target lookups") {
  rule_db rule_database;
  rule candidate = make_priority_target_rule();
  const variable priority_target("priority::rule_test_target");
  const variable general_target("rule_test_target");

  rule_database.add_rule(candidate);
  CHECK(rule_database.rules_by_target(priority_target).inSet(candidate));
  CHECK(rule_database.rules_by_target(general_target).inSet(candidate));

  rule_database.remove_rule(candidate);
  CHECK_FALSE(rule_database.rules_by_target(priority_target).inSet(candidate));
  CHECK_FALSE(rule_database.rules_by_target(general_target).inSet(candidate));
}

/// The source set must include the input value, both maps, the constraint, and
/// the conditional needed by the rule.
TEST_CASE("rules report every required input through sources") {
  rule rule_with_inputs(
    "qualifier(dependencies),"
    "source(source_map->input),"
    "target(target_map->output),"
    "constraint(rule_constraint),"
    "conditional(rule_condition)");

  CHECK(rule_with_inputs.sources().inSet(variable("input")));
  CHECK(rule_with_inputs.sources().inSet(variable("source_map")));
  CHECK(rule_with_inputs.sources().inSet(variable("target_map")));
  CHECK(rule_with_inputs.sources().inSet(variable("rule_constraint")));
  CHECK(rule_with_inputs.sources().inSet(variable("rule_condition")));
  CHECK(rule_with_inputs.targets().inSet(variable("output")));
}

/// Prepending a time level must move the map, source, constraint, conditional,
/// and target together.
TEST_CASE("prepend_rule preserves rule structure at the new time level") {
  rule original(
    "qualifier(rule_test),"
    "source(source_map->source_value),"
    "target(target_value),"
    "constraint(rule_constraint),"
    "conditional(rule_condition)");
  time_ident time_level("n", time_ident());

  rule prepended = prepend_rule(original, time_level);

  CHECK(prepended.get_info().maps().inSet(
    variable(time_level, variable("source_map"))));
  CHECK(prepended.sources().inSet(
    variable(time_level, variable("source_value"))));
  CHECK(prepended.get_info().constraints().inSet(
    variable(time_level, variable("rule_constraint"))));
  CHECK(prepended.sources().inSet(
    variable(time_level, variable("rule_condition"))));
  CHECK(prepended.targets().inSet(
    variable(time_level, variable("target_value"))));
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
