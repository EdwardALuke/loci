#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <list>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Loci;

namespace {

/// @brief Returns an inclusive `entitySet` interval `[first, last]`.
/// @param[in] first First entity in the interval.
/// @param[in] last Last entity in the interval.
/// @return An `entitySet` spanning the requested inclusive interval.
entitySet make_range(int first, int last) {
  return entitySet(interval(first, last));
}

/// @brief Allocates `values` on `domain` and fills it from `entries`.
/// @param[in,out] values Store to allocate and populate.
/// @param[in] domain Domain to allocate in the store.
/// @param[in] entries Values copied into the store in domain iteration order.
/// @throws std::logic_error When `entries.size()` does not match `domain.size()`.
void assign_store_values(store<int> &values, const entitySet &domain,
                         const std::vector<int> &entries) {
  if(static_cast<size_t>(domain.size()) != entries.size()) {
    throw std::logic_error("store seed data does not match the domain size");
  }

  values.allocate(domain);
  size_t index = 0;
  for(entitySet::const_iterator ei = domain.begin(); ei != domain.end();
      ++ei, ++index) {
    values[*ei] = entries[index];
  }
}

/// @brief Allocates `mapping` on `domain` and fills it from `entries`.
/// @param[in,out] mapping Map to allocate and populate.
/// @param[in] domain Domain to allocate in the map.
/// @param[in] entries Target entities copied into the map in domain iteration
///   order.
/// @throws std::logic_error When `entries.size()` does not match `domain.size()`.
void assign_map_values(Map &mapping, const entitySet &domain,
                       const std::vector<int> &entries) {
  if(static_cast<size_t>(domain.size()) != entries.size()) {
    throw std::logic_error("map seed data does not match the domain size");
  }

  mapping.allocate(domain);
  size_t index = 0;
  for(entitySet::const_iterator ei = domain.begin(); ei != domain.end();
      ++ei, ++index) {
    mapping[*ei] = entries[index];
  }
}

/// @brief Returns true when `rep` wraps a non-null representation pointer.
/// @param[in] rep Store representation handle to inspect.
/// @return True when `rep` is non-null.
bool has_rep(const storeRepP &rep) {
  return rep != static_cast<storeRep *>(0);
}

/// @brief Returns true when both handles refer to the same store representation.
/// @param[in] lhs First store representation handle.
/// @param[in] rhs Second store representation handle.
/// @return True when both handles are non-null and share the same underlying
///   representation pointer.
bool same_rep(const storeRepP &lhs, const storeRepP &rhs) {
  if(!has_rep(lhs) || !has_rep(rhs)) {
    return false;
  }
  return lhs.operator->() == rhs.operator->();
}

/// @brief Builds a rule descriptor string that `rule(...)` can parse.
/// @param[in] source_expr Source portion of the rule signature.
/// @param[in] target_expr Target portion of the rule signature.
/// @param[in] qualifier Qualifier portion of the rule signature.
/// @return A `rule` constructed from the assembled signature.
rule make_rule(const std::string &source_expr, const std::string &target_expr,
               const std::string &qualifier) {
  std::ostringstream signature;
  signature << "source(" << source_expr << "),target(" << target_expr
            << "),qualifier(" << qualifier << ")";
  return rule(signature.str());
}

/// @brief Collects `vars` into a `variableSet`.
/// @param[in] vars Variables to insert.
/// @return A `variableSet` containing every variable from `vars`.
variableSet make_var_set(const std::vector<variable> &vars) {
  variableSet set;
  for(size_t i = 0; i < vars.size(); ++i) {
    set += vars[i];
  }
  return set;
}

/// @brief Copies a `std::list<comm_info>` into an indexable vector.
/// @param[in] entries Communication records to copy.
/// @return A vector containing the same records in the same order.
std::vector<comm_info> to_vector(const std::list<comm_info> &entries) {
  return std::vector<comm_info>(entries.begin(), entries.end());
}

} // namespace

TEST_CASE("sched_db constructor mirrors fact_db variables and existence") {
  fact_db facts;
  store<int> cells;

  // Seed one typed fact so the constructor has scheduler state to mirror.
  assign_store_values(cells, make_range(4, 6), {10, 20, 30});
  facts.create_fact("cells", cells);

  sched_db sched(facts);

  // Verify the constructor registers the variable and its self-alias/synonym
  // bookkeeping.
  CHECK(sched.get_typed_variables().inSet(variable("cells")));
  CHECK(sched.get_aliases(variable("cells")).inSet(variable("cells")));
  CHECK(sched.get_synonyms(variable("cells")).inSet(variable("cells")));

  // Verify the mirrored fact domain becomes known existence and leaves the
  // remaining scheduler state empty/defaulted.
  CHECK(sched.variable_existence(variable("cells")) == make_range(4, 6));
  CHECK_FALSE(sched.is_a_Map(variable("cells")));
  CHECK(sched.get_variable_requests(variable("cells")) == EMPTY);
}

TEST_CASE("set_variable_type installs intensional facts and empty schedule state") {
  fact_db facts;
  sched_db sched(facts);
  store<int> derived;

  // Build a store representation that will be installed as an intensional fact.
  assign_store_values(derived, make_range(10, 11), {3, 5});

  sched.set_variable_type(variable("derived"), derived.Rep(), facts);

  // Verify the new variable is tracked by both sched_db and fact_db, and that
  // sched_db starts it with no derived existence yet.
  CHECK(sched.get_typed_variables().inSet(variable("derived")));
  CHECK(facts.get_intensional_facts().inSet(variable("derived")));
  CHECK_FALSE(facts.get_extensional_facts().inSet(variable("derived")));
  CHECK(sched.variable_existence(variable("derived")) == EMPTY);

  // Verify the installed fact points at the expected store contents.
  storeRepP rep = facts.get_variable("derived");
  REQUIRE(has_rep(rep));
  store<int> fetched(rep);
  CHECK(fetched.domain() == make_range(10, 11));
  CHECK(fetched[10] == 3);
  CHECK(fetched[11] == 5);
}

TEST_CASE("set_variable_type preserves map metadata for intensional variables") {
  fact_db facts;
  sched_db sched(facts);
  Map derived_map;
  assign_map_values(derived_map, make_range(1, 3), {10, 11, 11});

  sched.set_variable_type(variable("derived_map"), derived_map.Rep(), facts);

  CHECK(sched.get_typed_variables().inSet(variable("derived_map")));
  CHECK(facts.get_intensional_facts().inSet(variable("derived_map")));
  CHECK(sched.variable_existence(variable("derived_map")) == EMPTY);
  CHECK(sched.is_a_Map(variable("derived_map")));
  CHECK(sched.image(variable("derived_map"), make_range(1, 3)) ==
        make_range(10, 11));

  std::pair<entitySet, entitySet> preimage =
      sched.preimage(variable("derived_map"), make_range(11, 11));
  CHECK(preimage.first == make_range(2, 3));
  CHECK(preimage.second == make_range(2, 3));
}

TEST_CASE("alias_variable keeps alias bookkeeping in sync with fact_db") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(1, 3), {7, 8, 9});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.alias_variable(variable("cells"), variable("cells_alias"), facts);

  variableSet aliases = sched.get_aliases(variable("cells"));

  // Verify sched_db records both names in the shared alias metadata.
  CHECK(aliases.inSet(variable("cells")));
  CHECK(aliases.inSet(variable("cells_alias")));
  CHECK(sched.get_aliases(variable("cells_alias")) == aliases);
  CHECK(sched.get_antialiases(variable("cells_alias"))
            .inSet(variable("cells")));
  CHECK(sched.get_typed_variables().inSet(variable("cells_alias")));

  // Verify fact_db mirrors the alias relationship onto the same store
  // representation.
  CHECK(same_rep(facts.get_variable("cells"), facts.get_variable("cells_alias")));
}

TEST_CASE("alias_variable accepts either argument order when one name already exists") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(7, 9), {3, 4, 5});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.alias_variable(variable("cells_alias"), variable("cells"));

  variableSet aliases = sched.get_aliases(variable("cells"));

  // Verify the helper resolves the existing variable as the alias source rather
  // than treating the reversed argument order as an error.
  CHECK_FALSE(sched.errors_found());
  CHECK(aliases.inSet(variable("cells")));
  CHECK(aliases.inSet(variable("cells_alias")));
  CHECK(sched.get_aliases(variable("cells_alias")) == aliases);
  CHECK(sched.get_antialiases(variable("cells_alias"))
            .inSet(variable("cells")));
  CHECK(sched.get_typed_variables().inSet(variable("cells_alias")));
}

TEST_CASE("synonym_variable records canonical lookup behavior") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(2, 4), {1, 2, 3});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.synonym_variable(variable("cells"), variable("cells_shadow"), facts);

  // Verify the synonym resolves back to the canonical scheduler name.
  CHECK(sched.remove_synonym(variable("cells_shadow")) == variable("cells"));

  variableSet synonyms = sched.get_synonyms(variable("cells"));

  // Verify the synonym set and fact_db both treat the two names as one logical
  // variable.
  CHECK(synonyms.inSet(variable("cells")));
  CHECK(synonyms.inSet(variable("cells_shadow")));
  CHECK(sched.get_typed_variables().inSet(variable("cells_shadow")));
  CHECK(
      same_rep(facts.get_variable("cells"), facts.get_variable("cells_shadow")));
}

TEST_CASE("sched-only synonym_variable reports duplicate existing names" *
          doctest::may_fail()) {
  fact_db facts;
  store<int> cells;
  store<int> faces;
  assign_store_values(cells, make_range(1, 2), {1, 2});
  assign_store_values(faces, make_range(3, 4), {3, 4});
  facts.create_fact("cells", cells);
  facts.create_fact("faces", faces);

  sched_db sched(facts);
  sched.clear_errors();

  // The fact_db overload catches this; the sched-only overload currently logs
  // the error without setting the sticky scheduler error flag.
  sched.synonym_variable(variable("cells"), variable("faces"));

  CHECK(sched.errors_found());
}

TEST_CASE("synonyms resolve requests and direct existence through the canonical variable") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(2, 4), {1, 2, 3});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.synonym_variable(variable("cells"), variable("cells_shadow"), facts);
  sched.set_variable_existence(variable("cells_shadow"), make_range(20, 21));
  sched.variable_request(variable("cells_shadow"), make_range(21, 22));

  entitySet expected_existence = make_range(2, 4);
  expected_existence += make_range(20, 21);

  // Verify both names resolve through the same canonical scheduler entry.
  CHECK(sched.remove_synonym(variable("cells_shadow")) == variable("cells"));
  CHECK(sched.variable_existence(variable("cells")) == expected_existence);
  CHECK(sched.variable_existence(variable("cells_shadow")) == expected_existence);
  CHECK(sched.get_variable_requests(variable("cells")) == make_range(21, 22));
  CHECK(sched.get_variable_requests(variable("cells_shadow")) ==
        make_range(21, 22));
}

TEST_CASE("rotation groups are recorded and try_get helpers stay safe for unknown variables") {
  fact_db facts;
  store<int> cells;
  store<int> faces;
  assign_store_values(cells, make_range(1, 2), {10, 20});
  assign_store_values(faces, make_range(3, 4), {30, 40});
  facts.create_fact("cells", cells);
  facts.create_fact("faces", faces);

  sched_db sched(facts);
  variableSet rotations = make_var_set({variable("cells"), variable("faces")});
  sched.set_variable_rotations(rotations);

  // Verify every member of the rotation group sees the full recorded set.
  CHECK(sched.get_rotations(variable("cells")) == rotations);
  CHECK(sched.get_rotations(variable("faces")) == rotations);

  // Verify the unknown-safe accessors return empty results instead of aborting.
  CHECK(sched.try_get_synonyms(variable("ghost")) == variableSet(EMPTY));
  CHECK(sched.try_get_aliases(variable("ghost")) == variableSet(EMPTY));
  CHECK(sched.try_get_antialiases(variable("ghost")) == variableSet(EMPTY));
  CHECK(sched.try_get_rotations(variable("ghost")) == variableSet(EMPTY));
}

TEST_CASE("aliases share scheduling data but keep separate existence and request state") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(1, 2), {10, 20});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.alias_variable(variable("cells"), variable("cells_alias"), facts);
  sched.set_variable_existence(variable("cells"), make_range(30, 31));
  sched.variable_request(variable("cells"), make_range(31, 32));

  // Verify the alias still participates in the shared alias group.
  CHECK(sched.get_aliases(variable("cells_alias")).inSet(variable("cells")));
  CHECK(sched.get_aliases(variable("cells_alias"))
            .inSet(variable("cells_alias")));

  // Verify existence and requests remain separate because aliases keep distinct
  // sched_info records.
  CHECK(sched.variable_existence(variable("cells_alias")) == EMPTY);
  CHECK(sched.get_variable_requests(variable("cells_alias")) == EMPTY);
}

TEST_CASE("aliases of map variables share map metadata and cached map queries") {
  fact_db facts;
  Map cell_to_parent;
  assign_map_values(cell_to_parent, make_range(1, 4), {10, 10, 11, 13});
  facts.create_fact("cell_to_parent", cell_to_parent);

  sched_db sched(facts);
  sched.alias_variable(variable("cell_to_parent"),
                       variable("cell_to_parent_alias"), facts);

  // Alias variables share sched_data, so map-ness and cached map queries should
  // behave exactly like the canonical map variable.
  CHECK(sched.is_a_Map(variable("cell_to_parent_alias")));
  CHECK(sched.image(variable("cell_to_parent_alias"), make_range(1, 3)) ==
        make_range(10, 11));

  std::pair<entitySet, entitySet> preimage =
      sched.preimage(variable("cell_to_parent_alias"), make_range(10, 10));
  CHECK(preimage.first == make_range(1, 2));
  CHECK(preimage.second == make_range(1, 2));
}

TEST_CASE("requests are intersected with rule existence and can be cleared") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(10, 12), {4, 5, 6});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  rule producer = make_rule("upstream", "cells", "rule_one");

  sched.set_existential_info(variable("cells"), producer, make_range(10, 12));
  sched.variable_request(variable("cells"), make_range(11, 13));
  sched.add_extra_unit_request(variable("cells"), make_range(20, 21));

  // Verify the direct getters expose the stored existential and request state.
  CHECK(sched.get_existential_info(variable("cells"), producer) ==
        make_range(10, 12));
  CHECK(sched.get_variable_requests(variable("cells")) == make_range(11, 13));

  // Verify per-rule requests intersect the existential region, and that UNIT
  // extras are tracked separately.
  CHECK(sched.get_variable_request(producer, variable("cells")) ==
        make_range(11, 12));
  CHECK(sched.get_extra_unit_request(variable("cells")) == make_range(20, 21));

  sched.clear_variable_request();

  // Verify clearing requests also clears the extra UNIT-request bookkeeping.
  CHECK(sched.get_variable_requests(variable("cells")) == EMPTY);
  CHECK(sched.get_extra_unit_request(variable("cells")) == EMPTY);
}

TEST_CASE("time variables exist everywhere and ignore request updates") {
  sched_db sched;
  variable time_step("$time_step");

  // Sanity-check the fixture before probing sched_db behavior.
  REQUIRE(time_step.get_info().tvar);

  // Verify time variables behave as globally existing scheduler variables and
  // ignore request accumulation.
  CHECK(sched.variable_existence(time_step) == ~EMPTY);
  CHECK_NOTHROW(sched.variable_request(time_step, make_range(1, 3)));
  CHECK_FALSE(sched.errors_found());
}

TEST_CASE("conflicting existential rules are reported as scheduler errors") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(1, 3), {10, 11, 12});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  rule first = make_rule("left", "cells", "rule_left");
  rule second = make_rule("right", "cells", "rule_right");

  sched.clear_errors();
  sched.set_existential_info(variable("cells"), first, make_range(1, 2));

  // The first producer establishes ownership without conflict.
  CHECK_FALSE(sched.errors_found());

  sched.set_existential_info(variable("cells"), second, make_range(2, 3));

  // The overlapping second producer should trip the sticky error flag while
  // leaving each rule's recorded contribution intact.
  CHECK(sched.errors_found());
  CHECK(sched.get_existential_info(variable("cells"), first) ==
        make_range(1, 2));
  CHECK(sched.get_existential_info(variable("cells"), second) ==
        make_range(2, 3));
}

TEST_CASE("repeated existential writes from the same rule accumulate without conflicts") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(1, 4), {10, 11, 12, 13});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  rule producer = make_rule("left", "cells", "same_rule");
  rule unrelated = make_rule("right", "cells", "unrelated");

  sched.clear_errors();
  sched.set_existential_info(variable("cells"), producer, make_range(1, 2));
  sched.set_existential_info(variable("cells"), producer, make_range(2, 4));
  sched.variable_request(variable("cells"), make_range(3, 5));

  // Verify repeated writes from the same producer merge instead of conflicting.
  CHECK_FALSE(sched.errors_found());
  CHECK(sched.get_existential_info(variable("cells"), producer) ==
        make_range(1, 4));

  // Verify unrelated rules still report no request overlap.
  CHECK(sched.get_variable_request(producer, variable("cells")) ==
        make_range(3, 4));
  CHECK(sched.get_variable_request(unrelated, variable("cells")) == EMPTY);
}

TEST_CASE("map variables expose cached image and preimage queries") {
  fact_db facts;
  Map cell_to_parent;
  store<int> cells;
  assign_map_values(cell_to_parent, make_range(1, 4), {10, 10, 11, 13});
  assign_store_values(cells, make_range(1, 4), {1, 2, 3, 4});
  facts.create_fact("cell_to_parent", cell_to_parent);
  facts.create_fact("cells", cells);

  sched_db sched(facts);

  // Verify map-ness is detected from the shared scheduling data.
  CHECK(sched.is_a_Map(variable("cell_to_parent")));
  CHECK_FALSE(sched.is_a_Map(variable("cells")));

  // Verify forward image queries return the mapped codomain and remain stable
  // across repeated lookups.
  CHECK(sched.image(variable("cell_to_parent"), make_range(1, 3)) ==
        make_range(10, 11));
  CHECK(sched.image(variable("cell_to_parent"), make_range(1, 3)) ==
        make_range(10, 11));

  // Verify preimage queries preserve the underlying map semantics.
  std::pair<entitySet, entitySet> preimage =
      sched.preimage(variable("cell_to_parent"), make_range(10, 10));
  CHECK(preimage.first == make_range(1, 2));
  CHECK(preimage.second == make_range(1, 2));

  // Non-map variables should return the empty fallbacks for both APIs.
  CHECK(sched.image(variable("cells"), make_range(1, 3)) == EMPTY);
  std::pair<entitySet, entitySet> store_preimage =
      sched.preimage(variable("cells"), make_range(10, 10));
  CHECK(store_preimage.first == EMPTY);
  CHECK(store_preimage.second == EMPTY);
}

TEST_CASE("direct existence and shadow helpers update bookkeeping without adding rules") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(1, 2), {5, 6});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  entitySet expected_existence = make_range(1, 2);
  expected_existence += make_range(5, 6);

  // Start from a variable with no explicit existential producers.
  CHECK(sched.get_existential_rules(variable("cells")).size() == 0);

  sched.set_variable_existence(variable("cells"), make_range(5, 6));
  sched.set_variable_shadow(variable("cells"), make_range(20, 20));
  sched.variable_shadow(variable("cells"), make_range(21, 22));

  // Verify direct existence/shadow helpers update bookkeeping without
  // synthesizing existential producer records.
  CHECK(sched.variable_existence(variable("cells")) == expected_existence);
  CHECK(sched.get_existential_rules(variable("cells")).size() == 0);
  CHECK(sched.get_variable_shadow(variable("cells")) == make_range(20, 22));
}

TEST_CASE("print_summary reports the tracked container kind and scheduler state") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(3, 4), {7, 8});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.synonym_variable(variable("cells"), variable("cells_shadow"), facts);
  sched.variable_request(variable("cells"), make_range(4, 4));

  std::ostringstream summary;
  sched.print_summary(facts, summary);
  const std::string text = summary.str();

  // Keep this test resilient by checking for the key pieces of information the
  // summary claims to expose rather than matching exact formatting.
  CHECK(text.find("Summary of Existential deduction:") != std::string::npos);
  CHECK(text.find("Container = STORE") != std::string::npos);
  CHECK(text.find("cells") != std::string::npos);
  CHECK(text.find("cells_shadow") != std::string::npos);
  CHECK(text.find("request=") != std::string::npos);
}

TEST_CASE("duplication policies and flags propagate across synonym sets") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(5, 6), {8, 9});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.synonym_variable(variable("cells"), variable("cells_shadow"), facts);

  sched.add_policy(variable("cells"), sched_db::ALWAYS);
  sched.add_policy(variable("cells"), sched_db::MODEL_BASED);
  sched.set_duplicate_variable(variable("cells"), true);

  // Verify both the policy bitmask and the final duplication choice propagate
  // across the synonym set.
  CHECK(sched.is_policy(variable("cells"), sched_db::ALWAYS));
  CHECK(sched.is_policy(variable("cells_shadow"), sched_db::ALWAYS));
  CHECK(sched.is_policy(variable("cells"), sched_db::MODEL_BASED));
  CHECK(sched.is_policy(variable("cells_shadow"), sched_db::MODEL_BASED));
  CHECK(sched.is_duplicate_variable(variable("cells")));
  CHECK(sched.is_duplicate_variable(variable("cells_shadow")));
}

TEST_CASE("raw policy replacement affects synonyms but not aliases") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(5, 6), {8, 9});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.alias_variable(variable("cells"), variable("cells_alias"), facts);
  sched.synonym_variable(variable("cells"), variable("cells_shadow"), facts);

  sched.add_policy(variable("cells"), sched_db::ALWAYS);
  sched.add_policy(variable("cells"), sched_db::MODEL_BASED);
  sched.set_policy(variable("cells"), 1u);
  sched.set_duplicate_variable(variable("cells"), true);

  // Replacing the raw mask should leave only the NEVER bit on the canonical
  // scheduler record shared by synonyms.
  CHECK(sched.get_policy(variable("cells")) == 1u);
  CHECK(sched.get_policy(variable("cells_shadow")) == 1u);
  CHECK(sched.is_policy(variable("cells"), sched_db::NEVER));
  CHECK(sched.is_policy(variable("cells_shadow"), sched_db::NEVER));
  CHECK_FALSE(sched.is_policy(variable("cells"), sched_db::ALWAYS));
  CHECK_FALSE(sched.is_policy(variable("cells_shadow"), sched_db::MODEL_BASED));

  // Aliases keep their own sched_info, so neither the raw mask nor the final
  // duplicate-variable flag propagates to them automatically.
  CHECK(sched.get_policy(variable("cells_alias")) == 0u);
  CHECK_FALSE(sched.is_duplicate_variable(variable("cells_alias")));
}

TEST_CASE("possible duplicate worklist expands synonym groups without pulling in aliases") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(1, 2), {3, 4});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  sched.alias_variable(variable("cells"), variable("cells_alias"), facts);
  sched.synonym_variable(variable("cells"), variable("cells_shadow"), facts);

  sched.add_possible_duplicate_vars(make_var_set({variable("cells")}));

  variableSet possible = sched.get_possible_duplicate_vars();

  // The worklist helper expands through scheduler synonyms because they share
  // canonical state, but aliases remain separate scheduling entries.
  CHECK(possible.inSet(variable("cells")));
  CHECK(possible.inSet(variable("cells_shadow")));
  CHECK_FALSE(possible.inSet(variable("cells_alias")));
}

TEST_CASE("duplication metadata and timing models are retained") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(5, 6), {8, 9});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  rule producer = make_rule("upstream", "cells", "model_rule");
  std::map<rule, std::pair<double, double> > comp_info;
  comp_info[producer] = std::make_pair(1.5, 0.25);

  sched.set_proc_able_entities(variable("cells"), producer, make_range(5, 5));
  sched.set_proc_able_entities(variable("cells"), producer, make_range(6, 6));
  sched.set_my_proc_able_entities(variable("cells"), producer,
                                  make_range(5, 5));
  sched.set_reduce_proc_able_entities(variable("cells"), make_range(5, 6));
  sched.set_reduction_outputmap(variable("cells"), true);
  sched.add_model_info(2.0, 0.5, comp_info);
  sched.add_original_computation_time(variable("cells"), 4.0);
  sched.add_original_computation_time(variable("cells"), 1.0);
  sched.add_duplication_computation_time(variable("cells"), 2.5);
  sched.add_original_communication_time(variable("cells"), 3.25);
  sched.add_duplication_communication_time(variable("cells"), 1.75);

  // Verify the duplication bookkeeping helpers retain the stored entity sets
  // and reduction-output-map flag.
  CHECK(sched.get_proc_able_entities(variable("cells"), producer) ==
        make_range(5, 6));
  CHECK(sched.get_my_proc_able_entities(variable("cells"), producer) ==
        make_range(5, 5));
  CHECK(sched.get_reduce_proc_able_entities(variable("cells")) ==
        make_range(5, 6));
  CHECK(sched.is_reduction_outputmap(variable("cells")));

  // Verify the timing accumulators preserve the totals written through the
  // adder helpers.
  CHECK(sched.get_precalculated_original_computation_time(variable("cells")) ==
        doctest::Approx(5.0));
  CHECK(sched.get_precalculated_duplication_computation_time(variable("cells")) ==
        doctest::Approx(2.5));
  CHECK(sched.get_precalculated_original_communication_time(variable("cells")) ==
        doctest::Approx(3.25));
  CHECK(sched.get_precalculated_duplication_communication_time(variable("cells")) ==
        doctest::Approx(1.75));

  // Verify the installed communication and computation models can be read back.
  auto comm_model = sched.get_comm_model();
  double ts = 0.0;
  double tw = 0.0;
  comm_model.get_parameters(ts, tw);
  CHECK(ts == doctest::Approx(2.0));
  CHECK(tw == doctest::Approx(0.5));

  auto computation_model = sched.get_comp_model(producer);
  computation_model.get_parameters(ts, tw);
  CHECK(ts == doctest::Approx(1.5));
  CHECK(tw == doctest::Approx(0.25));

  // Missing rules should report the invalid sentinel model coefficients.
  auto missing_model =
      sched.get_comp_model(make_rule("upstream", "cells", "missing_rule"));
  missing_model.get_parameters(ts, tw);
  CHECK(ts < -99999.0);
  CHECK(tw < -99999.0);
}

TEST_CASE("reduction summaries replace prior values while per-rule sets accumulate") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(1, 4), {10, 11, 12, 13});
  facts.create_fact("cells", cells);

  sched_db sched(facts);
  rule producer = make_rule("upstream", "cells", "reduce_rule");

  sched.set_proc_able_entities(variable("cells"), producer, make_range(1, 2));
  sched.set_proc_able_entities(variable("cells"), producer, make_range(3, 4));
  sched.set_my_proc_able_entities(variable("cells"), producer, make_range(2, 2));
  sched.set_my_proc_able_entities(variable("cells"), producer, make_range(4, 4));
  sched.set_reduce_proc_able_entities(variable("cells"), make_range(1, 2));
  sched.set_reduce_proc_able_entities(variable("cells"), make_range(4, 4));
  sched.set_reduction_outputmap(variable("cells"), true);
  sched.set_reduction_outputmap(variable("cells"), false);

  // The per-rule caches accumulate across writes, while the reduction summary
  // and output-map flag are simple replacement setters.
  CHECK(sched.get_proc_able_entities(variable("cells"), producer) ==
        make_range(1, 4));
  CHECK(sched.get_my_proc_able_entities(variable("cells"), producer) ==
        make_range(2, 2) + make_range(4, 4));
  CHECK(sched.get_reduce_proc_able_entities(variable("cells")) ==
        make_range(4, 4));
  CHECK_FALSE(sched.is_reduction_outputmap(variable("cells")));
}

TEST_CASE("send entity caches only expose variables with mapped-output rules") {
  fact_db facts;
  store<int> cells;
  store<int> plain_cells;
  assign_store_values(cells, make_range(1, 2), {1, 2});
  assign_store_values(plain_cells, make_range(1, 2), {3, 4});
  facts.create_fact("cells", cells);
  facts.create_fact("plain_cells", plain_cells);

  sched_db sched(facts);
  // `parent_map->cells` gives the internal rule descriptor an output mapping.
  rule mapped = make_rule("upstream", "parent_map->cells", "mapped_rule");
  rule plain = make_rule("upstream", "plain_cells", "plain_rule");

  // Mark only one variable as coming from an output-mapped producer.
  sched.set_existential_info(variable("cells"), mapped, make_range(1, 2));
  sched.set_existential_info(variable("plain_cells"), plain, make_range(1, 2));

  // Seed both phase-specific send-entity caches with one mapped and one plain
  // variable.
  std::vector<std::pair<variable, entitySet> > barrier_updates;
  barrier_updates.push_back(
      std::make_pair(variable("cells"), make_range(10, 11)));
  barrier_updates.push_back(
      std::make_pair(variable("plain_cells"), make_range(20, 21)));
  sched.update_send_entities(barrier_updates, sched_db::BARRIER);

  std::vector<std::pair<variable, entitySet> > recurse_updates;
  recurse_updates.push_back(
      std::make_pair(variable("cells"), make_range(30, 31)));
  recurse_updates.push_back(
      std::make_pair(variable("plain_cells"), make_range(40, 41)));
  sched.update_send_entities(recurse_updates, sched_db::RECURSE_PRE);

  variableSet vars =
      make_var_set({variable("cells"), variable("plain_cells")});

  // Verify get_send_entities filters out variables whose producers do not map
  // output entities.
  std::vector<std::pair<variable, entitySet> > barrier =
      sched.get_send_entities(vars, sched_db::BARRIER);
  REQUIRE(barrier.size() == 1);
  CHECK(barrier[0].first == variable("cells"));
  CHECK(barrier[0].second == make_range(10, 11));

  std::vector<std::pair<variable, entitySet> > recurse =
      sched.get_send_entities(vars, sched_db::RECURSE_PRE);
  REQUIRE(recurse.size() == 1);
  CHECK(recurse[0].first == variable("cells"));
  CHECK(recurse[0].second == make_range(30, 31));

  // Rewriting the same cache entry should replace the stored entity set.
  std::vector<std::pair<variable, entitySet> > overwrite;
  overwrite.push_back(std::make_pair(variable("cells"), make_range(50, 50)));
  sched.update_send_entities(overwrite, sched_db::BARRIER);

  barrier = sched.get_send_entities(vars, sched_db::BARRIER);
  REQUIRE(barrier.size() == 1);
  CHECK(barrier[0].second == make_range(50, 50));
}

TEST_CASE("communication list caches return sends before receives in processor order") {
  fact_db facts;
  sched_db sched(facts);

  // Build an intentionally unsorted mix of send and receive records.
  comm_info recv_cells;
  recv_cells.v = variable("cells");
  recv_cells.processor = 2;
  recv_cells.recv_set = sequence(make_range(30, 31));

  comm_info send_cells;
  send_cells.v = variable("cells");
  send_cells.processor = 1;
  send_cells.send_set = make_range(10, 11);

  comm_info recv_flux;
  recv_flux.v = variable("flux");
  recv_flux.processor = 1;
  recv_flux.recv_set = sequence(make_range(40, 40));

  comm_info send_flux;
  send_flux.v = variable("flux");
  send_flux.processor = 0;
  send_flux.send_set = make_range(20, 20);

  std::list<comm_info> first_batch;
  first_batch.push_back(recv_cells);
  first_batch.push_back(send_cells);
  sched.update_comm_info_list(first_batch, sched_db::BARRIER_CLIST);

  std::list<comm_info> second_batch;
  second_batch.push_back(recv_flux);
  second_batch.push_back(send_flux);
  sched.update_comm_info_list(second_batch, sched_db::BARRIER_CLIST);

  variableSet vars = make_var_set({variable("cells"), variable("flux")});
  std::vector<comm_info> sorted =
      to_vector(sched.get_comm_info_list(vars, facts, sched_db::BARRIER_CLIST));

  // Verify repeated cache updates append records, and retrieval repacks the
  // combined cache so sends come first in processor order, followed by receives
  // in processor order.
  REQUIRE(sorted.size() == 4);
  CHECK(sorted[0].processor == 0);
  CHECK(sorted[0].send_set == make_range(20, 20));
  CHECK(sorted[1].processor == 1);
  CHECK(sorted[1].send_set == make_range(10, 11));
  CHECK(sorted[2].processor == 1);
  CHECK(sorted[2].recv_set == sequence(make_range(40, 40)));
  CHECK(sorted[3].processor == 2);
  CHECK(sorted[3].recv_set == sequence(make_range(30, 31)));
}

TEST_CASE("execution sequence cache keeps the first stored sequence") {
  sched_db sched;
  rule producer = make_rule("upstream", "cells", "exec_rule");

  // Missing entries should read back as EMPTY.
  CHECK(sched.get_exec_seq(producer) == EMPTY);

  // The first write populates the cache.
  sched.update_exec_seq(producer, make_range(1, 2));
  CHECK(sched.get_exec_seq(producer) == make_range(1, 2));

  // Later writes are ignored so the original execution sequence stays intact.
  sched.update_exec_seq(producer, make_range(5, 6));
  CHECK(sched.get_exec_seq(producer) == make_range(1, 2));
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
