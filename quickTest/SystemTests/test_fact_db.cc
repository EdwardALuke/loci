#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Loci;

namespace {

/// @brief Builds an inclusive entity range for compact test setup.
entitySet make_range(int first, int last) {
  return entitySet(interval(first, last));
}

/// @brief Allocates `values` on `domain` and fills it from `entries`.
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

/// @brief Wraps serialized facts in the braces expected by `fact_db::read()`.
std::string wrap_fact_db_input(const std::string &body) {
  return "{\n" + body + "\n}";
}

/// @brief Returns true when a fact lookup produced a valid backing rep.
bool has_rep(const storeRepP &rep) {
  return rep != static_cast<storeRep *>(0);
}

/// @brief Returns true when two fact lookups resolve to the same backing rep.
bool same_rep(const storeRepP &lhs, const storeRepP &rhs) {
  if(!has_rep(lhs) || !has_rep(rhs)) {
    return false;
  }
  return lhs.operator->() == rhs.operator->();
}

} // namespace

//----------------------------------------------------------------------------
// Potential Bugs
//----------------------------------------------------------------------------
//
// These are cases that appear to demonstrate behavior that is not intended.
//
// These stay in the suite as visible regressions without breaking the overall
// test run.

// remove_variable() erases the fact and typed-variable entry, but it leaves
// extensional_facts marked with the removed variable name.
TEST_CASE("known bug: remove_variable should clear extensional fact bookkeeping [known-bug]" *
          doctest::may_fail()) {
  fact_db facts;
  param<int> value;
  value = 11;

  // Seed an extensional fact so we can observe every bookkeeping structure
  // that `remove_variable()` is expected to clean up.
  facts.create_fact("doomed", value);
  REQUIRE(facts.get_extensional_facts().inSet(variable("doomed")));

  facts.remove_variable(variable("doomed"));

  // The extensional-fact set should no longer mention the removed variable.
  CHECK_FALSE(facts.get_extensional_facts().inSet(variable("doomed")));
}

TEST_CASE("fact_db constructor seeds EMPTY and UNIVERSE constraint facts") {
  fact_db facts;

  // Fetch the built-in constraint facts that the constructor should install.
  const storeRepP empty_rep = facts.get_variable("EMPTY");
  const storeRepP universe_rep = facts.get_variable("UNIVERSE");
  REQUIRE(has_rep(empty_rep));
  REQUIRE(has_rep(universe_rep));

  Constraint empty_fact(empty_rep);
  Constraint universe_fact(universe_rep);

  // Their values should match the canonical empty and universal sets.
  CHECK(*empty_fact == EMPTY);
  CHECK(*universe_fact == ~EMPTY);

  // The constructor should also record them as extensional facts.
  CHECK(facts.get_extensional_facts().inSet(variable("EMPTY")));
  CHECK(facts.get_extensional_facts().inSet(variable("UNIVERSE")));
}

TEST_CASE("remove_variable erases fact storage and typed-variable lookup") {
  fact_db facts;
  param<int> value;
  value = 11;

  // Seed a regular fact so removal has both storage and type metadata to clear.
  facts.create_fact("doomed", value);
  REQUIRE(facts.get_extensional_facts().inSet(variable("doomed")));

  facts.remove_variable(variable("doomed"));

  // The fact should disappear both from storage lookups and typed-variable
  // bookkeeping.
  CHECK_FALSE(has_rep(facts.get_variable("doomed")));
  CHECK_FALSE(facts.get_typed_variables().inSet(variable("doomed")));
}

TEST_CASE("create_fact installs store facts and updates max allocation") {
  fact_db facts;
  store<int> cells;
  assign_store_values(cells, make_range(5, 7), {10, 20, 30});

  // Install a populated store fact and read it back through the public lookup.
  facts.create_fact("cells", cells);

  const storeRepP fetched_rep = facts.get_variable("cells");
  REQUIRE(has_rep(fetched_rep));
  store<int> fetched(fetched_rep);

  // The fetched store should preserve both the domain and the payload values.
  CHECK(fetched.domain() == make_range(5, 7));
  CHECK(fetched[5] == 10);
  CHECK(fetched[6] == 20);
  CHECK(fetched[7] == 30);

  // `create_fact()` should advance the allocator watermark past the highest key.
  CHECK(facts.get_max_alloc(0) == 8);
}

TEST_CASE("set_variable_type returns fresh empty containers") {
  fact_db facts;
  store<int> field_template;
  assign_store_values(field_template, make_range(1, 2), {4, 9});

  // Register a store type, then ask for two fresh instances of that type.
  facts.set_variable_type("field", field_template);

  storeRepP first = facts.get_variable_type("field");
  storeRepP second = facts.get_variable_type("field");

  REQUIRE(has_rep(first));
  REQUIRE(has_rep(second));
  CHECK(isSTORE(first));
  CHECK(isSTORE(second));
  CHECK_FALSE(same_rep(first, second));

  store<int> first_field(first);
  store<int> second_field(second);

  // Each returned instance should start empty even though the template carried
  // data when it was registered.
  CHECK(first_field.domain() == EMPTY);
  CHECK(second_field.domain() == EMPTY);

  // Mutating one clone should not leak into the other clone.
  assign_store_values(first_field, make_range(9, 10), {1, 2});
  CHECK(second_field.domain() == EMPTY);
}

TEST_CASE("remove_namespace restores unscoped lookup behavior") {
  fact_db facts;
  param<int> rho;
  rho = 17;

  // Create a fact while a namespace is active so the implicit lookup is scoped.
  facts.set_namespace("fluid");
  facts.create_fact("rho", rho);

  const storeRepP scoped_rep = facts.get_variable("rho");
  REQUIRE(has_rep(scoped_rep));
  param<int> via_namespace(scoped_rep);
  CHECK(*via_namespace == 17);

  // After removing the namespace, unqualified lookup should stop resolving and
  // the explicit fully qualified name should still work.
  facts.remove_namespace();
  CHECK_FALSE(has_rep(facts.get_variable("rho")));

  const storeRepP explicit_rep = facts.get_variable("fluid@rho");
  REQUIRE(has_rep(explicit_rep));
  param<int> explicit_lookup(explicit_rep);
  CHECK(*explicit_lookup == 17);
}

TEST_CASE("erase_intensional_facts removes only intensional facts") {
  fact_db facts;
  param<int> extensional;
  param<int> intensional;
  extensional = 3;
  intensional = 5;

  // Seed one fact in each category so we can verify the cleanup is selective.
  facts.create_fact("extensional_value", extensional);
  facts.create_intensional_fact("intensional_value", intensional);

  CHECK(facts.get_extensional_facts().inSet(variable("extensional_value")));
  CHECK_FALSE(
      facts.get_extensional_facts().inSet(variable("intensional_value")));
  CHECK(facts.get_intensional_facts().inSet(variable("intensional_value")));

  facts.erase_intensional_facts();

  // Extensional facts should survive, while intensional facts should be fully
  // removed from storage and type bookkeeping.
  CHECK(has_rep(facts.get_variable("extensional_value")));
  CHECK_FALSE(has_rep(facts.get_variable("intensional_value")));
  CHECK(facts.get_typed_variables().inSet(variable("extensional_value")));
  CHECK_FALSE(facts.get_typed_variables().inSet(variable("intensional_value")));
  CHECK(facts.get_extensional_facts().inSet(variable("extensional_value")));
}

TEST_CASE("synonyms share the same backing fact and are removed with the base") {
  fact_db facts;
  param<int> base_value;
  base_value = 11;

  // Create a synonym and confirm both names appear in the typed-variable view.
  facts.create_fact("base", base_value);
  facts.synonym_variable(variable("base"), variable("alias"));

  variableSet vars = facts.get_typed_variables();
  CHECK(vars.inSet(variable("base")));
  CHECK(vars.inSet(variable("alias")));

  // Writing through the synonym should update the shared underlying fact.
  param<int> alias_value(facts.get_variable("alias"));
  alias_value = 23;

  param<int> fetched_base(facts.get_variable("base"));
  CHECK(*fetched_base == 23);

  // Removing the base should also invalidate the synonym.
  facts.remove_variable(variable("base"));
  CHECK_FALSE(has_rep(facts.get_variable("base")));
  CHECK_FALSE(has_rep(facts.get_variable("alias")));
  CHECK_FALSE(facts.get_typed_variables().inSet(variable("alias")));
}

TEST_CASE("synonym chains resolve to a single canonical fact") {
  fact_db facts;
  param<int> base_value;
  base_value = 31;

  // Build a two-link synonym chain that should still resolve to one canonical
  // stored fact.
  facts.create_fact("base", base_value);
  facts.synonym_variable(variable("base"), variable("alias"));
  facts.synonym_variable(variable("alias"), variable("alias2"));

  // Updating through the tail of the chain should be visible from every name.
  param<int> alias_two(facts.get_variable("alias2"));
  alias_two = 47;

  param<int> fetched_base(facts.get_variable("base"));
  param<int> fetched_alias(facts.get_variable("alias"));
  CHECK(*fetched_base == 47);
  CHECK(*fetched_alias == 47);
  CHECK(facts.get_typed_variables().inSet(variable("alias2")));
}

TEST_CASE("update_fact and replace_fact refresh values without losing fact identity") {
  fact_db facts;
  store<int> original;
  store<int> updated;
  store<int> replacement;
  assign_store_values(original, make_range(5, 7), {10, 20, 30});
  assign_store_values(updated, make_range(5, 6), {1, 2});
  assign_store_values(replacement, make_range(10, 12), {7, 8, 9});

  // Seed the fact, then update it in place with a new domain and payload.
  facts.create_fact("cells", original);
  facts.update_fact("cells", updated);

  store<int> updated_cells(facts.get_variable("cells"));

  // `update_fact()` should replace the payload while preserving the fact entry.
  CHECK(updated_cells.domain() == make_range(5, 6));
  CHECK(updated_cells[5] == 1);
  CHECK(updated_cells[6] == 2);
  CHECK(facts.get_max_alloc(0) == 8);

  // `replace_fact()` should swap in the replacement representation and keep
  // the bookkeeping attached to the same fact name.
  facts.replace_fact("cells", replacement);

  store<int> replaced_cells(facts.get_variable("cells"));
  CHECK(replaced_cells.domain() == make_range(10, 12));
  CHECK(replaced_cells[10] == 7);
  CHECK(replaced_cells[11] == 8);
  CHECK(replaced_cells[12] == 9);
  CHECK(facts.get_typed_variables().inSet(variable("cells")));
  CHECK(facts.get_extensional_facts().inSet(variable("cells")));
  CHECK(facts.get_max_alloc(0) == 13);
}

TEST_CASE("copying fact_db allows replace_fact without mutating the source database") {
  fact_db facts;
  param<int> count;
  count = 1;
  facts.create_fact("count", count);

  // The copy constructor should detach subsequent replacements from the source.
  fact_db copied(facts);
  param<int> copied_count;
  copied_count = 9;
  copied.replace_fact("count", copied_count);

  param<int> original_after_copy(facts.get_variable("count"));
  param<int> copy_after_replace(copied.get_variable("count"));
  CHECK(*original_after_copy == 1);
  CHECK(*copy_after_replace == 9);

  // The assignment operator should provide the same copy-on-write style
  // isolation for later replacements.
  fact_db assigned;
  assigned = facts;
  param<int> assigned_count;
  assigned_count = 27;
  assigned.replace_fact("count", assigned_count);

  param<int> original_after_assignment(facts.get_variable("count"));
  param<int> assigned_after_replace(assigned.get_variable("count"));
  CHECK(*original_after_assignment == 1);
  CHECK(*assigned_after_replace == 27);
}

TEST_CASE("getKeyDomain and serial distributed allocation track per-domain intervals") {
  fact_db facts;

  // Repeated lookup of the same domain name should reuse the same key-domain id.
  CHECK(facts.getKeyDomain("Main") == 0);
  CHECK(facts.getKeyDomain("Main") == 0);

  const int faces_kd = facts.getKeyDomain("Faces");
  CHECK(faces_kd == 1);
  CHECK(facts.getNumKeyDomains() == 2);

  // Distributed allocation should advance monotonically within that domain.
  const std::pair<entitySet, entitySet> first =
      facts.get_distributed_alloc(3, faces_kd);
  CHECK(first.first == make_range(0, 2));
  CHECK(first.second == make_range(0, 2));
  CHECK(facts.get_init_ptn(faces_kd)[0] == make_range(0, 2));

  // The second allocation should continue where the first one stopped.
  const std::pair<entitySet, entitySet> second =
      facts.get_distributed_alloc(2, faces_kd);
  CHECK(second.first == make_range(3, 4));
  CHECK(second.second == make_range(3, 4));
  CHECK(facts.get_init_ptn(faces_kd)[0] == make_range(0, 4));
}

TEST_CASE("rotate_vars cycles fact storage across variables") {
  fact_db facts;
  param<int> a;
  param<int> b;
  param<int> c;
  a = 1;
  b = 2;
  c = 3;

  // Seed three named facts whose values make the rotation easy to track.
  facts.create_fact("a", a);
  facts.create_fact("b", b);
  facts.create_fact("c", c);

  std::list<variable> rotation_order;
  rotation_order.push_back(variable("a"));
  rotation_order.push_back(variable("b"));
  rotation_order.push_back(variable("c"));

  // Rotate the backing reps according to the requested order.
  facts.rotate_vars(rotation_order);

  param<int> rotated_a(facts.get_variable("a"));
  param<int> rotated_b(facts.get_variable("b"));
  param<int> rotated_c(facts.get_variable("c"));

  // The values should shift cyclically with the rotation order.
  CHECK(*rotated_a == 3);
  CHECK(*rotated_b == 1);
  CHECK(*rotated_c == 2);
}

TEST_CASE("read creates facts from registered variable types") {
  fact_db facts;
  param<int> scalar_type;
  store<int> field_type;

  // Register the variable types that `read()` needs in order to instantiate
  // incoming facts.
  facts.set_variable_type("count", scalar_type);
  facts.set_variable_type("field", field_type);

  std::istringstream in(
      "{ count: 7 field: {([2,3]) 11 13} }");

  CHECK_NOTHROW(facts.read(in));

  // Read the newly created facts back through the public lookup API.
  const storeRepP count_rep = facts.get_variable("count");
  const storeRepP field_rep = facts.get_variable("field");
  REQUIRE(has_rep(count_rep));
  REQUIRE(has_rep(field_rep));

  param<int> count(count_rep);
  store<int> field(field_rep);

  // The parsed values should match the serialized input.
  CHECK(*count == 7);
  CHECK(field.domain() == make_range(2, 3));
  CHECK(field[2] == 11);
  CHECK(field[3] == 13);
}

TEST_CASE("read rejects missing fact definitions and malformed envelopes") {
  {
    // Unknown fact names should fail because no variable type was registered.
    fact_db facts;
    std::istringstream in("{ missing: 1 }");
    CHECK_THROWS(facts.read(in));
  }

  {
    // The input envelope itself is required by `fact_db::read()`.
    fact_db facts;
    std::istringstream in("count: 1");
    CHECK_THROWS(facts.read(in));
  }
}

TEST_CASE("write output can be wrapped and read back into a second fact_db") {
  fact_db source;
  param<int> count;
  store<int> field;
  count = 42;
  assign_store_values(field, make_range(3, 4), {8, 13});

  // Serialize a small database containing both scalar and store facts.
  source.create_fact("count", count);
  source.create_fact("field", field);

  std::ostringstream out;
  source.write(out);

  // The serialized form should mention both fact names.
  CHECK(out.str().find("count:") != std::string::npos);
  CHECK(out.str().find("field:") != std::string::npos);

  // Register matching types in the destination database before reading back.
  fact_db copy;
  param<int> scalar_type;
  store<int> field_type;
  copy.set_variable_type("count", scalar_type);
  copy.set_variable_type("field", field_type);

  std::istringstream in(wrap_fact_db_input(out.str()));
  CHECK_NOTHROW(copy.read(in));

  param<int> copied_count(copy.get_variable("count"));
  store<int> copied_field(copy.get_variable("field"));

  // The round-tripped database should reconstruct the original values.
  CHECK(*copied_count == 42);
  CHECK(copied_field.domain() == make_range(3, 4));
  CHECK(copied_field[3] == 8);
  CHECK(copied_field[4] == 13);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
