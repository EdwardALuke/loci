#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <list>
#include <sstream>
#include <string>

using namespace Loci;

namespace {

/// @brief Builds an inclusive entity range for compact test setup.
entitySet make_range(int first, int last) {
  return entitySet(interval(first, last));
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

template <typename Fn>
/// @brief Captures `std::cerr` while `fn` runs so expected failures stay quiet.
std::string capture_cerr(Fn fn) {
  std::ostringstream err;
  std::streambuf *old = std::cerr.rdbuf(err.rdbuf());
  try {
    fn();
  } catch(...) {
    std::cerr.rdbuf(old);
    throw;
  }
  std::cerr.rdbuf(old);
  return err.str();
}

} // namespace

//----------------------------------------------------------------------------
// Fact Database Behavior
//----------------------------------------------------------------------------
// These tests focus on public fact_db behavior: which facts exist, what values
// they expose, how namespaces and synonyms affect lookup, and whether read/write
// round trips preserve user-visible facts.
//
// Known-bug tests stay here when they describe a public fact-db invariant that
// should be fixed later, not just an internal implementation detail.

TEST_CASE("known bug: remove_variable should clear extensional fact membership [known-bug]" *
          doctest::may_fail()) {
  fact_db facts;
  param<int> value;
  value = 11;

  // `create_fact()` marks user-provided inputs as extensional facts, which is
  // the public view used to distinguish supplied data from derived data.
  facts.create_fact("doomed", value);
  REQUIRE(facts.get_extensional_facts().inSet(variable("doomed")));

  facts.remove_variable(variable("doomed"));

  // Once the fact is removed, it should no longer appear in any public fact set.
  CHECK_FALSE(has_rep(facts.get_variable("doomed")));
  CHECK_FALSE(facts.get_typed_variables().inSet(variable("doomed")));
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

  // Seed a regular fact so removal has lookup state and a typed-variable entry
  // to clear.
  facts.create_fact("doomed", value);
  REQUIRE(facts.get_extensional_facts().inSet(variable("doomed")));

  facts.remove_variable(variable("doomed"));

  // The fact should disappear from both lookup and the typed-variable view.
  CHECK_FALSE(has_rep(facts.get_variable("doomed")));
  CHECK_FALSE(facts.get_typed_variables().inSet(variable("doomed")));
}

TEST_CASE("create_fact installs scalar facts and marks them extensional") {
  fact_db facts;
  param<int> count;
  count = 30;

  // Install a user-provided fact and read it back through public lookup.
  facts.create_fact("count", count);

  const storeRepP fetched_rep = facts.get_variable("count");
  REQUIRE(has_rep(fetched_rep));
  param<int> fetched(fetched_rep);

  // The value should be available to query setup as an extensional input fact.
  CHECK(*fetched == 30);
  CHECK(facts.get_typed_variables().inSet(variable("count")));
  CHECK(facts.get_extensional_facts().inSet(variable("count")));
}

TEST_CASE("set_variable_type returns fresh containers for read-created facts") {
  fact_db facts;
  param<int> count_template;
  count_template = 4;

  // Register a type, then ask for two fresh instances of that type.
  facts.set_variable_type("count", count_template);

  storeRepP first = facts.get_variable_type("count");
  storeRepP second = facts.get_variable_type("count");

  REQUIRE(has_rep(first));
  REQUIRE(has_rep(second));
  CHECK(isPARAMETER(first));
  CHECK(isPARAMETER(second));
  CHECK_FALSE(same_rep(first, second));
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

  // Extensional facts should survive, while intensional facts should disappear
  // from lookup and the typed-variable view.
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
  param<int> original;
  param<int> updated;
  param<int> replacement;
  original = 10;
  updated = 20;
  replacement = 30;

  // Seed the fact, then update it in place with a new payload.
  facts.create_fact("count", original);
  facts.update_fact("count", updated);

  param<int> updated_count(facts.get_variable("count"));

  // `update_fact()` should replace the payload while preserving the fact entry.
  CHECK(*updated_count == 20);
  CHECK(facts.get_typed_variables().inSet(variable("count")));

  // `replace_fact()` should swap in the replacement representation and keep
  // the fact visible under the same name.
  facts.replace_fact("count", replacement);

  param<int> replaced_count(facts.get_variable("count"));
  CHECK(*replaced_count == 30);
  CHECK(facts.get_typed_variables().inSet(variable("count")));
  CHECK(facts.get_extensional_facts().inSet(variable("count")));
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

  // Register the variable types that `read()` needs in order to instantiate
  // incoming facts.
  facts.set_variable_type("count", scalar_type);

  std::istringstream in("{ count: 7 }");

  CHECK_NOTHROW(facts.read(in));

  // Read the newly created facts back through the public lookup API.
  const storeRepP count_rep = facts.get_variable("count");
  REQUIRE(has_rep(count_rep));

  param<int> count(count_rep);

  // The parsed values should match the serialized input.
  CHECK(*count == 7);
}

TEST_CASE("read rejects missing fact definitions and malformed envelopes") {
  {
    // Unknown fact names should fail because no variable type was registered.
    fact_db facts;
    std::istringstream in("{ missing: 1 }");
    CHECK_THROWS(capture_cerr([&]() { facts.read(in); }));
  }

  {
    // The input envelope itself is required by `fact_db::read()`.
    fact_db facts;
    std::istringstream in("count: 1");
    CHECK_THROWS(capture_cerr([&]() { facts.read(in); }));
  }
}

TEST_CASE("write output can be wrapped and read back into a second fact_db") {
  fact_db source;
  param<int> count;
  count = 42;

  // Serialize a small database containing a user-provided scalar fact.
  source.create_fact("count", count);

  std::ostringstream out;
  source.write(out);

  // The serialized form should mention the fact name.
  CHECK(out.str().find("count:") != std::string::npos);

  // Register matching types in the destination database before reading back.
  fact_db copy;
  param<int> scalar_type;
  copy.set_variable_type("count", scalar_type);

  std::istringstream in(wrap_fact_db_input(out.str()));
  CHECK_NOTHROW(copy.read(in));

  param<int> copied_count(copy.get_variable("count"));

  // The round-tripped database should reconstruct the original values.
  CHECK(*copied_count == 42);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
