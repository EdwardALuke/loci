#include <Loci.h>
#include <scheduler.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

using namespace Loci;

namespace {

/// @brief Redirects a stream into an in-memory buffer for the lifetime of the
///   capture object.
class StreamCapture {
  std::ostream &stream_;
  std::streambuf *saved_;
  std::ostringstream captured_;
public:
  explicit StreamCapture(std::ostream &stream)
      : stream_(stream), saved_(stream.rdbuf(captured_.rdbuf())) {}

  ~StreamCapture() { stream_.rdbuf(saved_); }

  /// @brief Returns the text captured so far.
  std::string str() const { return captured_.str(); }
};

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
void assign_store_values(store<int> &values, const entitySet &domain,
                         const std::vector<int> &entries) {
  REQUIRE(static_cast<size_t>(domain.size()) == entries.size());

  values.allocate(domain);
  size_t index = 0;
  for(entitySet::const_iterator ei = domain.begin(); ei != domain.end();
      ++ei, ++index) {
    values[*ei] = entries[index];
  }
}

/// @brief Reads `values` over `domain` into a flat vector.
/// @param[in] values Store to read.
/// @param[in] domain Domain iteration order to use for the readback.
/// @return The values stored across `domain`.
std::vector<int> read_store_values(const store<int> &values,
                                   const entitySet &domain) {
  std::vector<int> result;
  for(entitySet::const_iterator ei = domain.begin(); ei != domain.end(); ++ei) {
    result.push_back(values[*ei]);
  }
  return result;
}

/// @brief Wraps a concrete rule type in the external `rule` handle used by
///   `rule_db`.
template <class TRule>
rule make_rule() {
  return rule(rule_implP(new copy_rule_impl<TRule>));
}

/// @brief Returns true when `rep` wraps a non-null store representation.
/// @param[in] rep Store representation handle to inspect.
/// @return True when `rep` is non-null.
bool has_rep(const storeRepP &rep) {
  return rep != static_cast<storeRep *>(0);
}

/// @brief Returns a single-variable target set for scheduler entry points.
/// @param[in] name Variable name to query or schedule.
/// @return A `variableSet` containing only `name`.
variableSet make_target(const char *name) {
  variableSet target;
  target += variable(name);
  return target;
}

/// @brief Pointwise fixture that copies `src` into `dst` with a small offset.
class offset_copy_rule : public pointwise_rule {
  const_store<int> src;
  store<int> dst;
public:
  offset_copy_rule() {
    name_store("src", src);
    name_store("dst", dst);
    input("src");
    output("dst");
  }

  void calculate(Entity e) { dst[e] = src[e] + 10; }

  void compute(const sequence &seq) { do_loop(seq, this); }
};

} // namespace

TEST_CASE("create_execution_schedule returns null when no rule can produce the target") {
  rule_db rdb;
  fact_db facts;
  sched_db scheds;

  StreamCapture cout_capture(std::cout);
  StreamCapture cerr_capture(std::cerr);
  executeP schedule =
      create_execution_schedule(rdb, facts, scheds, make_target("missing"));

  // Without any supporting facts or rules, scheduler construction should fail
  // cleanly by returning a null schedule.
  CHECK(schedule == executeP(0));
  CHECK_FALSE(facts.get_typed_variables().inSet(variable("missing")));
}

TEST_CASE("create_execution_schedule builds and executes a simple pointwise rule") {
  rule_db rdb;
  fact_db facts;
  sched_db scheds;
  store<int> src;
  const entitySet domain = make_range(4, 6);

  assign_store_values(src, domain, {1, 3, 5});
  facts.create_fact("src", src);
  rdb.add_rule(make_rule<offset_copy_rule>());

  StreamCapture cout_capture(std::cout);
  StreamCapture cerr_capture(std::cerr);
  executeP schedule = create_execution_schedule(rdb, facts, scheds,
                                                make_target("dst"));

  // The scheduler should find a plan for `dst` once `src` and the copy rule
  // are present.
  REQUIRE(schedule != executeP(0));

  schedule->execute(facts, scheds);

  // Executing the schedule should materialize the requested output over the
  // same domain with the rule's transformed values.
  storeRepP rep = facts.get_variable("dst");
  REQUIRE(has_rep(rep));

  store<int> dst(rep);
  CHECK(dst.domain() == domain);
  CHECK(read_store_values(dst, domain) == std::vector<int>({11, 13, 15}));
}

TEST_CASE("makeQuery short-circuits extensional targets without disturbing their facts") {
  rule_db rdb;
  fact_db facts;
  store<int> src;
  const entitySet domain = make_range(8, 10);

  assign_store_values(src, domain, {2, 4, 6});
  facts.create_fact("src", src);

  StreamCapture cout_capture(std::cout);
  StreamCapture cerr_capture(std::cerr);
  const bool ok = makeQuery(rdb, facts, "src");

  // Queries that only name extensional facts should succeed without needing a
  // generated execution schedule.
  CHECK(ok);
  CHECK(facts.get_extensional_facts().inSet(variable("src")));

  // The original extensional store should remain intact after the fast path.
  storeRepP rep = facts.get_variable("src");
  REQUIRE(has_rep(rep));

  store<int> stored(rep);
  CHECK(stored.domain() == domain);
  CHECK(read_store_values(stored, domain) == std::vector<int>({2, 4, 6}));
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
