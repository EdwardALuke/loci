#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <initializer_list>
#include <sstream>
#include <vector>

using namespace Loci;

namespace {

  /*
    entitySet is the public name for intervalSet, while sequence preserves
    execution order, including reverse intervals.  These tests cover the public
    algebra and loop helpers that many scheduler and container paths rely on.
  */

  entitySet set_of(std::initializer_list<int> values) {
    entitySet result;
    for(std::initializer_list<int>::const_iterator i = values.begin();
        i != values.end(); ++i)
      result += *i;
    return result;
  }

  std::vector<int> collect(const entitySet &set) {
    std::vector<int> values;
    for(entitySet::const_iterator i = set.begin(); i != set.end(); ++i)
      values.push_back(*i);
    return values;
  }

  std::vector<int> collect(const sequence &seq) {
    std::vector<int> values;
    for(sequence::const_iterator i = seq.begin(); i != seq.end(); ++i)
      values.push_back(*i);
    return values;
  }

  void check_values(const std::vector<int> &actual,
                    std::initializer_list<int> expected) {
    REQUIRE(actual.size() == expected.size());
    size_t idx = 0;
    for(std::initializer_list<int>::const_iterator i = expected.begin();
        i != expected.end(); ++i, ++idx) {
      CAPTURE(idx);
      CHECK(actual[idx] == *i);
    }
  }

  struct Recorder {
    std::vector<int> seen;
    std::vector<int> segment_starts;
    std::vector<int> segment_counts;
    void calculate(Entity entity) { seen.push_back(entity); }
    void segment(Entity start, Entity count) {
      segment_starts.push_back(start);
      segment_counts.push_back(count);
    }
  };

} // namespace

TEST_CASE("entitySet unions normalize, sort, and coalesce intervals") {
  entitySet set = interval(5, 3);
  set += 7;
  set += interval(6, 6);

  CHECK(set == entitySet(interval(3, 7)));
  CHECK(set.size() == 5);
  CHECK(set.num_intervals() == 1);
  CHECK(set.Min() == 3);
  CHECK(set.Max() == 7);
  CHECK(set.inSet(3));
  CHECK(set.inSet(7));
  CHECK_FALSE(set.inSet(8));
  check_values(collect(set), {3, 4, 5, 6, 7});
}

TEST_CASE("entitySet set algebra handles sparse intervals and complements") {
  entitySet left = entitySet(interval(1, 3)) | interval(8, 10);
  entitySet right = entitySet(interval(3, 8)) | interval(12, 13);

  CHECK((left & right) == set_of({3, 8}));
  CHECK((left - interval(2, 8)) == set_of({1, 9, 10}));
  CHECK(((~EMPTY) & interval(2, 4)) == entitySet(interval(2, 4)));
  CHECK((entitySet(interval(2, 4)) - interval(3, 3)) == set_of({2, 4}));
}

TEST_CASE("entitySet shifts and symmetric difference preserve normalized intervals") {
  entitySet base = entitySet(interval(1, 3)) | interval(6, 8);

  CHECK((base >> 10) ==
        (entitySet(interval(11, 13)) | interval(16, 18)));
  CHECK(((base >> 10) << 10) == base);

  entitySet toggled = base ^ (entitySet(interval(3, 6)) | interval(9, 9));
  CHECK(toggled ==
        (entitySet(interval(1, 2)) | interval(4, 5) | interval(7, 9)));

  toggled ^= entitySet(interval(4, 5));
  CHECK(toggled == (entitySet(interval(1, 2)) | interval(7, 9)));
}

TEST_CASE("create_entitySet compresses adjacent runs from unsorted input") {
  std::vector<int> ids;
  ids.push_back(7);
  ids.push_back(2);
  ids.push_back(3);
  ids.push_back(5);
  ids.push_back(6);

  entitySet set = create_entitySet(ids.begin(), ids.end());

  REQUIRE(set.num_intervals() == 2);
  CHECK(set[0] == interval(2, 3));
  CHECK(set[1] == interval(5, 7));
}

TEST_CASE("entitySet stream format round trips finite and universe intervals") {
  entitySet finite = entitySet(interval(2, 4)) | interval(7, 8);
  std::ostringstream finite_out;
  finite_out << finite;

  entitySet finite_parsed;
  std::istringstream finite_in(finite_out.str());
  finite_in >> finite_parsed;

  CHECK_FALSE(finite_in.fail());
  CHECK(finite_parsed == finite);

  entitySet universe = ~EMPTY;
  std::ostringstream universe_out;
  universe_out << universe;

  entitySet universe_parsed;
  std::istringstream universe_in(universe_out.str());
  universe_in >> universe_parsed;

  CHECK_FALSE(universe_in.fail());
  CHECK(universe_parsed == universe);
}

TEST_CASE("sequence appends preserve reverse intervals during iteration") {
  sequence seq;
  seq += interval(1, 3);
  seq += interval(6, 4);

  CHECK(seq.size() == 6);
  REQUIRE(seq.num_intervals() == 2);
  CHECK(seq[0] == interval(1, 3));
  CHECK(seq[1] == interval(6, 4));
  check_values(collect(seq), {1, 2, 3, 6, 5, 4});

  sequence reversed = seq;
  reversed.Reverse();
  check_values(collect(reversed), {4, 5, 6, 3, 2, 1});

  CHECK(entitySet(seq) == entitySet(interval(1, 6)));
}

TEST_CASE("create_sequence keeps input order and merges adjacent runs") {
  std::vector<int> ids;
  ids.push_back(1);
  ids.push_back(2);
  ids.push_back(5);
  ids.push_back(4);
  ids.push_back(3);

  sequence seq = create_sequence(ids.begin(), ids.end());

  REQUIRE(seq.num_intervals() == 2);
  CHECK(seq[0] == interval(1, 2));
  CHECK(seq[1] == interval(5, 3));
  check_values(collect(seq), {1, 2, 5, 4, 3});
}

TEST_CASE("do_loop variants document sequence ordering semantics") {
  std::vector<int> seen_set;
  do_loop(entitySet(interval(2, 4)),
          [&](Entity entity) { seen_set.push_back(entity); });
  check_values(seen_set, {2, 3, 4});

  sequence seq;
  seq += interval(4, 2);

  std::vector<int> seen_seq;
  do_loop(seq, [&](Entity entity) { seen_seq.push_back(entity); });
  check_values(seen_seq, {4, 3, 2});

  Recorder normalized_recorder;
  do_loop(seq, &normalized_recorder);
  check_values(normalized_recorder.seen, {2, 3, 4});

  Recorder ordered_recorder;
  do_loop_seq(seq, &ordered_recorder);
  check_values(ordered_recorder.seen, {4, 3, 2});

  Recorder segments;
  do_segments(seq, &segments);
  check_values(segments.segment_starts, {2});
  check_values(segments.segment_counts, {3});
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
