#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <algorithm>
#include <initializer_list>
#include <sstream>
#include <utility>
#include <vector>

using namespace Loci;

namespace {

  entitySet set_of(std::initializer_list<int> values) {
    entitySet result;
    for(std::initializer_list<int>::const_iterator i = values.begin();
        i != values.end(); ++i)
      result += *i;
    return result;
  }

  std::vector<int> row_values(const multiMap &map, int row) {
    return std::vector<int>(map.begin(row), map.end(row));
  }

  void check_vector_values(const std::vector<int> &actual,
                           const std::vector<int> &expected) {
    REQUIRE(actual.size() == expected.size());
    for(size_t idx = 0; idx < expected.size(); ++idx) {
      CAPTURE(idx);
      CHECK(actual[idx] == expected[idx]);
    }
  }

  void check_values(const std::vector<int> &actual,
                    std::initializer_list<int> expected) {
    check_vector_values(actual, std::vector<int>(expected));
  }

  void check_values_unordered(std::vector<int> actual,
                              std::initializer_list<int> expected) {
    std::vector<int> sorted_expected(expected);
    std::sort(actual.begin(), actual.end());
    std::sort(sorted_expected.begin(), sorted_expected.end());
    check_vector_values(actual, sorted_expected);
  }

  void check_row(const multiMap &map, int row,
                 std::initializer_list<int> expected) {
    check_values(row_values(map, row), expected);
  }

  void check_row_unordered(const multiMap &map, int row,
                           std::initializer_list<int> expected) {
    check_values_unordered(row_values(map, row), expected);
  }

  void fill_map(Map &map,
                std::initializer_list<std::pair<int, int> > entries) {
    entitySet domain;
    for(std::initializer_list<std::pair<int, int> >::const_iterator i =
          entries.begin(); i != entries.end(); ++i)
      domain += i->first;

    map.allocate(domain);
    for(std::initializer_list<std::pair<int, int> >::const_iterator i =
          entries.begin(); i != entries.end(); ++i)
      map[i->first] = i->second;
  }

  void fill_row_sizes(store<int> &sizes,
                      std::initializer_list<std::pair<int, int> > rows) {
    entitySet domain;
    for(std::initializer_list<std::pair<int, int> >::const_iterator i =
          rows.begin(); i != rows.end(); ++i)
      domain += i->first;

    sizes.allocate(domain);
    for(std::initializer_list<std::pair<int, int> >::const_iterator i =
          rows.begin(); i != rows.end(); ++i)
      sizes[i->first] = i->second;
  }

  void fill_sample_multimap(multiMap &map) {
    store<int> sizes;
    fill_row_sizes(sizes, {{0, 2}, {1, 1}, {2, 0}});

    // One shared fixture covers duplicate images and an intentionally empty row.
    map.allocate(sizes);
    map[0][0] = 10;
    map[0][1] = 11;
    map[1][0] = 10;
  }

} // namespace

TEST_CASE("image_section returns unique mapped values") {
  const int dense_values[] = {3, 5, 4, 5};
  const int sparse_values[] = {-10, 100, -10};

  CHECK(image_section(dense_values, dense_values + 4) == set_of({3, 4, 5}));
  CHECK(image_section(sparse_values, sparse_values + 3) == set_of({-10, 100}));
  CHECK(image_section(dense_values, dense_values) == EMPTY);
}

TEST_CASE("Map image, preimage, inverse, and compose use the active domains") {
  Map map;
  fill_map(map, {{0, 10}, {1, 11}, {2, 10}});

  CHECK(map.image(set_of({0, 1, 99})) == set_of({10, 11}));

  std::pair<entitySet, entitySet> preimage = map.preimage(set_of({10}));
  CHECK(preimage.first == set_of({0, 2}));
  CHECK(preimage.second == set_of({0, 2}));

  multiMap inverse;
  inverseMap(inverse, map, set_of({10, 11}), map.domain());

  CHECK(inverse.domain() == set_of({10, 11}));
  check_row_unordered(inverse, 10, {0, 2});
  check_row_unordered(inverse, 11, {1});

  dMap codomain_remap;
  codomain_remap[10] = 100;
  codomain_remap[11] = 101;

  MapRepP(map.Rep())->compose(codomain_remap, map.domain());

  CHECK(map[0] == 100);
  CHECK(map[1] == 101);
  CHECK(map[2] == 100);
}

TEST_CASE("Map and multiMap image handle sparse multi-interval domains") {
  Map map;
  fill_map(map, {{0, 10}, {2, 20}, {4, 10}, {6, 30}, {8, 20}});

  CHECK(map.image(set_of({0, 2, 4, 6, 8})) == set_of({10, 20, 30}));

  store<int> sizes;
  fill_row_sizes(sizes, {{0, 1}, {2, 1}, {4, 1}, {6, 1}, {8, 1}});

  multiMap row_map;
  row_map.allocate(sizes);
  row_map[0][0] = 10;
  row_map[2][0] = 20;
  row_map[4][0] = 10;
  row_map[6][0] = 30;
  row_map[8][0] = 20;

  CHECK(MapRepP(row_map.Rep())->image(set_of({0, 2, 4, 6, 8})) ==
        set_of({10, 20, 30}));
}

TEST_CASE("inverseMap keeps requested image rows and filters preimage rows") {
  Map map;
  fill_map(map, {{0, 10}, {1, 10}, {2, 11}});

  multiMap inverse;
  inverseMap(inverse, map, set_of({10, 11, 12}), set_of({0, 2, 99}));

  CHECK(inverse.domain() == set_of({10, 11, 12}));
  check_row_unordered(inverse, 10, {0});
  check_row_unordered(inverse, 11, {2});
  check_row_unordered(inverse, 12, {});

  multiMap row_map;
  fill_sample_multimap(row_map);

  multiMap row_inverse;
  inverseMap(row_inverse, row_map, set_of({10, 11, 12}), set_of({0}));

  CHECK(row_inverse.domain() == set_of({10, 11, 12}));
  check_row_unordered(row_inverse, 10, {0});
  check_row_unordered(row_inverse, 11, {0});
  check_row_unordered(row_inverse, 12, {});
}

TEST_CASE("Map stream round trip preserves sparse domains and values") {
  Map map;
  fill_map(map, {{-2, 20}, {0, 30}, {3, 20}});

  std::stringstream stream;
  stream << map;

  Map parsed;
  stream >> parsed;

  CHECK(parsed.domain() == set_of({-2, 0, 3}));
  CHECK(parsed[-2] == 20);
  CHECK(parsed[0] == 30);
  CHECK(parsed[3] == 20);
}

TEST_CASE("MapRemap moves a Map through separate domain and range remaps") {
  Map map;
  fill_map(map, {{0, 100}, {1, 101}, {2, 100}});

  dMap domain_remap;
  domain_remap[0] = 20;
  domain_remap[1] = 21;
  domain_remap[2] = 22;

  dMap range_remap;
  range_remap[100] = 7;
  range_remap[101] = 8;

  storeRepP remapped_rep =
    MapRepP(map.Rep())->MapRemap(domain_remap, range_remap);
  Map remapped(remapped_rep);

  CHECK(remapped.domain() == set_of({20, 21, 22}));
  CHECK(remapped[20] == 7);
  CHECK(remapped[21] == 8);
  CHECK(remapped[22] == 7);
}

TEST_CASE("MapRemap drops Map entries without a range remap") {
  Map map;
  fill_map(map, {{0, 100}, {1, 999}, {2, 101}});

  dMap domain_remap;
  domain_remap[0] = 20;
  domain_remap[1] = 21;
  domain_remap[2] = 22;

  dMap range_remap;
  range_remap[100] = 7;
  range_remap[101] = 8;

  storeRepP remapped_rep =
    MapRepP(map.Rep())->MapRemap(domain_remap, range_remap);
  Map remapped(remapped_rep);

  CHECK(remapped.domain() == set_of({20, 22}));
  CHECK(remapped[20] == 7);
  CHECK(remapped[22] == 8);
}

TEST_CASE("multiMap image, preimage, inverse, and compose handle row values") {
  multiMap map;
  fill_sample_multimap(map);

  CHECK(MapRepP(map.Rep())->image(interval(0, 2)) == set_of({10, 11}));

  std::pair<entitySet, entitySet> preimage =
    MapRepP(map.Rep())->preimage(set_of({10}));
  // An empty row is vacuously contained in the codomain, but does not overlap it.
  CHECK(preimage.first == set_of({1, 2}));
  CHECK(preimage.second == set_of({0, 1}));

  multiMap inverse;
  inverseMap(inverse, map, set_of({10, 11}), map.domain());

  CHECK(inverse.domain() == set_of({10, 11}));
  check_row_unordered(inverse, 10, {0, 1});
  check_row_unordered(inverse, 11, {0});

  dMap codomain_remap;
  codomain_remap[10] = 100;
  codomain_remap[11] = 101;

  MapRepP(map.Rep())->compose(codomain_remap, set_of({0, 1}));
  check_row(map, 0, {100, 101});
  check_row(map, 1, {100});
  check_row(map, 2, {});

  dMap partial_remap;
  partial_remap[100] = 1000;

  MapRepP(map.Rep())->compose(partial_remap, set_of({0}));
  check_row(map, 0, {1000, -1});
}

TEST_CASE("multiMap stream round trip preserves row sizes and values") {
  multiMap map;
  fill_sample_multimap(map);

  std::stringstream stream;
  stream << map;

  multiMap parsed;
  stream >> parsed;

  CHECK(parsed.domain() == set_of({0, 1, 2}));
  check_row(parsed, 0, {10, 11});
  check_row(parsed, 1, {10});
  check_row(parsed, 2, {});
}

TEST_CASE("multiMap MapRemap preserves row sizes while remapping keys") {
  multiMap map;
  fill_sample_multimap(map);

  dMap domain_remap;
  domain_remap[0] = 20;
  domain_remap[1] = 21;
  domain_remap[2] = 22;

  dMap range_remap;
  range_remap[10] = 100;
  range_remap[11] = 101;

  storeRepP remapped_rep =
    MapRepP(map.Rep())->MapRemap(domain_remap, range_remap);
  multiMap remapped(remapped_rep);

  CHECK(remapped.domain() == set_of({20, 21, 22}));
  check_row(remapped, 20, {100, 101});
  check_row(remapped, 21, {100});
  check_row(remapped, 22, {});
}

TEST_CASE("multiMap MapRemap marks row values missing from the range remap") {
  multiMap map;
  fill_sample_multimap(map);

  dMap domain_remap;
  domain_remap[0] = 20;
  domain_remap[1] = 21;
  domain_remap[2] = 22;

  dMap range_remap;
  range_remap[10] = 100;

  storeRepP remapped_rep =
    MapRepP(map.Rep())->MapRemap(domain_remap, range_remap);
  multiMap remapped(remapped_rep);

  CHECK(remapped.domain() == set_of({20, 21, 22}));
  check_row(remapped, 20, {100, -1});
  check_row(remapped, 21, {100});
  check_row(remapped, 22, {});
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
