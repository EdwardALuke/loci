#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <algorithm>
#include <initializer_list>
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

  void check_values(const std::vector<int> &actual,
                    std::initializer_list<int> expected) {
    REQUIRE(actual.size() == expected.size());
    size_t idx = 0;
    for(std::initializer_list<int>::const_iterator i = expected.begin();
        i != expected.end(); ++i, ++idx)
      CHECK(actual[idx] == *i);
  }

  void check_row(const multiMap &map, int row,
                 std::initializer_list<int> expected) {
    check_values(row_values(map, row), expected);
  }

} // namespace

TEST_CASE("Map image, preimage, inverse, and compose use the active domains") {
  Map map;
  map.allocate(interval(0, 2));
  map[0] = 10;
  map[1] = 11;
  map[2] = 10;

  CHECK(map.image(set_of({0, 1, 99})) == set_of({10, 11}));

  std::pair<entitySet, entitySet> preimage = map.preimage(set_of({10}));
  CHECK(preimage.first == set_of({0, 2}));
  CHECK(preimage.second == set_of({0, 2}));

  multiMap inverse;
  inverseMap(inverse, map, set_of({10, 11}), map.domain());

  CHECK(inverse.domain() == set_of({10, 11}));
  check_row(inverse, 10, {2, 0});
  check_row(inverse, 11, {1});

  dMap codomain_remap;
  codomain_remap[10] = 100;
  codomain_remap[11] = 101;

  MapRepP(map.Rep())->compose(codomain_remap, map.domain());

  CHECK(map[0] == 100);
  CHECK(map[1] == 101);
  CHECK(map[2] == 100);
}

TEST_CASE("MapRemap moves a Map through separate domain and range remaps") {
  Map map;
  map.allocate(interval(0, 2));
  map[0] = 100;
  map[1] = 101;
  map[2] = 100;

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

TEST_CASE("multiMap image, preimage, inverse, and compose handle row values") {
  store<int> sizes;
  sizes.allocate(interval(0, 2));
  sizes[0] = 2;
  sizes[1] = 1;
  sizes[2] = 0;

  multiMap map(sizes);
  map[0][0] = 10;
  map[0][1] = 11;
  map[1][0] = 10;

  CHECK(MapRepP(map.Rep())->image(interval(0, 2)) == set_of({10, 11}));

  std::pair<entitySet, entitySet> preimage =
    MapRepP(map.Rep())->preimage(set_of({10}));
  CHECK(preimage.first == set_of({1, 2}));
  CHECK(preimage.second == set_of({0, 1}));

  multiMap inverse;
  inverseMap(inverse, map, set_of({10, 11}), map.domain());

  CHECK(inverse.domain() == set_of({10, 11}));
  check_row(inverse, 10, {1, 0});
  check_row(inverse, 11, {0});

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

TEST_CASE("multiMap MapRemap preserves row sizes while remapping keys") {
  store<int> sizes;
  sizes.allocate(interval(0, 2));
  sizes[0] = 2;
  sizes[1] = 1;
  sizes[2] = 0;

  multiMap map(sizes);
  map[0][0] = 10;
  map[0][1] = 11;
  map[1][0] = 10;

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

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
