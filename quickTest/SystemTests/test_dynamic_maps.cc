#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

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

  std::vector<int> static_row(const multiMap &map, int row) {
    return std::vector<int>(map.begin(row), map.end(row));
  }

  std::vector<int> dynamic_row(const dmultiMap &map, int row) {
    const std::vector<int, malloc_alloc<int> > &values = map[row];
    return std::vector<int>(values.begin(), values.end());
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
    check_values(static_row(map, row), expected);
  }

  void check_row(const dmultiMap &map, int row,
                 std::initializer_list<int> expected) {
    check_values(dynamic_row(map, row), expected);
  }

} // namespace

TEST_CASE("dMap insertion, domain edits, freeze, and get_map stay consistent") {
  dMap map;
  map[5] = 50;
  map[6] = 50;
  map[7] = 70;

  CHECK(map.domain() == set_of({5, 6, 7}));
  CHECK(map.image(map.domain()) == set_of({50, 70}));

  std::pair<entitySet, entitySet> preimage = map.preimage(set_of({50}));
  CHECK(preimage.first == set_of({5, 6}));
  CHECK(preimage.second == set_of({5, 6}));

  multiMap inverse;
  inverseMap(inverse, map, set_of({50, 70}), map.domain());
  check_row(inverse, 50, {6, 5});
  check_row(inverse, 70, {7});

  map.erase(set_of({6}));
  CHECK(map.domain() == set_of({5, 7}));

  NPTR<dMapRepI> rep(map.Rep());
  rep->guarantee_domain(set_of({8}));
  map[8] = 80;
  rep->invalidate(set_of({7, 8}));

  CHECK(map.domain() == set_of({7, 8}));
  CHECK(map[7] == 70);
  CHECK(map[8] == 80);

  storeRepP frozen_rep = map.Rep()->freeze();
  Map frozen(frozen_rep);
  CHECK(frozen.domain() == set_of({7, 8}));
  CHECK(frozen[7] == 70);
  CHECK(frozen[8] == 80);

  storeRepP static_rep = MapRepP(map.Rep())->get_map();
  multiMap static_map(static_rep);
  CHECK(static_map.domain() == set_of({7, 8}));
  check_row(static_map, 7, {70});
  check_row(static_map, 8, {80});
}

TEST_CASE("dMap gather, scatter, pack, and remap use dynamic domains") {
  dMap source;
  source[0] = 100;
  source[1] = 101;
  source[2] = 102;

  dMap gather_map;
  gather_map[20] = 2;
  gather_map[21] = 0;
  gather_map[22] = 1;

  dMap gathered;
  gathered.allocate(interval(20, 22));
  storeRepP source_rep = source.Rep();
  gathered.Rep()->gather(gather_map, source_rep, interval(20, 22));

  CHECK(gathered[20] == 102);
  CHECK(gathered[21] == 100);
  CHECK(gathered[22] == 101);

  dMap scattered;
  scattered.allocate(interval(30, 32));

  dMap scatter_map;
  scatter_map[20] = 31;
  scatter_map[21] = 30;
  scatter_map[22] = 32;

  storeRepP gathered_rep = gathered.Rep();
  scattered.Rep()->scatter(scatter_map, gathered_rep, interval(20, 22));

  CHECK(scattered[30] == 100);
  CHECK(scattered[31] == 102);
  CHECK(scattered[32] == 101);

  int size = source.Rep()->pack_size(source.domain());
  std::vector<char> buffer(size);
  int position = 0;
  source.Rep()->pack(&buffer[0], position, size, source.domain());

  dMap unpacked;
  unpacked.allocate(interval(40, 42));
  sequence unpack_order;
  unpack_order += interval(42, 40);

  position = 0;
  unpacked.Rep()->unpack(&buffer[0], position, size, unpack_order);
  CHECK(unpacked[42] == 100);
  CHECK(unpacked[41] == 101);
  CHECK(unpacked[40] == 102);

  dMap domain_remap;
  domain_remap[0] = 50;
  domain_remap[1] = 51;
  domain_remap[2] = 52;

  dMap range_remap;
  range_remap[100] = 7;
  range_remap[102] = 9;

  storeRepP remapped_rep =
    MapRepP(source.Rep())->MapRemap(domain_remap, range_remap);
  dMap remapped(remapped_rep);

  CHECK(remapped.domain() == set_of({50, 52}));
  CHECK(remapped[50] == 7);
  CHECK(remapped[52] == 9);
}

TEST_CASE("dmultiMap insertion, inverse, compose, and get_map preserve rows") {
  dmultiMap map;
  map[0].push_back(10);
  map[0].push_back(11);
  map[1].push_back(10);
  map[2];

  CHECK(map.domain() == set_of({0, 1, 2}));
  CHECK(MapRepP(map.Rep())->image(map.domain()) == set_of({10, 11}));

  std::pair<entitySet, entitySet> preimage =
    MapRepP(map.Rep())->preimage(set_of({10}));
  CHECK(preimage.first == set_of({1, 2}));
  CHECK(preimage.second == set_of({0, 1}));

  dmultiMap inverse;
  inverseMap(inverse, map, set_of({10, 11}), map.domain());

  CHECK(inverse.domain() == set_of({10, 11}));
  check_row(inverse, 10, {0, 1});
  check_row(inverse, 11, {0});

  dMap remap;
  remap[10] = 100;
  remap[11] = 101;
  MapRepP(map.Rep())->compose(remap, set_of({0, 1}));

  check_row(map, 0, {100, 101});
  check_row(map, 1, {100});
  check_row(map, 2, {});

  storeRepP static_rep = MapRepP(map.Rep())->get_map();
  multiMap static_map(static_rep);
  check_row(static_map, 0, {100, 101});
  check_row(static_map, 1, {100});
  check_row(static_map, 2, {});
}

TEST_CASE("dmultiMap domain edits and MapRemap filter incomplete range maps") {
  dmultiMap map;
  map[0].push_back(10);
  map[0].push_back(11);
  map[1].push_back(12);
  map[2];

  NPTR<dmultiMapRepI> rep(map.Rep());
  rep->guarantee_domain(set_of({3}));
  map[3].push_back(13);
  rep->erase(set_of({3}));
  rep->invalidate(set_of({0, 1, 2}));

  CHECK(map.domain() == set_of({0, 1, 2}));
  check_row(map, 2, {});

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

  CHECK(remapped.domain() == set_of({20, 22}));
  check_row(remapped, 20, {100, 101});
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
