#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <initializer_list>
#include <sstream>
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

  void fill(store<int> &values, int offset) {
    entitySet dom = values.domain();
    for(entitySet::const_iterator i = dom.begin(); i != dom.end(); ++i)
      values[*i] = offset + *i;
  }

  std::vector<int> read_values(const store<int> &values,
                               const entitySet &domain) {
    std::vector<int> result;
    for(entitySet::const_iterator i = domain.begin(); i != domain.end(); ++i)
      result.push_back(values[*i]);
    return result;
  }

} // namespace

TEST_CASE("store allocates sparse domains and exposes indexed values") {
  entitySet dom = entitySet(interval(2, 4)) | interval(8, 8);

  store<int> values;
  values.allocate(dom);
  fill(values, 100);

  CHECK(values.domain() == dom);
  CHECK(values[2] == 102);
  CHECK(values[4] == 104);
  CHECK(values[8] == 108);
}

TEST_CASE("store reallocation preserves values on the overlapping domain") {
  store<int> values;
  values.allocate(interval(1, 4));
  fill(values, 100);

  values.allocate(interval(3, 6));

  CHECK(values.domain() == entitySet(interval(3, 6)));
  CHECK(values[3] == 103);
  CHECK(values[4] == 104);

  values[5] = 205;
  values[6] = 206;
  values.allocate(interval(4, 5));

  CHECK(values.domain() == entitySet(interval(4, 5)));
  CHECK(values[4] == 104);
  CHECK(values[5] == 205);

  values.allocate(EMPTY);
  CHECK(values.domain() == EMPTY);

  values.allocate(interval(9, 9));
  values[9] = 909;
  CHECK(values.domain() == entitySet(interval(9, 9)));
  CHECK(values[9] == 909);
}

TEST_CASE("store shift renumbers the domain without moving stored values") {
  store<int> values;
  values.allocate(set_of({2, 4}));
  values[2] = 20;
  values[4] = 40;

  values.Rep()->shift(10);

  CHECK(values.domain() == set_of({12, 14}));
  CHECK(values[12] == 20);
  CHECK(values[14] == 40);
}

TEST_CASE("store new_store creates an independent empty store of the same type") {
  store<int> source;
  source.allocate(interval(1, 2));
  fill(source, 10);

  storeRepP fresh(source.Rep()->new_store(interval(7, 8)));
  store<int> created(fresh);
  fill(created, 100);

  CHECK(isSTORE(fresh));
  CHECK(created.domain() == entitySet(interval(7, 8)));
  CHECK(read_values(created, interval(7, 8)) == std::vector<int>({107, 108}));
  CHECK(read_values(source, interval(1, 2)) == std::vector<int>({11, 12}));
}

TEST_CASE("storeRep copy updates only the requested context") {
  store<int> source;
  source.allocate(interval(0, 4));
  fill(source, 10);

  store<int> target;
  target.allocate(interval(0, 4));
  fill(target, -100);

  storeRepP source_rep = source.Rep();
  target.Rep()->copy(source_rep, set_of({1, 3}));

  CHECK(target[0] == -100);
  CHECK(target[1] == 11);
  CHECK(target[2] == -98);
  CHECK(target[3] == 13);
  CHECK(target[4] == -96);
}

TEST_CASE("store pack_size reports the packable intersection") {
  store<int> source;
  source.allocate(interval(2, 4));
  source[2] = 22;
  source[3] = 33;
  source[4] = 44;

  entitySet packed;
  const entitySet requested = set_of({1, 2, 4, 9});
  int size = source.Rep()->pack_size(requested, packed);

  CHECK(packed == set_of({2, 4}));
  CHECK(size == static_cast<int>(2 * sizeof(int)));

  std::vector<char> buffer(size);
  int position = 0;
  source.Rep()->pack(&buffer[0], position, size, packed);
  CHECK(position == size);

  store<int> target;
  target.allocate(interval(20, 21));
  sequence unpack_order;
  unpack_order += interval(20, 21);

  position = 0;
  target.Rep()->unpack(&buffer[0], position, size, unpack_order);
  CHECK(position == size);

  CHECK(read_values(target, interval(20, 21)) == std::vector<int>({22, 44}));
}

TEST_CASE("store one-argument pack_size ignores out-of-domain entities" *
          doctest::may_fail()) {
  store<int> source;
  source.allocate(interval(2, 4));

  const entitySet requested = set_of({1, 2, 4, 9});

  CHECK(source.Rep()->pack_size(requested) ==
        static_cast<int>(2 * sizeof(int)));
}

TEST_CASE("store gather and scatter move values through dMap images") {
  store<int> source;
  source.allocate(interval(10, 12));
  source[10] = 1000;
  source[11] = 1100;
  source[12] = 1200;

  dMap gather_map;
  gather_map.allocate(interval(0, 2));
  gather_map[0] = 12;
  gather_map[1] = 10;
  gather_map[2] = 11;

  store<int> gathered;
  gathered.allocate(interval(0, 2));
  storeRepP source_rep = source.Rep();
  gathered.Rep()->gather(gather_map, source_rep, interval(0, 2));

  CHECK(gathered[0] == 1200);
  CHECK(gathered[1] == 1000);
  CHECK(gathered[2] == 1100);

  store<int> scattered;
  scattered.allocate(interval(20, 22));

  dMap scatter_map;
  scatter_map.allocate(interval(0, 2));
  scatter_map[0] = 21;
  scatter_map[1] = 20;
  scatter_map[2] = 22;

  storeRepP gathered_rep = gathered.Rep();
  scattered.Rep()->scatter(scatter_map, gathered_rep, interval(0, 2));

  CHECK(scattered[20] == 1000);
  CHECK(scattered[21] == 1200);
  CHECK(scattered[22] == 1100);
}

TEST_CASE("store remap creates a new store over the map image") {
  store<int> source;
  source.allocate(interval(1, 3));
  source[1] = 7;
  source[2] = 8;
  source[3] = 9;

  dMap remap;
  remap.allocate(interval(1, 3));
  remap[1] = 30;
  remap[2] = 20;
  remap[3] = 31;

  store<int> result(source.Rep()->remap(remap));

  CHECK(result.domain() == set_of({20, 30, 31}));
  CHECK(result[20] == 8);
  CHECK(result[30] == 7);
  CHECK(result[31] == 9);
}

TEST_CASE("store stream output can be read back into another store") {
  store<int> source;
  source.allocate(set_of({3, 5, 6}));
  source[3] = 13;
  source[5] = 15;
  source[6] = 16;

  std::ostringstream output;
  source.Print(output);

  store<int> parsed;
  std::istringstream input(output.str());
  parsed.Input(input);

  CHECK(parsed.domain() == source.domain());
  CHECK(read_values(parsed, parsed.domain()) == std::vector<int>({13, 15, 16}));
}

TEST_CASE("store pack and unpack respects sequence ordering") {
  store<int> source;
  source.allocate(interval(2, 4));
  source[2] = 22;
  source[3] = 33;
  source[4] = 44;

  int size = source.Rep()->pack_size(interval(2, 4));
  std::vector<char> buffer(size);
  int position = 0;
  source.Rep()->pack(&buffer[0], position, size, interval(2, 4));
  CHECK(position == size);

  store<int> target;
  target.allocate(interval(10, 12));

  sequence unpack_order;
  unpack_order += interval(12, 10);

  position = 0;
  target.Rep()->unpack(&buffer[0], position, size, unpack_order);
  CHECK(position == size);

  CHECK(target[12] == 22);
  CHECK(target[11] == 33);
  CHECK(target[10] == 44);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
