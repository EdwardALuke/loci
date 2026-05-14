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

  void fill(store<int> &values, int offset) {
    entitySet dom = values.domain();
    for(entitySet::const_iterator i = dom.begin(); i != dom.end(); ++i)
      values[*i] = offset + *i;
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

  store<int> target;
  target.allocate(interval(10, 12));

  sequence unpack_order;
  unpack_order += interval(12, 10);

  position = 0;
  target.Rep()->unpack(&buffer[0], position, size, unpack_order);

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
