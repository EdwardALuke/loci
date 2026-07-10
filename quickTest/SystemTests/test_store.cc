#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

/// Reallocating a store must preserve values for entities shared by the old
/// and new domains.
TEST_CASE("store reallocation preserves values on the overlapping domain") {
  store<int> values;
  values.allocate(interval(1, 4));
  values[1] = 101;
  values[2] = 102;
  values[3] = 103;
  values[4] = 104;

  const entitySet new_domain = interval(3, 6);
  values.allocate(new_domain);

  CHECK(values.domain() == new_domain);
  CHECK(values[3] == 103);
  CHECK(values[4] == 104);
}

/// Remapping a store must move each value to the image of its original entity.
TEST_CASE("store remap carries values to the remapped domain") {
  const entitySet source_domain = interval(1, 3);
  store<int> source;
  source.allocate(source_domain);
  source[1] = 7;
  source[2] = 8;
  source[3] = 9;

  dMap entity_remap;
  entity_remap.allocate(source_domain);
  entity_remap[1] = 12;
  entity_remap[2] = 10;
  entity_remap[3] = 11;

  store<int> remapped(source.Rep()->remap(entity_remap));

  const entitySet remapped_domain = interval(10, 12);
  CHECK(remapped.domain() == remapped_domain);
  CHECK(remapped[10] == 8);
  CHECK(remapped[11] == 9);
  CHECK(remapped[12] == 7);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
