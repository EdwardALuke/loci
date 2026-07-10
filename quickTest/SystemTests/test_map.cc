#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

/// A Map must carry selected domain entities to their image and recover the
/// domain entities associated with a selected image.
TEST_CASE("Map image and preimage preserve the represented relation") {
  Map face_to_cell;
  const entitySet faces = interval(0, 2);
  face_to_cell.allocate(faces);
  face_to_cell[0] = 10;
  face_to_cell[1] = 11;
  face_to_cell[2] = 10;

  entitySet selected_faces;
  selected_faces += 0;
  selected_faces += 1;
  const entitySet cells = interval(10, 11);
  CHECK(face_to_cell.image(selected_faces) == cells);

  const std::pair<entitySet, entitySet> preimage =
    face_to_cell.preimage(interval(10, 10));
  entitySet faces_of_cell_ten;
  faces_of_cell_ten += 0;
  faces_of_cell_ten += 2;
  CHECK(preimage.first == faces_of_cell_ten);
  CHECK(preimage.second == faces_of_cell_ten);
}

/// inverseMap must group every domain entity under the image entity it reaches.
TEST_CASE("inverseMap groups domain entities by image entity") {
  Map face_to_cell;
  const entitySet faces = interval(0, 2);
  face_to_cell.allocate(faces);
  face_to_cell[0] = 10;
  face_to_cell[1] = 11;
  face_to_cell[2] = 10;

  const entitySet cells = interval(10, 11);
  multiMap cell_to_face;
  inverseMap(cell_to_face, face_to_cell, cells, faces);

  CHECK(cell_to_face.domain() == cells);

  entitySet faces_of_cell_ten;
  for(const Entity *face = cell_to_face.begin(10);
      face != cell_to_face.end(10); ++face)
    faces_of_cell_ten += *face;
  entitySet expected_faces_of_cell_ten;
  expected_faces_of_cell_ten += 0;
  expected_faces_of_cell_ten += 2;
  CHECK(faces_of_cell_ten == expected_faces_of_cell_ten);

  entitySet faces_of_cell_eleven;
  for(const Entity *face = cell_to_face.begin(11);
      face != cell_to_face.end(11); ++face)
    faces_of_cell_eleven += *face;
  const entitySet expected_faces_of_cell_eleven = interval(1, 1);
  CHECK(faces_of_cell_eleven == expected_faces_of_cell_eleven);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
