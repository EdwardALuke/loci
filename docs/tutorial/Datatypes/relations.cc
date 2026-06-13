#include <Loci.h>

#include <iostream>
using std::cout ;
using std::endl ;

int main(int argc, char *argv[])
{

  Loci::Init(&argc,&argv) ;
  
  // Create a set of entities over which we will contain values using stores
  // Note: the numbering of entities doesn't have to be contiguous, however
  // Loci's store allocates the "gaps", that is it allocates from Min()
  // to Max(), so it is more memory efficient to allocate in contiguous blocks
  entitySet alloc_set = interval(1,6) ;


  //**************************************************************************
  //* Constraints:
  //* A constraint is the simplest form of relation.  It simply asserts that
  //* some subset of entites have a given property.
  //**************************************************************************

  // Here we define a constraint "red" that identifies that the entities
  // numbered '1-3' have a distinct property (i.e. color red)
  constraint red ;
  *red = entitySet(interval(1,3)) ;

  //**************************************************************************
  //* Map:
  //* A Map provides a one-by-one correspondence from one set of entities
  //* to another set of entities
  //**************************************************************************

  // Here we create a map that points to the rightmost entitity in the
  // interval [1,6]
  Map shift_map ;
  shift_map.allocate(alloc_set) ;

  // Here create the relation on the domain [1,6] that shift_map[i] = i+1
  for(int i=1;i<=6;++i)
    shift_map[Entity(i)] = Entity(i+1) ;

  // We can also obtain the domain and range of any given map by using the
  // domain and image member functions. (The range of a map is the image
  // of its domain.) 
  entitySet map_domain = shift_map.domain() ;
  entitySet map_range  = shift_map.image(map_domain) ;
  cout << "The domain of shift map is " << map_domain << endl; 
  cout << "The range of shift map is " << map_range << endl ;


  //**************************************************************************
  //* MapVec:
  //* A MapVec provides a one-to-many correspondence from one set of entities
  //* to another set of entities
  //**************************************************************************

  // Here we use a MapVec to create a set of triangles as follows.  Triangles
  // are numbered [1-6] while nodes are numbered [7-14] as shown in the figure:
  //
  //   7    8    9   10
  //  *____*____*____*
  //  | 1 /| 3 /| 5 /|
  //  | / 2| /4 | /6 |
  //  *____*____*____*
  //  11   12   13   14

  // Allocate triangles [1-6]
  MapVec<3> triangle_nodes ;
  triangle_nodes.allocate(alloc_set) ;

  // Loop over triangles and assign nodes
  for(int i=1;i<=6;++i) {
    int offsets1[3] = {1, 2, 5} ;
    int offsets2[3] = {5, 2, 6} ;
    for(int j=0;j<3;++j)
      if(i%2 == 0) {
        triangle_nodes[Entity(i)][j] = 6+offsets2[j]+i/2-1 ;
      } else {
        triangle_nodes[Entity(i)][j] = 6+offsets1[j]+(i-1)/2 ;
      }
  }

  // Output the domain and range of the resulting triangle_node map
  map_domain = triangle_nodes.domain() ;
  map_range = triangle_nodes.image(map_domain) ;
  cout << "The domain of triangle_nodes is " << map_domain << endl; 
  cout << "The range of triangle_nodes is " << map_range << endl ;

  
  //**************************************************************************
  //* multiMap:
  //* A MapVec provides a one-to-many correspondence from one set of entities
  //* to another set of entities where the number of relations vary in size
  //* from entity to entity.  This relation parallels the multiStore
  //* container.
  //**************************************************************************

  // Here we create a CRS for an upper triangular matrix over entities [1,6]
  multiMap upper ;

  store<int> count ;
  count.allocate(alloc_set) ;
  for(int i=1;i<=6;++i) {
    count[Entity(i)] = 7-i ;
  }
  upper.allocate(count) ;
  entitySet::const_iterator ei ;
  for(ei=alloc_set.begin(); ei!= alloc_set.end(); ++ei) {
    for(int i=0;i<upper.vec_size(*ei);++i)
      upper[*ei][i] = i ;
  }

  // Output MultiMap
  cout << "upper = " << endl ;
  for(ei=alloc_set.begin(); ei!= alloc_set.end(); ++ei) {
    cout << *ei ;
    for(int i=0;i<upper.vec_size(*ei);++i)
      cout << " " << upper[*ei][i] ;
    cout << endl;
  }
  
  // There are dynamic hash based coutnerparts to all of the above Maps
  //  dMap dynamic_map ;
  //  dMapVec<3> dynamic_map3 ;
  //  dmultiMap dynamic_multiMap ;

  dMap d1 ;

  for(int i=0;i<10;++i) {
    int i2 = 1000-i ;
    int i3 = 9-i ;
    d1[-i] = i2 ;
    d1[i2] = -i3 ;
  }
  cout << "d1 = " << d1 << endl ;
  entitySet dom = d1.domain() ;
  cout << "d1.domain() = " << dom << endl ;

  Loci::storeRepP rp = d1.Rep() ;
  Loci::MapRepP mp = Loci::MapRepP(rp) ;

  entitySet image = mp->image(dom) ;
  std::pair<entitySet,entitySet> pimage = mp->preimage(dom) ;
  cout << "d1 image(domain) = " << image << endl ;
  cout << "d1 preimage(image(domain)) = " << pimage.first << endl ;
  dMap d2 ;

  entitySet map_entities ;
  map_entities += interval(-151,70735) ;
  map_entities += interval(70886,106292) ;
  
  for(entitySet::const_iterator ei = map_entities.begin();ei!=map_entities.end();++ei)
    d2[*ei] = *ei ;

  cout << "d2.domain() =" << d2.domain() << endl ;
  //  Map d3 ;
  //  d3 = mp->remap(d2) ;
  //  cout << "d3 = " << d3 << endl ;
  constraint x ;
  x = ~EMPTY ;
  constraint y ;
  y = x.Rep()->remap(d2) ;
  cout << "y.dom = " << y.Rep()->domain() << endl ;

  Loci::Finalize() ;

  return 0 ;
}
