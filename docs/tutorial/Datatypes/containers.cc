#include <Loci.h>

#include <iostream>
using std::cout ;
using std::endl ;

#include <vector>
using std::vector ;

#ifdef NO_CSTDLIB
#include <string.h>
#else
#include <cstring>
using std::strncpy ;
#endif

// This is an atomic element class.  It is defined here to illustrate that
// blackbox containers can hold arbitrary datatypes.  Loci will make no
// attempt to manage the data stored inside of a blackbox.  This includes
// interprocessor communication as well as serialized input and output.
class Element {
  int mAtomicNumber;
  double mAtomicWeight;
  char mSymbol[3];
public:
  Element() {}
  Element(const int n, const double w, const char * s)
    : mAtomicNumber(n), mAtomicWeight(w) {
    strncpy(mSymbol, s, 2);
  }
};
  
int main(int argc, char *argv[]) {

  Loci::Init(&argc,&argv) ;
  
  // Create a set of entities over which we will contain values using stores
  // Note, the numbering of entities doesn't have to be contiguous; 
  // the Loci store, however, allocates the "gaps", that is it allocates from 
  // Min()
  // to Max(), so it is more memory efficient to allocate in contiguous blocks.
  int alloc_entities[6] = {1,2,3,14,15,16} ;
  entitySet alloc_set = create_entitySet(alloc_entities,alloc_entities+6) ;


  //**************************************************************************
  //* Parameters:
  //* The simplest container is a parameter.  It associates *one* value  to
  //* a set of entities.
  //**************************************************************************

  // To create a parameter use the param<> template with the value type as
  // the template argument.  For example
  param<float> float_param ;

  // The set_entitySet() method assigns the set of entities to which the
  // value will be associated.  By default, the parameter value is associated
  // to all entities.  If a parameter value should only be associated with a
  // subset of entities (for example, a boundary condition) then the
  // set_entitySet method can be used to constrain where the parameter
  // can be applied in a Loci schedule.  For example, to constrain the
  // parameter to only apply to alloc_set entities we perform:
  float_param.set_entitySet(alloc_set) ;

  // There are two possible ways to access the value of a parameter.  One
  // is to use the '*' dereference operator and the other is to use the
  // [] array operator.  The array operator includes a check on debug
  // compiles that ensures that the parameter is defined for the provided
  // entity.  For example, to assign a value to the param use:
  *float_param = 3.1415927 ;

  
  const float radius = 1.0 ;
  // To access a float parameter use the [] operator, for example
  const float area = 2.0*float_param[Entity(1)]*radius ;

  *float_param = area ;
  
  //**************************************************************************
  //* Blackbox:
  //* Works just like a parameter, except it can hold any data type and is not
  //* automatically synchronized between processors.
  //**************************************************************************

  // We cannot pass arguments directly to the constructor of the datatypes
  // stored in the param or blackbox containers since the actual object being
  // instantiated is the container class.  Assuming the copy constructor or
  // assignment operator is valid for the stored datatype, the following
  // syntax may be used to rapidly assign values to the stored object.
  blackbox<Element> hydrogen;

  hydrogen.set_entitySet(alloc_set);

  *hydrogen = Element(1, 1.00794, "H");


  //**************************************************************************
  //* Stores:
  //* The store associates entities to values.
  //**************************************************************************

  // Here we create a container (store) that will contain floats.
  // We can allocate for any type by replacing the type in the template
  // argument.
  store<float> float_store ;

  // We allocate floats for entities labeled in the entitySet alloc_set
  // After this step, we can place values in the store.
  float_store.allocate(alloc_set) ;


  // Here we assign values to the entities in alloc_set and store the
  // resulting values in float_store.  We store consecutively 3/n where n
  // is assigned consecutively to entities starting from 1.
  entitySet::const_iterator ei ;
  int number = 1 ;
  for( ei = alloc_set.begin(); ei != alloc_set.end(); ++ei)
    float_store[*ei] = 3./float(number++) ;

  // Print out contents of store (as pairs of entity identifier and value)
  // Outputs:
  // store contains =  (1,3) (2,1.5) (3,1) (14,0.75) (15,0.6) (16,0.5)
  cout << "store contains = " ;
  for( ei = alloc_set.begin(); ei != alloc_set.end(); ++ei)
    cout << " (" << *ei << "," << float_store[*ei] << ")" ;
  cout << endl ;

  //**************************************************************************
  //* storeVec:
  //* Associate runtime sizeable vectors of values to entities (one size for
  //* all entities.
  //**************************************************************************
  // storeVec has the same interface as the store container with the addition
  // of a setVecSize method.  For example:
  storeVec<float>  vecstore ;
  vecstore.allocate(alloc_set) ;
  int vector_size = 2 ;
  vecstore.setVecSize(vector_size) ; // create 2 length vectors for each entity.

  // Accessing elements of the vector
  number = 1;
  for( ei = alloc_set.begin(); ei != alloc_set.end(); ++ei)
    for(int i=0;i<vector_size;++i) {
      // Accessing vector element i associated with entity *ei
      vecstore[*ei][i] = 3./float(number++) + float(i);
    }
  
  
  //**************************************************************************
  //* storeMat:
  //* Associate runtime sizeable square matrices of values to entities
  // (one size for all entities).
  //**************************************************************************

  storeMat<float> matstore ;
  matstore.allocate(alloc_set) ;
  matstore.setVecSize(vector_size) ; // Create 2x2 matrices
  number = 1;
  for( ei = alloc_set.begin(); ei != alloc_set.end(); ++ei)
    for(int i=0;i<vector_size;++i)
      for(int j=0;j<vector_size;++j) {
        // Accessing vector element i associated with entity *ei
        matstore[*ei][i][j] = 0.0 ;
      }

  //**************************************************************************
  //* multiStore:
  //* Associate runtime sizeable vectors of values to entities
  // (different size for each entity).
  //**************************************************************************
  multiStore<float>  mstore; 

  // To allocate a multiStore, we need a list of the sizes of the vector for
  // each entity.  We use a store<int> to describe these sizes.
  store<int> counts ;
  counts.allocate(alloc_set) ;

  // Create vectors of increasing size
  int size = 1 ;
  for( ei = alloc_set.begin(); ei != alloc_set.end(); ++ei)
    counts[*ei] =  size++ ;

  // Allocate vectors for entities in counts
  mstore.allocate(counts) ;

  // We can then access the vectors of the multiStore in a similar way to
  // storeVec
  for( ei = alloc_set.begin(); ei != alloc_set.end(); ++ei)
    for(int i=0;i<mstore.vec_size(*ei);++i) {
      mstore[*ei][i] = i* mstore.vec_size(*ei) ;
    }


  
  //**************************************************************************
  //* dstore:
  //* A Dynamic Store.  This container uses a hash-table to associate values
  //* with entities.  It has automatic incremental allocation, unlike store,
  //* but also has greater cost of accessing individual values.
  //**************************************************************************

  dstore<float> dynamic_store ;
  number = 1;
  for( ei = alloc_set.begin(); ei != alloc_set.end(); ++ei)
    // Dynamic stores automatically allocate when assigned.
    dynamic_store[*ei] = 3./float(number++) ;

  cout << "dstore contains = " ;
  // Get the entitySet of all allocated entities in dstore
  entitySet domain = dynamic_store.domain() ;
  // loop over all entities in the dynamic store and output.
  for( ei = domain.begin(); ei != domain.end(); ++ei)
    cout << " (" << *ei << "," << dynamic_store[*ei] << ")" ;
  cout << endl ;

  // Containers can contain complex types for more flexible data structures.
  // For example, to create a container that consists of variable numbers of
  // values per entities, use a vector<T> as the type of one of the Loci
  // containers.
  dstore<vector<int> > dynamic_store2 ;

  // Assign values to random entities
  for(int i=0;i<10;++i)
    dynamic_store2[rand()%113].push_back(i) ;
  domain = dynamic_store2.domain() ;
  // Loop over random entities and add another element to even entities
  // which already contain an even number.
  for( ei = domain.begin(); ei != domain.end(); ++ei)
    if((dynamic_store2[*ei].front()%2) == 0)
      dynamic_store2[*ei].push_back(*ei) ;

  cout << "random dstore contains = " ;
  // Get the entitySet of all allocated entities in dstore
  domain = dynamic_store2.domain() ;
  // loop over all entities in the dynamic store and output.

  for( ei = domain.begin(); ei != domain.end(); ++ei) {
    cout << " (" << *ei ;
    for(size_t i=0;i<dynamic_store2[*ei].size();++i) 
      cout << "," << dynamic_store2[*ei][i] ;
    cout << ")" ;
  }
  cout << endl ;

  Loci::Finalize() ;
}
