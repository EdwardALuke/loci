// Example program illustrating how to manipulate Entity, entitySet,
// and sequence data structures in Loci

// Every Loci program includes the Loci header file.
#include <Loci.h>

// includes for standard C++ functions
#include <iostream>
using std::cout ;
using std::endl ;
#include <vector>
using std::vector ;

int main(int argc,char *argv[])
{

  Loci::Init(&argc,&argv) ;
  
  ////////////////////////////////////////////////////////////////////////////
  // The Entity class 
  ////////////////////////////////////////////////////////////////////////////
  // An Entity is labeled by an integer.
  // The collections of entities discussed below are represented by
  // collections of integers.
  Entity e10 = Entity(10) ;  // First, the entity labeled by integer 10
  cout << "e10 = " << e10 << endl ;

  ////////////////////////////////////////////////////////////////////////////
  // The interval (A collection of consecutively labeled entities)
  ////////////////////////////////////////////////////////////////////////////
  // For example, onedigit is an interval of entity labels consisting of
  // only one decimal digit.
  interval onedigit = interval(0,9) ;
  cout << "onedigit = " << "[0,9]" << endl ; 

  ////////////////////////////////////////////////////////////////////////////
  // Class entitySet provides facilities for general sets of entities.
  ////////////////////////////////////////////////////////////////////////////
  // Initialization
  // An entitySet can be initialized to an interval.
  entitySet A = onedigit ;
  entitySet B = interval(14,100) ;
  entitySet C = interval(5,15) ;
  entitySet F = interval(10,20) ;
  cout << "A = " << A << endl ; 
  cout << "B = " << B << endl ; 
  cout << "C = " << C << endl ;
  cout << "F = " << F << endl ; 


  ////////////////////////////////////////////////////////////////////////////
  // For efficiency, an entitySet is stored as an ordered set; 
  // more precisely, as an ordered set of ordered intervals.
  // The class for an ordered set of ordered intervals is
  // intervalSet.  Lower-level methods and operators take intervalSet
  // and interval arguments.  

  ////////////////////////////////////////////////////////////////////////////
  // Adjunction
  ////////////////////////////////////////////////////////////////////////////
  // We can also add an individual entity to any existing entitySet using the
  // += operator, for example we can include the entity e10 in set B:
  B += e10 ;

  //  A = ([0,9])
  //  B = ([10,10][14,100])
  //  C = ([5,15])
  cout << "A = " << A << endl ;
  cout << "B = " << B << endl ;
  cout << "C = " << C << endl ;
  cout << "F = " << F << endl ;

  /////////////////////////////////////////////////////////////////////////////
  // num_intervals
  /////////////////////////////////////////////////////////////////////////////
  // The num_intervals() method returns the number of intervals in the
  // internal representation of an entitySet as an intervalSet.
  cout << "B = " << B << endl ;
  // B.num_intervals() = 2
  cout << "B.num_intervals() = " << B.num_intervals() << endl ; 

  ////////////////////////////////////////////////////////////////////////////
  // Neither order nor duplication matters to an entitySet.
  // For example
  entitySet E = B + C ;
  //  E = B union C =  ([5,100])
  // [gives the set ([5,100]) without duplicating 14 and 15]
  cout << "E = B union C =  " << E << endl ;

  ////////////////////////////////////////////////////////////////////////////
  // Distinguished constants
  ////////////////////////////////////////////////////////////////////////////
  // EMPTY is a distinguished constant for the empty set.
  cout << "EMPTY = " << EMPTY << endl ;

  
  // UNIVERSE_MAX and UNIVERSE_MIN are distinguished constants for the
  // largest and smallest integers that may be used to label entities.
  cout << "UNIVERSE_MAX = " << Loci::UNIVERSE_MAX << endl ;
  cout << "UNIVERSE_MIN = " << Loci::UNIVERSE_MIN << endl ;

  // Another useful derived constant is ~EMPTY or not EMPTY, the '~' operator
  // provides the set complement with respect to the universal set, therefore
  // ~EMPTY is the universal set (Set of all possible entities).

  cout << "Universal Set = " << ~EMPTY << endl ;
  
  ////////////////////////////////////////////////////////////////////////////
  // Set membership
  ////////////////////////////////////////////////////////////////////////////
  // Method inSet tests whether an entity is an element of an
  // entitySet.
  // The argument to inSet is the integer label of the entity.
  cout << "A.inSet(9) = " ; 
  if ( A.inSet(9) )
      cout << "TRUE." ;
  cout << endl ;
  cout << "A.inSet(10) = " ; 
  if ( ! A.inSet(10) )
      cout << "FALSE." ;  
  cout << endl ;

  ////////////////////////////////////////////////////////////////////////////
  // Set operations
  ////////////////////////////////////////////////////////////////////////////
  // entitySet supports union, intersection, relative complement, 
  // and complement relative
  // to the set of permitted identifiers.
  //
  // For example, entitySet D becomes the union of A and B
  // That is D = ([0-9],[14-100])
  entitySet D = A + B ;
  //  D = ([0,10][14,100])
  cout << "D = " << D << endl ;
  // A & C  gives the intersection of A and C (the interval [5-9])
  cout << "A intersect C = " << (A & C) << endl ;
  // A - C gives the set difference (A take away C) (the interval [0-4])
  cout << "relative complement of A in C = " << (A - C) << endl ;
  cout << "A intersect (B union F) = " << (A & (B + F)) << endl ;  
  cout << "A union (B intersect F) = " << (A + (B & F)) << endl ;

  ////////////////////////////////////////////////////////////////////////////
  // Explicit membership list
  // We can also create an entitySet from an arbitrary list 
  // of entity identifiers.
  // For example,
  int values[] = {10,15,12,1,14,16,17} ;
  entitySet vset = create_entitySet(values,values+7) ;
  // vset = ([1,1][10,10][12,12][14,17])
  cout << "vset = " << vset << endl ;

  ////////////////////////////////////////////////////////////////////////////
  // Creating an entitySet from a vector of integers.
  // create_entitySet works with std::vector and begin() and end().
  // (See a standard C++ book for more details on using vector<> and other
  // STL containers.)
  vector<int> vec ;
  for(int i=10;i>0;--i)
    vec.push_back(i) ;
  entitySet vecset = create_entitySet(vec.begin(),vec.end()) ;
  // vecset = ([1,10])
  cout << "vecset = " <<vecset << endl ;

  ////////////////////////////////////////////////////////////////////////////
  // Iteration over an entitySet
  // We can also iterate (loop) over an entitySet in a fashion similar to
  // how we iterate over standard C++ containers.  For example, to iterate
  // over all of the entities in vset we would write a loop as follows.
  // First we create an iterator ei for entity sets.
  entitySet::const_iterator ei ;
  // Then we use the iterator to loop over a given entitySet using the
  // begin() and end() methods.  For example to loop over the entities
  // contained in vset we write:
  cout << "looping over vset:" ;
  for(ei = vset.begin(); ei != vset.end(); ++ei) 
    cout << " " << *ei ;
  cout << endl ;
  // above outputs:
  // looping over vset: 1 10 12 14 15 16 17

  ////////////////////////////////////////////////////////////////////////////
  // Min, Max, size, set membership
  // Other useful methods include Min() and Max() which can be used to find
  // the largest and smallest integer labels of entities contained in a given
  // entitySet.  For example, vset.Max() == 17 and vset.Min() == 1

  // vset.Min() == 1, vset.Max() == 17
  cout << "vset.Min() = " << vset.Min() << endl ;
  cout << "vset.Max() = " << vset.Max() << endl ;

  // Similarly, we can check the number of entities contained within an
  // entitySet by using the size() method.  For example

  // vset.size() == 7
  cout << "vset contains " << vset.size() << " entities" << endl ;

  // We can also check to see if a particular entity label is in a
  // entitySet using the inSet() method.  For example
  if(A.inSet(5))
    cout << "entity labeled 5 is in entitySet A" << endl ;

  ////////////////////////////////////////////////////////////////////////////
  // Equal, less_than, greater_than, Union (interval), 
  // Union (entitySet), Intersection (interval), 
  // Intersection (entitySet), [absolute] Complement, 
  // Print, Input
  //
  // The entitySet.Equal method tests whether the sets are equal.
  if ( A.Equal((A & C) + (A - C)) )
    cout << "A = ((A & C) + (A - C))." << endl ;
  else
    cout << "A != ((A & C) + (A - C))." << endl ;

  //  The entitySet.less_than and greater_than methods 
  //  provide a linear ordering of entitySets.
  //  The linear ordering is lexical.
  if ( A.less_than (A & C) )
    cout << "A&C <= A." << endl ;


  if ( A.less_than (A) )
    cout << "A <= A." << endl ;


  entitySet G = A - C ;
  if ( G.greater_than (A) )
    cout << "A >= A - C." << endl ;

  // Methods for union and intersection, for intervals and for
  // intervalSets.  (Remember that an entitySet is currently
  // implemented as an intervalSet.) 
  // Note that these methods modify the object to which they belong.  
  cout << "G = " << G << "." << endl ;
  G.Union (interval(13,32)) ;
  cout << "G Union [13,32] = " << G << "." << endl ;
  cout << "G = " << G << "." << endl ;
  G.Union (D) ;
  cout << "G Union D = " << G << "." << endl ;
  cout << "G = " << G << "." << endl ;
  G.Intersection (interval(13,32)) ;
  cout << "G Intersection [13,32] = " << G << "." << endl ;
  cout << "G = " << G << "." << endl ;
  G.Intersection (D) ;
  cout << "G Intersection D = " << G << "." << endl ;
  cout << "G = " << G << "." << endl ; 

  // Method for absolute complement of an entitySet.
  // Note that the complement replaces the given set.
  // G = ([14,32]).
  G.Complement() ;
  // G = ([#,13][33,#]).
  // The first occurrence of "#" stands for Loci::UNIVERSE_MIN;
  // the second occurrence of "#" stands for Loci::UNIVERSE_MAX.
  cout << "Complement of G = " << G << endl ; 
  cout << "G = " << G << endl ; 

  ////////////////////////////////////////////////////////////////////////////
  // Sequences
  //
  // Sequences provide a way of storing ordered lists of entities.
  // Usually, users don't need to worry about creating sequences
  // directly in Loci, but instead Loci creates sequences to describe
  // the order in which calculations will be carried out.  However,
  // for completeness we will give some examples here of how to create
  // sequences, and of how to iterate over sequences.

  // We can create an arbitrary sequence from a list of integers
  // in a fashion similar to the create_entitySet function.  For
  // example:
  int listvals[] = {10,15,12,1,14,16,17} ;
  
  sequence vseq = create_sequence(listvals,listvals+7) ;

  // vseq = ([10,10][15,15][12,12][1,1][14,14][16,17])
  cout << "vseq = " << vseq << endl ;

  // Note that we can also create a sequence from an entitySet, however
  // in this case the sequence will contain the contents of the entitySet
  // in increasing order.  For example
  sequence Aseq = A ;

  // Aseq = ([0-9])
  cout << "Aseq = " << Aseq << endl ;

  // Also we can append to sequences, for example:
  sequence Cseq = Aseq + vseq;

  // Cseq = ([0,10][15,15][12,12][1,1][14,14][16,17])
  cout << "Cseq = " << Cseq << endl ;

  // Similarly we can reverse the order of a sequence using the Reverse()
  // method.  For example
  Cseq.Reverse() ;

  // Cseq = ([17,16][14,14][1,1][12,12][15,15][10,0])
  cout << "reversed Cseq = " << Cseq << endl ;

  Loci::Finalize() ;
  return 0 ;
}
