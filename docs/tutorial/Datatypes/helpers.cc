#include <Loci.h>

#include <iostream>
using std::cout ;
using std::endl ;

int main(int argc, char *argv[])
{

  Loci::Init(&argc,&argv) ;
  
  // Loci provides a helper class for creating Arrays as first class objects
  // (Loci containers must contain first class objects, so if you want a
  // store to contain an array, use this class.
  // The array is a template class parameterized on the type and size of the
  // array
  Array<double,5> A ;  // Analogous to double A[5] 

  // Otherwise Arrays act much like C arrays.  The one exception is that
  // unlike C arrays, if you pass these arrays to a function they will be
  // passed by value not by reference!

  for(int i=0;i<5;++i)
    A[i] = double(i) ;

  // Also Assignment works on the entire array
  Array<double,5> B = A ;

  // Other tricks are also added (mostly for compatibility with Loci reduction
  // rules to be covered later.

  B += A ;  // Does element by element addition of array A to array B
  // Operators +,-,*, and / are similarly overloaded

  // Also Array supports STL iterator style access, e.g.
  double total = 0;
  for(Array<double,5>::const_iterator ii=B.begin();ii!=B.end();++ii)
    total += *ii ;

  cout << "total = " << total << endl ;

  // 3d vectors are supported as types, templated by the type of x,y, and z
  vector3d<double> v1(0.0,1.0,2.0),v2(1.0,2.0,0.0),v3(2.0,1.0,0.0) ;

  // 3d vectors overload all of the expected algebraic operations in the usual
  // way including scalar vector multiply.  For example, the vector found by
  // "averaging" the vectors v1, v2, and v3 is computed as follows.
  vector3d<double> avg = (v1+v2+v3)/3.0 ;

  // Also 3dvectors support cross and dot products
  vector3d<double> cross_product = cross(v1,v2) ;
  // Also dot and cross products can be nested...
  double           dotcross = dot(cross(v1,v2),v3) ;
  // L2 norm is also provided
  double           areav1xv2 = norm(cross_product) ;

  // Of course interaction with scalar values acts as expected.
  // Normalize cross product
  cross_product /= areav1xv2 ;

  cout << "average vector = " << avg << endl ;
  cout << "cross_product = " << cross_product << endl ;
  cout << "areav1xv2 = " << areav1xv2 << endl ;
  cout << "dotcross = " << dotcross << endl ;
  
  // 2d vectors are supported as types, templated by the type of x and y
  vector2d<double> v4(1.0,1.0),v5(-1.0,1.0),v6(2.0,3.0) ;

  // 2d vectors also support cross and dot products similar to 3d vectors
  // with the exception that cross products in 2d space are scalars (the z
  // component of a 3d cross product).
  double twodcross = cross(v4,v5) ;
  double twoddot   = dot(v5,v6) ;

  // Of course all operators are overloaded as in the 3d vector case.
  vector2d<double> twodavg = (v4+v5+v6)/3.0 ;

  cout << "twodcross = " << twodcross << endl ;
  cout << "twoddot =   " << twoddot << endl ;
  cout << "twodavg =   " << twodavg << endl ;

  Loci::Finalize() ;
  return 0 ;
}
