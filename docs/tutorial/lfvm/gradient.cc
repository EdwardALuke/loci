#include <Loci.h>

/////////////////////////////////////////////////////////////////////////////
//
// Here we are providing a generic or parametric rule for computing gradients
// at faces in the mesh.   This routine assumes that the two cell centroids on
// either side of the face are aligned with the normal of the face and it
// is computing the gradient in the normal direction.  In an actual model
// we would perform a more sophisticated gradient computation.

/////////////////////////////////////////////////////////////////////////////
// Our generic rule computes the variable gradDotN parameterized by one
// argument, X.  We specify this by naming the variable gradDotN(X).
// The way that this rule will actually be invoked is by providing the
// variable that we want the gradient of as the argument.  For example
// in this heat transfer code we are interested in gradients of temperature,
// so we will be querying for gradDotN(T) in other parts of the code.
// When we query for gradDotN(T), the X's in the rule below will be
// substituted with T.  i.e. grad(X)<-(cl-cr)->(centroid,X) will become
// grad(T)<-(cl,cr)->(centroid,T).

class calc_gradXdotn : public pointwise_rule {
  // Here we create the containers like any other rule
  const_Map cl,cr ;
  const_store<vector2d<double> > centroid ;
  const_store<double> X ;
  store<double> gradDotN ;
public:
  calc_gradXdotn() {
    // here we name the containers.  However, the trick is that we name
    // the variable that is part of the parametric computations X and
    // we name the output with parenthesis X to indicate that the rule
    // is parameterized by X
    name_store("cl",cl) ;
    name_store("cr",cr) ;
    name_store("centroid",centroid) ;
    // Here is our parameteric input
    name_store("X",X) ;
    // Here is our output parameterized on X
    name_store("gradDotN(X)",gradDotN) ;
    input("(cl,cr)->(centroid,X)") ;
    // Note, we use the same name in the output as we used in the name_store
    output("gradDotN(X)") ;
  }
  void calculate(Entity e) {
    // Here the calculation is trivial, we assume the function is linear
    // between the two cells and compute the resulting gradient.
    gradDotN[e] = (X[cr[e]]-X[cl[e]])/norm(centroid[cr[e]]-centroid[cl[e]]) ;
  }
  // We perform the loop just like before
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

/////////////////////////////////////////////////////////////////////////////
//
// We don't have to do anything special to register the parametric rule.
// we just use register_rule like any other rule.
register_rule<calc_gradXdotn> register_calc_gradXdotn ;
