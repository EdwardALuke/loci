#include <Loci.h>


///////////////////////////////////////////////////////////////////////////////
// In this file we will use a reduction rule to compute the areas of triangles.
// Here, instead of computing the areas of the triangles directly, we instead
// compute the contribution to area that each edge will make and sum these
// results into the final area.  To compute the contribution from any given
// edge, we compute the area of the triangle formed between the two edge
// nodes and the centroid of the cell.  Notice that this approach is more
// general than the triangle based computation since it is able to compute
// the area of a general polyhedral cell.

///////////////////////////////////////////////////////////////////////////////
// Reductions are defined with respect to some associative operator that has
// an identity.  For example, summation whose identity is zero, product whose
// identity is one, maximum whose identity is the smallest possible input.
// **** It is important that the operator be associative if you expect your
// **** computations to be deterministic, that is, if you expect to get the
// **** same result every time you perform a computation.

///////////////////////////////////////////////////////////////////////////////
// For the area problem we are interested in performing a reduction with
// respect to the addition operator which is associative(*) and has an
// identity of zero.
// (*) Actually floating point addition is only approximately associative.
//     We usually neglect this fact since the approximation is often good
//     enough.  However, Loci makes sure that when a parallel computation
//     is performed, all processors use the same association so that we can
//     guarantee that the result obtained on different processors will be
//     identical.
//
// Here we define the identity of the area reduction using a unit_rule.
// The unit rule will assign the identity value to the entities that have
// area.  Also, the unit rule defines what entities will have the attribute
// area, so we will use the constraint to specify what entities should
// have the property of area.

// We create a unit rule by inheriting from unit_rule defined in Loci.h
class area_unit : public unit_rule {
  // Here we create the container for area
  store<double> area ;
public:
  area_unit() {
    // We always name our containers.  
    name_store("Area",area) ;
    // In this case, the rule has no inputs, only the output which is the
    // identity value for area (i.e. zero)
    output("Area") ;
    // Here we define what entities can have the property of area.  In this
    // case we restrict the property to things we call cells.  In general,
    // the constraint is like an input that we don't use, the attributes
    // given in the constraint will not actually be computed. the attribute,
    // however, must be computable (or present in the fact_db) for
    // this rule to apply.
    constraint("cells") ;
  }
  // Here we calculate the identity, which is simply the value zero
  void calculate(Entity e) {
    area[e] = 0.0 ;
  }
  // Standard interface to the execution schedule..
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;


///////////////////////////////////////////////////////////////////////////////
// Here we add up the contributions to the left side of any given face. This
// is an application of reduction and we have to use inherit from the templated
// apply_rule class to create this rule.  The apply_rule template takes two
// arguments, 1) the type of the container for which the reduction will occur
// and 2) the operator used in the reduction. Here we are using the Summation
// operator for double, however we have 4 builtin operator templates that are
// provided for your convienience (Summation, Product, Maximum, Minimum).
// The user can also specify his own operators by creating a class and
// providing a two argument operator().  Examples will be given in other files.
//
//                                  / container\    /  operator          \    .
//                                  vvvvvvvvvvvv   vvvvvvvvvvvvvvvvvvvvvvv    .
class area_left : public apply_rule<store<double>, Loci::Summation<double> > {
  //  Here we define the containers for the inputs and outputs of our reduction
  // computation  
  const_store<vector2d<double> > centroid, pos ;
  const_MapVec<2> edge_nodes ;
  const_Map cl ;
  // Area is our output and input (we will add our contribution to what already
  // is there)
  store<double> area ;
public:
  area_left() {
    // Here we name our containers as usual
    name_store("centroid",centroid) ;
    name_store("pos",pos) ;
    name_store("edge_nodes",edge_nodes) ;
    name_store("cl",cl) ;
    name_store("Area",area) ;
    // We input the edge nodal positions and the left cell centroid
    // We also input the area (since we are adding to it)
    input("edge_nodes->pos,cl->centroid,cl->Area") ;
    // Here we output the area. Note that we are expected to write only to the
    // same area that we read from. (i.e. += type of operations)
    output("cl->Area") ;
  }
  void calculate(Entity e) {
    // Here we do the actual computations of the area of the face contribution.
    // Note that we are taking advantage of the geometric properties to ensure
    // that we should arrive at positive areas, however, it is possible that an
    // edge could provide a negative contribution for some edges if the
    // polyhedron isn't convex.
    const vector2d<double> dv1 = pos[edge_nodes[e][0]]-pos[edge_nodes[e][1]] ;
    const vector2d<double> dv2 = pos[edge_nodes[e][0]]-centroid[cl[e]] ;
    const double area_edge = 0.5*cross(dv1,dv2) ;
    // Once we have computed the area contribution from the edge, we join this
    // partial result to the running accumlation.  Since we are computing a
    // summation here we could have just as easily written
    // area[cl[e]] += area_edge
    // however, the join method is guaranteed to use the operator as
    // defined in the apply_rule definition above, and as such is a more
    // conistent approach to computing the reduction.  (We can change the
    // operator later in one place.)
    join(area[cl[e]],area_edge) ;
  }
  // Here is the standard interface to the execution schedule..
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

///////////////////////////////////////////////////////////////////////////////
// So far we have only accumulated the contribution from edges to their left
// side.  However, most edges have two sides.  Here we have an additional rule
// that accumulates to the right side.  We weren't able to do this computation
// in one rule (adding to the left and right sides) since some edges (e.g.
// boundary edges) have only one side.  However, Loci guarantees that all
// apply rules will be computed before array is used, but the order that
// these rules will be applied is not guaranteed.

// The following rule is a mirror of the above with the exception that we are
// accessing the right side cell instead of the left (and there is a sign
// change in the area computation to account for the topological difference
// of the two sides.
class area_right : public apply_rule<store<double>, Loci::Summation<double> > {
  const_store<vector2d<double> > centroid, pos ;
  const_MapVec<2> edge_nodes ;
  const_Map cr ;
  store<double> area ;
public:
  area_right() {
    name_store("centroid",centroid) ;
    name_store("pos",pos) ;
    name_store("edge_nodes",edge_nodes) ;
    name_store("cr",cr) ;
    name_store("Area",area) ;
    input("edge_nodes->pos,cr->centroid,cr->Area") ;
    output("cr->Area") ;
  }
  void calculate(Entity e) {
    const vector2d<double> dv1 = pos[edge_nodes[e][0]]-pos[edge_nodes[e][1]] ;
    const vector2d<double> dv2 = pos[edge_nodes[e][0]]-centroid[cr[e]] ;
    const double area_edge = -0.5*cross(dv1,dv2) ;
    join(area[cr[e]],area_edge) ;
    warn(area_edge<0.0) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

  
  
// Here we insert the rules into the global rule list.
register_rule<area_unit> register_area_unit ;
register_rule<area_left> register_area_left ;
register_rule<area_right> register_area_right ;
