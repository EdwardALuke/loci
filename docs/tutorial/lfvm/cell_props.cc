// All Loci Programs include this
#include <Loci.h>

/////////////////////////////////////////////////////////////////////////
// This is a pointwise rule used to compute the centroid of a triangle.
// For triangles, the centroid is simply the average of the three nodal
// positions.
// Since this computation is performed triangle by triangle, it is a
// pointwise rule.
class triangle_centroid : public pointwise_rule {
  ///////////////////////////////////////////////////////////////////////
  // Here we define the containers that hold the input and output
  // variables.  Input variables are stored in const containers
  // to enforce their input property.  
  const_store<vector2d<double> > pos ;
  const_MapVec<3> triangle_nodes ;
  store<vector2d<double> > centroid ;
public:
  triangle_centroid() {
    ////////////////////////////////////////////////////////////////////
    // The constructor of this rule has two responsibilities:
    // 1) Give a symbolic name to the containers used in the rule
    //    (This facilitates binding these containers to values
    //     stored in the fact database durring execution scheduling.
    // 2) Provide documentation on how this rule accesses values.
    //    (This documentation is used to build a correct execution
    //     schedule, as well as communication schedule in the
    //     parallel case, so it is important that it accurately
    //     reflect the actual computations provided.

    ////////////////////////////////////////////////////////////////////
    // Here we provide the symbolic names for all containers used in
    // this computation.
    name_store("pos",pos) ;
    name_store("triangle_nodes",triangle_nodes) ;
    name_store("centroid",centroid) ;

    ////////////////////////////////////////////////////////////////////
    // Here we document the input.  We are indicating that we will
    // access values "pos" using indirection through the map
    // triangle_nodes.  Loci will guarantee that this computation
    // will only occur for entities for which both triangle_nodes
    // is defined and all three nodes pointed to by the triangle_nodes
    // map have a value called "pos".
    input("triangle_nodes->pos") ;

    ////////////////////////////////////////////////////////////////////
    // Here we are declaring that the value we produce as a result of
    // this computation is called "centroid".
    output("centroid") ;
  }

  //////////////////////////////////////////////////////////////////////
  // Here we provide a routine that will calculate the centroid of a
  // single triangle.  This routine usually is inlined out of
  // performance concerns,  although it isn't actually necessary.
  void calculate(Entity e) {
    // Note that we are accessing Entity e in the centroid container
    // but we are accessing the pos container through the contents of
    // triangle_nodes. (consistent with our rule signature of
    // centroid<-triangle_nodes->pos
    centroid[e] = (1./3.)*(pos[triangle_nodes[e][0]]+
                           pos[triangle_nodes[e][1]]+
                           pos[triangle_nodes[e][2]]) ;
  }
  //////////////////////////////////////////////////////////////////////
  // This is the actual interface to the computations that the rule
  // provides.  The execution schedule will call this function when it
  // computes centroids.  It gets a sequence of entities for which
  // a centroid computation is requested as an argument.  It uses the
  // do_loop template function, provided by Loci, to loop over the
  // sequence and call calculate from above.  Most rules in Loci
  // will define a compute method that looks identical to below.
  virtual void compute(const sequence &seq) {
    // do_loop is a template function that iterates over "seq" calling
    // this->calculate() for each entity in sequence seq.  While you
    // could iterate over seq manually with a sequence::const_iterator,
    // the do_loop template is optimized for looping over sequences.
    do_loop(seq,this) ;
  }
} ;                          

////////////////////////////////////////////////////////////////////////
// Here we create a global object that will register this rule in the
// global_rule_list before main() is called.  Therefore, no outside
// references to this class is required.
register_rule<triangle_centroid> register_triangle_centroid ;



// Here is another example of a pointwise rule.  A rule that computes
// areas of triangles:
class triangle_area : public pointwise_rule {
  // Here are the containers used in the computation
  const_store<vector2d<double> > pos ;
  const_MapVec<3> triangle_nodes ;
  store<double> area ;
public:
  triangle_area() {
    // Here we provide a symbolic (string) name for the containers
    // that Loci uses to access them.
    name_store("pos",pos) ;
    name_store("triangle_nodes",triangle_nodes) ;
    name_store("area",area) ;
    // Here we provide documentation of inputs and outputs.
    input("triangle_nodes->pos") ;
    output("area") ;
  }
  // Here is the calculate method that computes the area of the triangle
  // using the cross product of two edge vectors.
  void calculate(Entity e) {
    area[e] = 0.5*cross(pos[triangle_nodes[e][2]]-pos[triangle_nodes[e][1]],
                        pos[triangle_nodes[e][0]]-pos[triangle_nodes[e][1]]) ;
    // Here is a debugging statement (opposite of assert statement).  The
    // statement is only present if DEBUG is #define'ed.
    warn(area[e]<0.0) ;
  }
  // Here is the usual compute method.
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

// Here again, we register the rule in the global_rule_list.
register_rule<triangle_area> register_triangle_area ;
