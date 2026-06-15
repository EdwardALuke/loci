#include <Loci.h>
#include "prototypes.h"

#include <iostream>
using std::cout ;
using std::endl ;
using std::istream ;
using std::cerr ;

#include <vector>
using std::vector ;

#include <algorithm>
using std::swap ;

#include <string>
using std::string ;

//*************************************************************************/
//* kill_comments(istream &s)                                             */
//* This routine skips over any whitespace or comments from the input     */
//* stream.  Upon return the stream s points to the first character of    */
//* non-commented input.                                                  */
//*************************************************************************/

void kill_comments(std::istream &s) {
  Loci::parse::kill_white_space(s) ;
  while(s.peek() == '#') {
    while((s.get()) != '\n' && !s.eof())
      /* NULL STATEMENT */ ;
    Loci::parse::kill_white_space(s) ;
  }
}


// Read in a grid consisting of triangles from ifile, store the resulting mesh
// in the fact database, facts.
void read_triangles(istream &ifile, fact_db &facts) {
  kill_comments(ifile) ;
  int num_triangles=0, num_nodes=0 ;
  ifile >> num_triangles ;
  kill_comments(ifile) ;
  ifile >> num_nodes ;
  // Once we know how many triangles and nodes we have, then we
  // need to allocate space for these.  We can ask the fact database
  // to provide us disjoint sets containing the allocated entities.
  entitySet node_alloc = facts.get_allocation(num_nodes) ;
  entitySet triangle_alloc = facts.get_allocation(num_triangles) ;

  // The allocated entities will be in contiguous blocks.  Since the
  // triangles will refer to nodes that are implicitly numbered from zero, we
  // will have to adjust the numbering to match our allocation.  We do
  // this by adding the appropriate offset computed below.
  int node_offset = node_alloc.Min() ;

  // The triangles are represented as a map that relates each triangle to
  // its three nodes.
  MapVec<3> triangle_nodes ;
  triangle_nodes.allocate(triangle_alloc) ;

  cout << "node_alloc = " << node_alloc << endl ;
  cout << "triangle_alloc = " << triangle_alloc << endl ;

  // create a constraint that represents all cells in the problem.  Useful in
  // some edge based computations.
  constraint cells ;
  *cells = triangle_alloc ;
  
  // Positions of the nodes are simply 2-D vectors associated with every node
  // entity.
  store<vector2d<double> > pos ;
  pos.allocate(node_alloc) ;

  kill_comments(ifile) ;

  // Loop over all nodes and read in the positions.
  entitySet::const_iterator ei ;
  for(ei=node_alloc.begin();ei!=node_alloc.end();++ei)
    ifile >> pos[*ei] ;

  // Now loop over all triangles and read in the three corresponding nodes
  // Remember:  Adjust the node reference so that it refers to the entities
  // allocated for nodes.
  kill_comments(ifile) ;
  for(ei=triangle_alloc.begin();ei!=triangle_alloc.end();++ei) {
    for(int i=0;i<3;++i) {
      ifile >> triangle_nodes[*ei][i] ;
      triangle_nodes[*ei][i] += node_offset ;
      
    }
  }
  

  // For the purposes of this code, we assume that triangles are stored in
  // a counter-clockwise manner.  If they are not, reverse the order
  for(ei=triangle_alloc.begin();ei!=triangle_alloc.end();++ei) {
    vector2d<double>
      p0 = pos[triangle_nodes[*ei][0]],
      p1 = pos[triangle_nodes[*ei][1]],
      p2 = pos[triangle_nodes[*ei][2]] ;
    if(cross(p0-p1,p2-p1) > 0) {
      cerr << "triangle #" << *ei-triangle_alloc.Min() << " not following right hand rule." << endl ;
      swap(triangle_nodes[*ei][0],triangle_nodes[*ei][1]) ;
    }
  }
                             

  // Add these facts to the fact database
  facts.create_fact("triangle_nodes",triangle_nodes)  ;
  facts.create_fact("cells",cells) ;
  facts.create_fact("pos",pos) ;

  // Now we read in the boundary entities
  kill_comments(ifile) ;
  int num_boundaries ;
  ifile >> num_boundaries ;
  for(int i=0;i<num_boundaries;++i) {
    // For each boundary, we read in a boundary name, and the corresponding
    // nodes that are associated with the boundary.
    // Remember:  Adjust the node reference so that it refers to the entities
    // allocated for nodes.
    kill_comments(ifile) ;
    string boundary_name ;
    ifile >> boundary_name ;
    int num_boundary_nodes ;
    ifile >> num_boundary_nodes ;
    vector<int> bnodes(num_boundary_nodes) ;
    for(int j=0;j<num_boundary_nodes;++j) {
      ifile >> bnodes[j] ;
      bnodes[j] += node_offset ;
    }
    entitySet boundary_nodes = create_entitySet(bnodes.begin(),bnodes.end()) ;
    constraint boundary_constraint ;
    *boundary_constraint = boundary_nodes ;
    string fact_name = string("BC_") + boundary_name ;
    cout << "constraint = " << fact_name
         << " nodes = " << boundary_nodes << endl ;

    // We add a fact to the fact database which is a constraint named for the
    // identified boundary.
    facts.create_fact(fact_name,boundary_constraint) ;
  }
}
