#include <Loci.h>
#include "prototypes.h"

#include <vector>
using std::vector ;


// A utility data structure for holding edge information
struct edge_info {
  int nd1,nd2 ; // The two nodes defining an edge
  int cl,cr ;   // The cells on the left and right side 
  edge_info(int n1,int n2,int c1,int c2) : nd1(n1),nd2(n2),cl(c1),cr(c2) {}
} ;

void setup_edges(fact_db &facts) {

  // First, get the data structure that defines the triangle nodes
  MapVec<3> triangle_nodes ;
  triangle_nodes = facts.get_fact("triangle_nodes") ;

  // The set of triangles are the domain of the triangle_nodes map
  // while the set of nodes used are the range of this map.
  entitySet triangle_set = triangle_nodes.domain() ;
  entitySet node_set = triangle_nodes.image(triangle_set) ;

  // Find an inverse mapping (one that goes from nodes to triangles
  // The domain of this inverse map will be all nodes and the range
  // will be all triangles.
  multiMap nodes2tri ;
  Loci::inverseMap(nodes2tri,triangle_nodes,node_set,triangle_set) ;


  // We are now going to find edges.  To keep from getting duplicate edges
  // we mark each triangle as it is searched.  If an edge is formed between
  // triangles, and one of them is marked then this edge is already
  // recorded, so don't record a second time.
  entitySet triangles_searched = EMPTY ;

  // Keep boundary edges separate from internal edges.  This will be useful
  // since boundaries tend to be applied to boundary edges
  vector<edge_info> edge_data ;
  vector<edge_info> boundary_edge_data ;

  // Loop over every triangle.  For each triangle, loop over edges and search
  // all neighboring triangles.  If there is a triangle that shares this edge,
  // then record this triangle in the vector of found edges (provided that it
  // is not already in the vector).  If no triangle shares the edge, then
  // the edge is on the boundary.
  entitySet::const_iterator ei ;
  for(ei=triangle_set.begin();ei!=triangle_set.end();++ei) {
    triangles_searched += *ei ;
    // Local enumeration of the edges of the triangle.
    // ledge is the left node, redge is the right node.
    const int ledge[3] = {0,1,2}, redge[3] = {1,2,0} ;
    for(int edge=0;edge<3;++edge) {
      const int nd1 = min(triangle_nodes[*ei][ledge[edge]],
                          triangle_nodes[*ei][redge[edge]]) ;
      const int nd2 = max(triangle_nodes[*ei][ledge[edge]],
                          triangle_nodes[*ei][redge[edge]]) ;
      // Sanity check, nodes of an edge should be distinct.  This statement
      // will only produce a warning message if DEBUG is defined.
      WARN(nd1 == nd2) ;
      // Search for neighboring triangle that shares this edge.
      for(int st=0;st<nodes2tri.vec_size(nd1);++st) {
        int stri = nodes2tri[nd1][st] ;
        if(stri == *ei)
          continue ;
        for(int es=0;es<3;++es) {
          const int tnd1 = min(triangle_nodes[stri][ledge[es]],
                              triangle_nodes[stri][redge[es]]) ;
          const int tnd2 = max(triangle_nodes[stri][ledge[es]],
                              triangle_nodes[stri][redge[es]]) ;
          if(tnd1 == nd1 && tnd2 == nd2) {
            if(!triangles_searched.inSet(stri)) // If not already found
              // insert new edge into vector of internal edges
              // Note, by preserving the correct order of the triangle
              // nodes, we preserve the correct ordering of left and right
              // sides of the edge.
              edge_data.push_back(edge_info(triangle_nodes[*ei][ledge[edge]],
                                            triangle_nodes[*ei][redge[edge]],
                                            *ei,stri)) ;
            goto found_triangle ;
          }
        }
      }
      
      // If no other triangle found then this must be a boundary edge.
      boundary_edge_data.push_back(edge_info(triangle_nodes[*ei][ledge[edge]],
                                             triangle_nodes[*ei][redge[edge]],
                                             *ei,-1)) ;
    found_triangle:
      ;
    }
  }

  std::cout << "found " << edge_data.size() << " internal edges, "
	    << boundary_edge_data.size() << " on boundary" << std::endl ;

  // First allocate some space
  // For:
  // 1) ghost cells (these are cells created so that edges on boundaries
  //    have cells on both sides.  Sometimes Ghost cells are used to
  //    implement boundary conditions in finite volume schemes.
  // 2) boundary edges
  // 3) internal edges

  entitySet ghost_cells = facts.get_allocation(boundary_edge_data.size()) ;
  entitySet boundary_edges = facts.get_allocation(boundary_edge_data.size()) ;
  entitySet internal_edges = facts.get_allocation(edge_data.size()) ;
  entitySet all_edges = boundary_edges+internal_edges ;

  constraint boundary_edge_constraint ;
  *boundary_edge_constraint = boundary_edges ;
  facts.create_fact("boundary_edges",boundary_edge_constraint) ;
  
  // Create data structures for the edges.  These comprise the
  // nodes that define the edge and the left and right cells (edge_nodes,
  // cl, and cr respectively).
  MapVec<2> edge_nodes ;
  Map cl,cr ;

  edge_nodes.allocate(all_edges) ;
  cl.allocate(all_edges) ;
  cr.allocate(all_edges) ;

  // gci is an iterator to ghost_cells, we use it to assign
  // a unique ghost cell number to every boundary edge.
  entitySet::const_iterator gci = ghost_cells.begin() ;

  // Loop over boundary edges and assign datastructure values
  for(ei=boundary_edges.begin();ei!=boundary_edges.end();++ei) {
    WARN(boundary_edge_data.size() == 0) ;
    edge_nodes[*ei][0] = boundary_edge_data.back().nd1 ;
    edge_nodes[*ei][1] = boundary_edge_data.back().nd2 ;
    cl[*ei] = boundary_edge_data.back().cl ;
    cr[*ei] = *gci ;
    gci++ ;
    // After each edge is inserted into the datastructure, remove it from the
    // list
    boundary_edge_data.pop_back() ;
  }
  //Loop over internal edges and assign data structure values
  for(ei=internal_edges.begin();ei!=internal_edges.end();++ei) {
    WARN(edge_data.size() == 0) ;
    
    edge_nodes[*ei][0] = edge_data.back().nd1 ;
    edge_nodes[*ei][1] = edge_data.back().nd2 ;
    cl[*ei] = edge_data.back().cl ;
    cr[*ei] = edge_data.back().cr ;
    edge_data.pop_back() ;
  }


  facts.create_fact("edge_nodes",edge_nodes) ;
  facts.create_fact("cl",cl) ;
  facts.create_fact("cr",cr) ;

  // Sanity Check.  Make sure left sides point inside

  store<vector2d<double> > pos ;
  pos = facts.get_fact("pos") ;
  for(ei=all_edges.begin();ei!=all_edges.end();++ei) {
    const int l = cl[*ei] ;
    const int nd0 = triangle_nodes[l][0] ;
    const int nd1 = triangle_nodes[l][1] ;
    const int nd2 = triangle_nodes[l][2] ;

    const int end0 = edge_nodes[*ei][0] ;
    const int end1 = edge_nodes[*ei][1] ;
    
    const vector2d<double> t_center = (1./3.)*(pos[nd0]+pos[nd1]+pos[nd2]) ;
    const vector2d<double> e_center = 0.5*(pos[end0]+pos[end1]) ;

    double test = cross(pos[end1]-pos[end0],t_center-e_center) ;
    if(test<0) {
      std::cerr << "edge " << *ei << " points in the wrong direction" << std::endl ;
    }
  }
}
