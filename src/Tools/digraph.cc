//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Tools/digraph.h>
#include <Tools/stream.h>
#include <Tools/debug.h>

#include <vector>
using std::vector ;
#include <map>
using std::map ;

namespace Loci {

  digraph::digraph() {}
  digraph::~digraph() {}

  digraph::digraphRep::digraphRep() {}

  digraph::digraphRep::~digraphRep() {}

  void digraph::digraphRep::add_edges(const vertexSet &ns, int j) 
    {
      vertexSet::const_iterator i ;
      for(i=ns.begin();i!=ns.end();++i)
        graph[*i] += j ;
      source_vertices += ns ;
      all_vertices += ns ;
      all_vertices += j ;
    }

  void digraph::digraphRep::remove_edge(int i, int j) 
  {
    graph_matrix::const_iterator mi ;
    mi=graph.find(i) ;
    if(mi==graph.end())
      return ;
    
    graph[i] -= j ;
    if(graph[i] == EMPTY)
      source_vertices -= i ;
  }
  
  void digraph::digraphRep::remove_edges(int i, const vertexSet &ns)
  {
    graph_matrix::const_iterator mi ;
    mi=graph.find(i) ;
    if(mi==graph.end())
      return ;
    graph[i] -= ns ;
    if(graph[i] == EMPTY)
      source_vertices -= i ;
  }
  
  void digraph::digraphRep::remove_vertex(int i)
  {
    if(source_vertices.inSet(i))
      source_vertices -= i ;
    // erase the record in the graph
    graph.erase(i) ;
    all_vertices -= i ;

    // remove edges point to this vertex
    graph_matrix::iterator ii ;
    for(ii=graph.begin();ii!=graph.end();++ii) {
      ii->second -= i ;
      // if there are no edges go out from a vertex,
      // take it off from the source_vertices set
      if(ii->second == EMPTY) {
        source_vertices -= ii->first ;
      }
    }
  }

  void digraph::digraphRep::remove_dangling_vertices()
  {
    vertexSet dangling_vertices = all_vertices - source_vertices
      - get_target_vertices() ;
    if(dangling_vertices != EMPTY) {
      for(vertexSet::const_iterator vi=dangling_vertices.begin();
	  vi!=dangling_vertices.end();++vi)
	graph.erase(*vi) ;
      all_vertices -= dangling_vertices ;
    }
  }
  
  void digraph::digraphRep::add_graph(const digraphRep &gr) {
    graph_matrix::const_iterator ii ;
    for(ii=gr.graph.begin();ii!=gr.graph.end();++ii)
      add_edges(ii->first,ii->second) ;
  }

  void digraph::digraphRep::subtract_graph(const digraphRep &gr) {
    // first get all dangling nodes in the graph
    vertexSet old_dangling_vertices = all_vertices - source_vertices 
      - get_target_vertices() ;
    // then remove any edges that are shared in both graph
    graph_matrix::const_iterator ii ;
    for(ii=gr.graph.begin();ii!=gr.graph.end();++ii) 
      remove_edges(ii->first,ii->second) ;
    // then we get all the dangling nodes again
    vertexSet new_dangling_vertices = all_vertices - source_vertices
      - get_target_vertices() ;
    // we compute the difference of dangling nodes, and remove them
    vertexSet remove = new_dangling_vertices - old_dangling_vertices ;

    vertexSet::const_iterator vi ;
    for(vi=remove.begin();vi!=remove.end();++vi) {
      graph.erase(*vi) ;
      all_vertices -= *vi ;
    }
  }

  void digraph::digraphRep::subgraph(const vertexSet &ns) {
    graph_matrix::iterator ii ;
    source_vertices = source_vertices & ns ;
    all_vertices = all_vertices & ns ;

    vertexSet empty_sources ;
    for(ii=graph.begin();ii!=graph.end();) 
      if(!ns.inSet(ii->first)) {
        // erase records that are no longer the graph's data
        graph.erase(ii++) ;
      }
      else {
        vertexSet isect = ii->second & ns ;
        if(ii->second != isect)
          ii->second = isect ;
        if(isect == EMPTY)
          empty_sources += ii->first ;

        // increment ii
        ++ii ;
      }
    source_vertices -= empty_sources ;
  }

  int digraph::max_vertex() const {
    vertexSet allv = get_all_vertices() ;
    return allv==EMPTY?0:allv.Max() ;
  }


  // Perform topological sort on the name_tags based on the graph described by
  // registered function dependencies.  This algorithm is described in
  // "Intoduction to Algorithms" by T Cormen, C Leiserson, and R Rivest.

  // Unnamed namespace causes Intel compiler problems
  //  namespace {
    enum vertex_color { WHITE, GRAY, BLACK } ;
    typedef map<int,vertex_color> vertex_color_map ;
    typedef map<int,sequence> priority_graph ;
    typedef digraph::vertexSet vertexSet ;
  //  }



  vertexSet dfs_visit(int vertex, const digraph &pg,
                    sequence &vertex_list,vertex_color_map &vertex_color)
    {
      vertex_color[vertex] = GRAY ;
      const vertexSet &edges = pg[vertex] ;

      vertexSet visited ;
      visited += vertex ;
      for(vertexSet::const_iterator ii=edges.begin();ii!=edges.end();++ii) {
        switch(vertex_color[*ii]) {
        case GRAY:
          // Backedge, indicates recursion
          //            back_edges.add_edge(vertex,*ii) ;
          //            back_edges_exist = true ;
          break ;
        case BLACK:
          // Already Searched so skip
          break ;
        case WHITE:
          // Not visited, continue depth first search
          visited += dfs_visit(*ii,pg,vertex_list,vertex_color) ;
          break ;
        default:
#ifdef DEBUG
          cerr << "Default case in dfs_visit should not occur" << endl ;
#endif
          break ;
        }
      }
      vertex_color[vertex] = BLACK ;
      vertex_list += vertex ;
      return visited ;
    }

  component_sort::component_sort(const digraph &dg)
    {
      if(dg.get_source_vertices() == EMPTY)
        return ;
      vertex_color_map vertex_color ;
      vertexSet all_vertices = dg.get_source_vertices() + dg.get_target_vertices() ;
      vertexSet::const_iterator ii ;
      for(ii = all_vertices.begin() ; ii != all_vertices.end() ; ++ii)
        vertex_color[*ii] = WHITE ;

      sequence order ;
      for(ii = all_vertices.begin() ; ii != all_vertices.end() ; ++ii) {
        switch(vertex_color[*ii]) {
        case WHITE:
          dfs_visit(*ii,dg,order,vertex_color) ;
          break ;
        case BLACK:
          break ;
        default:
          cerr << "invalid vertex color in switch" << endl ;
        }
      }
      vertex_color.clear() ;

      sequence rorder = order.Reverse() ;
      vector<int> omap ;  // order map
      map<int,int> imap ; // inverse
      for(sequence::const_iterator is=rorder.begin();is!=rorder.end();++is) {
        imap[*is] = omap.size() ;
        omap.push_back(*is) ;
      }

      digraph dgt ;
      const vertexSet &sources = dg.get_source_vertices() ;
      for(ii=sources.begin();ii!=sources.end();++ii) {
        const vertexSet &edges = dg[*ii] ;
        vertexSet::const_iterator jj ;
        warn(imap.find(*ii) == imap.end()) ;
        for(jj=edges.begin();jj!=edges.end();++jj) {
          warn(imap.find(*jj) == imap.end()) ;
          dgt.add_edge(imap[*jj],imap[*ii]) ;
        }
      }
      for(size_t i=0;i<omap.size();++i)
        vertex_color[i] = WHITE ;

      all_vertices = interval(0,omap.size()-1) ;
      for(ii = all_vertices.begin() ; ii != all_vertices.end() ; ++ii)
        vertex_color[*ii] = WHITE ;

      order = EMPTY ;
      vector<vertexSet> comp_vec ;
      for(ii = all_vertices.begin() ; ii != all_vertices.end() ; ++ii) {
        vertexSet n ;
        switch(vertex_color[*ii]) {
        case WHITE:
          comp_vec.push_back(dfs_visit(*ii,dgt,order,vertex_color)) ;
          break ;
        case BLACK:
          break ;
        default:
          cerr << "invalid vertex color in switch" << endl ;
        }
      }
      vector<vertexSet>::const_iterator ci ;
      for(ci=comp_vec.begin();ci!=comp_vec.end();++ci) {
        vertexSet n ;
        for(ii=ci->begin();ii!=ci->end();++ii) {
          if(*ii < int(omap.size()))
            n += omap[*ii] ;
          else 
            cerr << "confused , *ii = " << *ii << endl ;
        }
        components.push_back(n) ;
      }

      for(sequence::const_iterator is=order.begin();is!=order.end();++is) {
        vertex_list += omap[*is] ;
      }  
    }

  digraph::vertexSet visit_vertices(digraph dg,digraph::vertexSet begin) {
    
    digraph::vertexSet visit = begin ;
    digraph::vertexSet visited ;
    digraph::vertexSet::const_iterator ni ;
    
    // keep visiting vertices until no new vertices are found
    while(visit != EMPTY) {
      digraph::vertexSet newvertices ;
      // visit all the vertices that this graph leads to
      for(ni=visit.begin();ni!=visit.end();++ni)
        newvertices += dg[*ni] ;
      // update visit, but don't re-visit a vertex
      visit = newvertices - visited ;
      visited = visited + newvertices ;
    }
    return visited ;
  }
  
}
