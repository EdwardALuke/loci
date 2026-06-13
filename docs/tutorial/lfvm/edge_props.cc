#include <Loci.h>

class edge_length : public pointwise_rule {
  const_store<vector2d<double> > pos ;
  const_MapVec<2> edge_nodes ;
  store<double> length ;
public:
  edge_length() {
    name_store("pos",pos) ;
    name_store("edge_nodes",edge_nodes) ;
    name_store("length",length) ;
    input("edge_nodes->pos") ;
    output("length") ;
  }
  void calculate(Entity e) {
    length[e] = norm(pos[edge_nodes[e][0]]-pos[edge_nodes[e][1]]) ;
  }
  void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<edge_length> register_edge_length ;

class edge_position : public pointwise_rule {
  const_store<vector2d<double> > pos ;
  const_MapVec<2> edge_nodes ;
  store<vector2d<double> > epos ;
public:
  edge_position() {
    name_store("pos",pos) ;
    name_store("edge_nodes",edge_nodes) ;
    name_store("epos",epos) ;
    input("edge_nodes->pos") ;
    output("epos") ;
  }
  void calculate(Entity e) {
    epos[e] = double(0.5)*(pos[edge_nodes[e][0]]+pos[edge_nodes[e][1]]) ;
  }
  void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<edge_position> register_edge_position ;

