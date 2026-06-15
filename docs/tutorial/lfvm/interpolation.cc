#include <Loci.h>

class face_temp : public pointwise_rule {
  const_store<double> T ;
  const_Map cl,cr ;
  store<double> Tface ;
public:
  face_temp() {
    name_store("T",T) ;
    name_store("cl",cl) ;
    name_store("cr",cr) ;
    name_store("Tface",Tface) ;
    input("(cl,cr)->T") ;
    output("Tface") ;
  }
  void calculate(Entity e) {
    Tface[e] = 0.5*(T[cl[e]]+T[cr[e]]) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<face_temp> register_face_temp ;

class face_temp_boundary : public pointwise_rule {
  const_store<double> T ;
  const_Map cl ;
  store<double> Tface ;
public:
  face_temp_boundary() {
    name_store("T",T) ;
    name_store("cl",cl) ;
    name_store("Tface",Tface) ;
    input("(cl)->T") ;
    output("Tface") ;
    constraint("boundary_edges") ;
  }
  void calculate(Entity e) {
    Tface[e] =T[cl[e]] ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<face_temp_boundary> register_face_temp_boundary ;

struct nodal_temp_totals {
  double Tsum ;
  int nsum ;
  nodal_temp_totals() {}
  nodal_temp_totals(double T) : Tsum(T),nsum(1) {}
  nodal_temp_totals(double T, int n) : Tsum(T),nsum(n) {}
  double Tavg() const { return Tsum/double(nsum) ; }
    
} ;


std::ostream & operator<<(std::ostream &s, const nodal_temp_totals &v) {
  s << v.Tsum << ' ' << v.nsum ;
  return s ;
}
std::istream & operator>>(std::istream &s, nodal_temp_totals &v) {
  s >> v.Tsum >> v.nsum ;
  return s ;
}

namespace Loci {
  template<> struct data_schema_traits<nodal_temp_totals> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(nodal_temp_totals()) ;
      LOCI_INSERT_TYPE(ct,nodal_temp_totals,Tsum) ;
      LOCI_INSERT_TYPE(ct,nodal_temp_totals,nsum) ;
      return DatatypeP(ct) ;
    }
  } ;
}

struct nodal_temp_totals &operator+=(nodal_temp_totals &v1,
                                        const nodal_temp_totals &v2) {
  v1.Tsum += v2.Tsum ;
  v1.nsum += v2.nsum ;
  return v1 ;
}   

class nodal_sum_unit : public unit_rule {
  store<nodal_temp_totals> nodal_temp_sum ;
public:
  nodal_sum_unit() {
    name_store("nodal_temp_sum",nodal_temp_sum) ;
    output("nodal_temp_sum") ;
    constraint("pos") ;
  }
  void calculate(Entity e) {
    nodal_temp_sum[e] = nodal_temp_totals(0,0) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }    
} ;

class face2nodesum :
  public apply_rule<store<nodal_temp_totals>,
  Loci::Summation<nodal_temp_totals> > {
  const_MapVec<2> edge_nodes ;
  const_store<double> Tface ;
  store<nodal_temp_totals> nodal_temp_sum ;
public:
  face2nodesum() {
    name_store("edge_nodes",edge_nodes) ;
    name_store("Tface",Tface) ;
    name_store("nodal_temp_sum",nodal_temp_sum) ;
    input("Tface,edge_nodes->nodal_temp_sum") ;
    output("edge_nodes->nodal_temp_sum") ;
  }
  void calculate(Entity e) {
    const Entity nd1 = edge_nodes[e][0] ;
    const Entity nd2 = edge_nodes[e][1] ;
    join(nodal_temp_sum[nd1],nodal_temp_totals(Tface[e])) ;
    join(nodal_temp_sum[nd2],nodal_temp_totals(Tface[e])) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

class nodal_temp : public pointwise_rule {
  const_store<nodal_temp_totals> nodal_temp_sum ;
  store<double> Tnode ;
public:
  nodal_temp() {
    name_store("nodal_temp_sum",nodal_temp_sum) ;
    name_store("Tnode",Tnode) ;
    input("nodal_temp_sum") ;
    output("Tnode") ;
  }
  void calculate(Entity e) {
      Tnode[e] = nodal_temp_sum[e].Tavg() ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<nodal_sum_unit> register_nodal_sum_unit ;
register_rule<face2nodesum>   register_face2nodesum ;
register_rule<nodal_temp>     register_nodal_temp ;
