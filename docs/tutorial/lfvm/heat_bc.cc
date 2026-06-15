#include <Loci.h>

class heating_bc : public pointwise_rule {
  const_param<double> bc_heat_density ;
  store<double> heat_flux ;
public:
  heating_bc() {
    name_store("bc_heat_density",bc_heat_density) ;
    name_store("heat_flux",heat_flux) ;
    constraint("boundary_edges,edge_nodes->heat_flux_boundary") ;
    input("bc_heat_density") ;
    output("heat_flux") ;
  }
  void calculate(Entity e) {
    heat_flux[e] = bc_heat_density[e] ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

class adiabatic_bc : public pointwise_rule {
  store<double> heat_flux ;
public:
  adiabatic_bc() {
    name_store("heat_flux",heat_flux) ;
    constraint("adiabatic_bc") ;
    output("heat_flux") ;
  }
  void calculate(Entity e) {
    heat_flux[e] = 0.0 ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

class Tset_bc : public pointwise_rule {
  const_store<double> k,T ;
  const_param<double> T_initial ;
  const_store<vector2d<double> > epos,centroid ;
  const_Map cl ;
  store<double> heat_flux ;
public:
  Tset_bc() {
    name_store("k",k) ;
    name_store("T",T) ;
    name_store("T_initial",T_initial) ;
    name_store("cl",cl) ;
    name_store("epos",epos) ;
    name_store("centroid",centroid) ;
    name_store("heat_flux",heat_flux) ;
    input("k,T_initial,epos,cl->(centroid,T)") ;
    constraint("boundary_edges,edge_nodes->set_temperature_boundary") ;
    output("heat_flux") ;
  }
  void calculate(Entity e) {
    vector2d<double> dv = centroid[cl[e]]-epos[e] ;
    double gradT = (T_initial[e]-T[cl[e]])/norm(dv) ;
    heat_flux[e] = gradT*k[e] ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<heating_bc> register_heating_bc ;
register_rule<adiabatic_bc> register_adiabatic_bc ;
register_rule<Tset_bc> register_Tset_bc ;
