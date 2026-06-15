#include <Loci.h>

class calcT : public pointwise_rule {
  const_store<double> e0 ;
  const_param<double> heat_capacity ;
  store<double> T ;
public:
  calcT() {
    name_store("e0",e0) ;
    name_store("heat_capacity",heat_capacity) ;
    name_store("T",T) ;
    input("e0,heat_capacity") ;
    output("T") ;
  }
  void calculate(Entity e) {
    T[e] = e0[e]/heat_capacity[e] ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&calcT::calculate) ;
  }
} ;

register_rule<calcT> register_calcT ;

class compute_k : public pointwise_rule {
  const_param<double> conductivity ;
  store<double> k ;
public:
  compute_k() {
    name_store("conductivity",conductivity) ;
    name_store("k",k) ;
    input("conductivity") ;
    output("k") ;
  }
  void calculate(Entity e) {
    k[e] = conductivity[e] ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&compute_k::calculate) ;
  }
} ;

register_rule<compute_k> register_compute_k ;

class calc_heat_flux : public pointwise_rule {
  const_store<double> k ;
  const_store<double> gradT_dot_n ;
  store<double> heat_flux ;
public:
  calc_heat_flux() {
    name_store("k",k) ;
    name_store("gradDotN(T)",gradT_dot_n) ;
    name_store("heat_flux",heat_flux) ;
    input("k,gradDotN(T)") ;
    output("heat_flux") ;
  }
  void calculate(Entity e) {
    heat_flux[e] = gradT_dot_n[e]*k[e] ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&calc_heat_flux::calculate) ;
  }
} ;

register_rule<calc_heat_flux> register_calc_heat_flux ;

