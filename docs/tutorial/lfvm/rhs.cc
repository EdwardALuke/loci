#include <Loci.h>

// Space integration
class calc_flux_sum_unit : public unit_rule {
  store<double> flux_sum ;
public:
  calc_flux_sum_unit() {
    name_store("flux_sum",flux_sum) ;
    constraint("area") ;
    output("flux_sum") ;
  }
  void calculate(Entity e) {
    flux_sum[e] = 0 ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

class calc_flux_sum_left :
  public apply_rule<store<double>, Loci::Summation<double> > {
    const_Map cl ;
    const_store<double> heat_flux ;
    const_store<double> length ;
    store<double> flux_sum ;
public:
    calc_flux_sum_left() {
      name_store("cl",cl) ;
      name_store("heat_flux",heat_flux) ;
      name_store("length",length) ;
      name_store("flux_sum",flux_sum) ;
      input("heat_flux,length,cl->flux_sum") ;
      output("cl->flux_sum") ;
      constraint("cl->flux_sum") ;
    }
    void calculate(Entity e) {
      join(flux_sum[cl[e]],heat_flux[e]*length[e]) ;
    }
    virtual void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
} ;

class calc_flux_sum_right :
  public apply_rule<store<double>, Loci::Summation<double> > {
    const_Map cr ;
    const_store<double> heat_flux ;
    const_store<double> length ;
    store<double> flux_sum ;
public:
    calc_flux_sum_right() {
      name_store("cr",cr) ;
      name_store("heat_flux",heat_flux) ;
      name_store("flux_sum",flux_sum) ;
      name_store("length",length) ;
      input("heat_flux,length,cr->flux_sum") ;
      output("cr->flux_sum") ;
      constraint("cr->flux_sum") ;
    }
    void calculate(Entity e) {
      join(flux_sum[cr[e]],-heat_flux[e]*length[e]) ;
    }
    virtual void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
} ;

register_rule<calc_flux_sum_unit> register_calc_flux_sum_unit ;
register_rule<calc_flux_sum_left> register_calc_flux_sum_left ;
register_rule<calc_flux_sum_right> register_calc_rhs_right ;


class calc_rhs : public pointwise_rule {
  const_store<double> flux_sum ;
  const_param<double> density ;
  const_store<double> area ;
  store<double> rhs ;
public:
  calc_rhs() {
    name_store("flux_sum",flux_sum) ;
    name_store("density",density) ;
    name_store("area",area) ;
    name_store("rhs",rhs) ;
    input("flux_sum,density,area") ;
    output("rhs") ;
  }
  void calculate(Entity e) {
    rhs[e] = flux_sum[e]/(area[e]*density[e]) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<calc_rhs> register_calc_rhs ;

