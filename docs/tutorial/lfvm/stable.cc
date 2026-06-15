#include <Loci.h>

class stable_sum_unit : public unit_rule {
  store<double> stable_sum ;
public:
  stable_sum_unit() {
    name_store("stable_sum",stable_sum) ;
    constraint("area") ;
    output("stable_sum") ;
  }
  void calculate(Entity e) {
    stable_sum[e] = 0 ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&stable_sum_unit::calculate) ;
  }
} ;

class stable_sum_left : public apply_rule<store<double>, Loci::Summation<double> >{
  store<double> stable_sum ;
  const_store<vector2d<double> > centroid, epos ;
  const_Map cl ;
  const_store<double> k ;
  const_store<double> length ;
public:
  stable_sum_left() {
    name_store("stable_sum",stable_sum) ;
    name_store("centroid",centroid) ;
    name_store("epos",epos) ;
    name_store("cl",cl) ;
    name_store("k",k) ;
    name_store("length",length) ;
    input("cl->centroid,epos,k,length") ;
    input("cl->stable_sum") ;
    output("cl->stable_sum") ;
  }
  void calculate(Entity e) {
    const int left = cl[e] ;
    const double dx = norm(epos[e]-centroid[left]) ;
    const double stable_contrib = length[e]*k[e]/(dx*dx) ;
    join(stable_sum[left],stable_contrib) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&stable_sum_left::calculate) ;
  }
} ;
    
class stable_sum_right : public apply_rule<store<double>, Loci::Summation<double> >{
  store<double> stable_sum ;
  const_store<vector2d<double> > centroid, epos ;
  const_Map cr ;
  const_store<double> k ;
  const_store<double> length ;
public:
  stable_sum_right() {
    name_store("stable_sum",stable_sum) ;
    name_store("centroid",centroid) ;
    name_store("epos",epos) ;
    name_store("cr",cr) ;
    name_store("k",k) ;
    name_store("length",length) ;
    input("cr->centroid,epos,k,length") ;
    input("cr->stable_sum") ;
    output("cr->stable_sum") ;
  }
  void calculate(Entity e) {
    const int right = cr[e] ;
    const double dx = norm(epos[e]-centroid[right]) ;
    const double stable_contrib = length[e]*k[e]/(dx*dx) ;
    join(stable_sum[right],stable_contrib) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&stable_sum_right::calculate) ;
  }
} ;


register_rule<stable_sum_unit> register_stable_sum_unit ;
register_rule<stable_sum_left> register_stable_sum_left ;
register_rule<stable_sum_right> register_stable_sum_right ;

class stable_time_unit : public unit_rule {
  param<double> time_step ;
public:
  stable_time_unit() {
    name_store("time_step",time_step) ;
    constraint("UNIVERSE") ;
    output("time_step") ;
  }
  virtual void compute(const sequence &seq) {
    *time_step = 1e9 ;
  }
} ;

class calc_stable_time : public apply_rule<param<double>, Loci::Minimum<double> > {
  const_store<double> stable_sum ;
  const_store<double> area ;
  const_param<double> density ;
  const_param<double> heat_capacity ;
  param<double> time_step ;
public:
  calc_stable_time() {
    name_store("stable_sum",stable_sum) ;
    name_store("area",area) ;
    name_store("density",density) ;
    name_store("heat_capacity",heat_capacity) ;
    name_store("time_step",time_step) ;
    input("stable_sum,area,density,heat_capacity") ;
    input("time_step") ;
    output("time_step") ;
  }
  void calculate(Entity e) {
    const double dt = 2*(density[e]*heat_capacity[e]*area[e])/stable_sum[e];
    join(time_step[e],dt) ;
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&calc_stable_time::calculate) ;
  }
} ;

register_rule<stable_time_unit> register_stable_time_unit ;
register_rule<calc_stable_time> register_calc_stable_time ;
