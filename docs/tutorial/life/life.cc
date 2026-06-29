#include <Loci.h>

// build rule
class life_build : public pointwise_rule {
  store<int> life_0 ;
  const_store<int> initial_life ;
public:
  life_build() {
    name_store("initial_life", initial_life) ;
    name_store("life{I=0}",life_0 ) ;
    input("initial_life") ; 
    output("life{I=0}") ;
  }
  void calculate(Entity e) { 
    life_0[e] = initial_life[e] ; 
  }
  virtual void compute(const sequence &seq){
    do_loop(seq,this) ;
  }
} ;

register_rule<life_build> register_life_build ;


// advance rule
class life_step : public pointwise_rule {
  const_Map N, S, E, W ;
  const_store<int> life_n ;
  store<int> life_np1 ;
public:
  life_step() {
    name_store("N",N) ;
    name_store("S",S) ;
    name_store("E",E) ;
    name_store("W",W) ;
    name_store("life{I}",life_n) ;
    name_store("life{I+1}",life_np1) ;
    input("N->life{I},E->N->life{I},W->N->life{I}") ;
    input("S->life{I},E->S->life{I},W->S->life{I}") ;
    input("E->life{I}") ;
    input("W->life{I}") ;
    input("life{I}") ;
    output("life{I+1}") ;
  }
  inline void calculate(Entity e) {
    int sum = 
      life_n[N[e]]+life_n[N[E[e]]]+life_n[N[W[e]]]+
      life_n[S[e]]+life_n[S[E[e]]]+life_n[S[W[e]]]+
      life_n[E[e]]+
      life_n[W[e]] ;
    if(life_n[e]==0) {
      life_np1[e] = (sum == 3) ;
    } else {
      life_np1[e] = (sum == 2 || sum == 3) ;
    }
  }
  void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<life_step> register_life_step ;


// collapse rule
class collapse_rule : public pointwise_rule {
  const_store<int> life ;
  store<int> final_life ;
public:
  collapse_rule() {
    name_store("life{I}",life) ;
    name_store("final_life",final_life) ;
    input("life{I}") ;
    output("final_life") ;
    conditional("life_finished{I}") ;
  }
  void calculate(Entity e) { 
    final_life[e] = life[e];
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

register_rule<collapse_rule> register_collapse_rule ;

class collapse_condition : public singleton_rule {
  const_param<int> I, max_iteration ;
  param<bool> life_finished ;
public:
  collapse_condition() {
    name_store("$I",I) ;
    name_store("max_iteration",max_iteration) ;
    name_store("life_finished",life_finished) ;
    input("$I,max_iteration") ;
    output("life_finished") ; 
  }
  virtual void compute(const sequence &seq) {
    *life_finished = (*I >= *max_iteration) ;
  }
} ;

register_rule<collapse_condition> register_collapse_condition ;
