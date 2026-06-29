#include<Loci.h>

/////////////////////////////////////////////////////////////////////////////
//  Here we are implementing the time integration routine to perform an
//   explicit, first order forward Euler time integration.  We achieve
//   this integration by computing the following formula
//
//  e0^{n+1} = e0^{n} + \delta t rhs^{n}
//
//  Where area*e0 is the total energy contained within the cell and
//  rhs is provided by the diffusion source terms (numerical integration of
//  the spatial terms).

/////////////////////////////////////////////////////////////////////////////
// To perform this time integration we must iterate.  We do this by computing
// first computing the initial value for e0.  Here we assume that our initial
// conditions are a single constant temperature that is given in the parameter
// T_initial.  This rule that initializes the iteration is called a build
// rule.
// Since this is a point by point assignment, we use a pointwise_rule type to
// specify the build rule.  In general, build rules will inherit from either
// the pointwise_rule or singleton_rule classes.
class ic_rule : public pointwise_rule {
  // Here are the containers that we need to initialize our computations
  const_param<double> T_initial ;
  const_param<double> heat_capacity ;
  store<double> e0 ;
public:
  ic_rule() {
    // Here we name the containers
    name_store("T_initial",T_initial) ;
    name_store("heat_capacity",heat_capacity) ;
    // Here we name the output variable.  Note that the {n=0} is added to
    // the end of the variable name.  This specifies that e0 exists within
    // iteration n, and is the first iteration (when n=0).
    name_store("e0{n=0}",e0) ;
    // We input the initial temperature and the heat capacity to compute
    // the initial total energy.
    input("T_initial,heat_capacity") ;
    // Here, our output is the initial state of the variable we are
    // integrating at timestep 0.
    output("e0{n=0}") ;
    // We are only interested in computing the time integration for entities
    // that have area, here we make sure we can assign e0 to every entity
    // that has an area (i.e. all cells)
    constraint("area") ;
  }
  void calculate(Entity e) {
    // Here the computation is simple.  Note that this is just the product
    // of two parameters.  We are accessing these containers using the []
    // operator; however, we could have just as easily used the * operator.
    e0[e] = T_initial[e]*heat_capacity[e] ;
  }
  virtual void compute(const sequence &seq) {
    // As usual for a pointwise rule, we loop over the entities for which
    // we are interested in computing.
    do_loop(seq,this) ;
  }
} ;

/////////////////////////////////////////////////////////////////////////////
//  Register the initial condition rule in the global_rule_list
register_rule<ic_rule> register_ic_rule ;

/////////////////////////////////////////////////////////////////////////////
//
// Now that we have computed the value of e0 at the first timestep, we need
// to provide a rule that tells us how to advance in time.  (A rule that
// specifies this is called an advance rule.)  Here we compute the
// value of e0 at time level n+1 using an explicit time integration method
// described above.
class advance_heat : public pointwise_rule {
  // We name the variables that we are about to use.
  const_store<double> e0n ;
  const_param<double> time_step ;
  const_store<double> rhs ;
  store<double> e0np1 ;
public:
  advance_heat() {
    // Here we name our containers.  We are using the notation {n}
    // and {n+1} to indicate values we are inputing from time level
    // {n} and using to compute values at iteration level {n+1}.  This rule
    // actually defines the iteration.
    name_store("e0{n}",e0n) ;
    name_store("time_step{n}",time_step) ;

    // Note here also, that we are using values of rhs at time level n,
    // however, the rule that computes rhs has no iteration specification.
    // Loci automatically promotes all rules that have no time specification
    // to the level of iteration where their dependencies require.  For
    // example, since rhs is computed using e0, it will be computed
    // as part of the iteration process.  However, the only thing the
    // user needs to specify about iteration is the iteration rules
    // in this file.  The advantage of leaving most rules independent
    // of the iteration specification is that the computations can be
    // reused in other contexts (for example, an implicit integrator).
    name_store("rhs{n}",rhs) ;
    name_store("e0{n+1}",e0np1) ;
    input("time_step{n},rhs{n},e0{n}") ;
    output("e0{n+1}") ;
    // We usually constrain our advance rules so that they neither create
    // nor destroy attributes from the previous iteration.
    constraint("e0{n}") ;
  }
  void calculate(Entity e) {
    // Here is the explicit Euler time integration
    e0np1[e] = e0n[e] + time_step[e]*rhs[e] ;
  }
  virtual void compute(const sequence &seq) {
    // Here is the standard loop interface
    do_loop(seq,this) ;
  }
} ;

/////////////////////////////////////////////////////////////////////////////
//
// Here we register the advance rule for performing the euler time integration
register_rule<advance_heat> register_advance_heat ;


/////////////////////////////////////////////////////////////////////////////
//
// Here we are computing the final result and terminating the iteration.
// Only the values computed in the collapse rule(s) will be seen outside
// of the iteration.  Here we are collapsing values associated with the
// node temperatures to the variable "solution".  Notice, that we don't
// have to collapse the same variable that we are advancing. We use the
// boolean parameter simulation finished to determine when the
// collapse takes place.  The actual computation of simulation_finished
// is provided by another rule that follows.
class collapse_time : public pointwise_rule {
  // Here we allocate the containers that we require.
  const_store<double> Tnode ;
  store<double> solution ;
public:
  collapse_time() {
    // Here we name the containers.  Note, this is a collapse rule.
    // We determine this by noting that it reads in values at a time
    // level {n} while outputing values at stationary time (solution).
    name_store("Tnode{n}",Tnode) ;
    name_store("solution",solution) ;
    input("Tnode{n}") ;
    output("solution") ;
    // Here is how we say that this rule is evaluated under the conditions
    // that param<bool> simulation_finished evaluates to true
    conditional("simulation_finished{n}") ;
  }
  // The calculation for this rule is trivial...  copy
  void calculate(Entity e) {
    solution[e] = Tnode[e] ;
  }

  // Obligatory compute loop that translates above local computation to
  // an aggregate one.
  virtual void compute(const sequence &seq) {
    do_loop(seq,this) ;
  }
} ;

// Here again we register the rule so it is seen in the global_rule_list
register_rule<collapse_time> register_collapse_time ;

/////////////////////////////////////////////////////////////////////////////
//
// We can have more than one collapse rule.  Here we are computing the
// time that the solution represents in real simulation space.
// Note, that all collapse rules are required to share the same collapse
// condition.

class collapse_info : public singleton_rule {
  const_param<int> n ;
  const_param<double> time_step ;
  param<double> solution_time ;
public:
  collapse_info() {
    name_store("$n{n}",n) ;
    name_store("time_step{n}",time_step) ;
    name_store("solution_time",solution_time) ;
    input("$n{n},time_step{n}") ;
    output("solution_time") ;
    conditional("simulation_finished{n}") ;
  }

  virtual void compute(const sequence &seq) {
    *solution_time = double(*n)*(*time_step) ;
  }
} ;


register_rule<collapse_info> register_collapse_info ;

/////////////////////////////////////////////////////////////////////////////
//
// Here we compute the condition that determines when iteration is complete.
// In this case we have a maximum number of iterations or a maximum simulated
// time that determines iteration termination.  This is a singleton rule
// since it accesses only parameters.
class collapse_condition : public singleton_rule {
  const_param<int> n, max_iteration ;
  const_param<double> max_time ;
  const_param<double> time_step ;
  param<bool> simulation_finished ;
public:
  collapse_condition() {
    //****** Here is some special magic:
    // We can chose to access that value of our iteration variable
    // (in this case n) by placing a $ before it in its name.
    // By this $n contains the current iteration number
    name_store("$n",n) ;
    name_store("max_iteration",max_iteration) ;
    name_store("time_step",time_step) ;
    name_store("max_time",max_time) ;
    name_store("simulation_finished",simulation_finished) ;
    input("$n,max_iteration,max_time,time_step") ; 
    output("simulation_finished") ;
  }
  // Note since this is a singleton rule, we have no do_loop or
  // calculate method.
  virtual void compute(const sequence &seq) {
    *simulation_finished = (*n >= *max_iteration) 
      ||((double(*n) * *time_step) >= *max_time);
  }
} ;

register_rule<collapse_condition> register_collapse_condition ;
