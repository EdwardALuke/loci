//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include "Config/conf.h"
#include "comp_tools.h"
#include "loci_globs.h"

#include <mpi.h>

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <list>
using std::list ;
#include <map>
using std::map ;

#include <Tools/hash_map.h>

using std::pair ;
using std::make_pair ;

using std::string ;
using std::ostream ;
using std::ostringstream ;

#include "dist_tools.h"
#include "loci_globs.h"


//#define VERBOSE

namespace Loci {
  extern bool profile_memory_usage ;
  extern bool collect_memory_info ;

  extern double LociAppPeakMemory ;
  extern double LociAppAllocRequestBeanCounting ;
  extern double LociAppFreeRequestBeanCounting ;
  extern double LociAppPeakMemoryBeanCounting ;
  extern double LociAppLargestAlloc ;
  extern variable LociAppLargestAllocVar ;
  extern double LociAppLargestFree ;
  extern variable LociAppLargestFreeVar ;
  extern double LociAppPMTemp ;

  
   
  // Create a schedule for traversing a directed acyclic graph.  This schedule
  // may be concurrent, or many vertices of the graph may be visited at each
  // step of the schedule  If the graph contains cycles, the schedule may
  // not include all of the vertices in the graph.
  vector<digraph::vertexSet> schedule_dag(const digraph &g,
                                          digraph::vertexSet start_vertices,
                                          digraph::vertexSet only_vertices) {

    digraph gt = g.transpose() ;

    vector<digraph::vertexSet> schedule ;
    // First schedule any vertices that have no edges leading into them and have
    // not been scheduled previously (in start vertices)

    digraph::vertexSet working = g.get_source_vertices() -
      (g.get_target_vertices()+start_vertices) ;
    if(working != EMPTY)
      schedule.push_back(working) ;

    // visited vertices are all vertices that have already been scheduled
    digraph::vertexSet visited_vertices = start_vertices + working ;
    // In the beginning our working set are all scheduled vertices
    working = visited_vertices ;
    while(working != EMPTY) {
      // While we have vertices to work on, compute additional vertices that
      // can be scheduled
      digraph::vertexSet new_vertices ;
      digraph::vertexSet::const_iterator ni ;
      // loop over working set and create a list of candidate vertices
      for(ni=working.begin();ni != working.end(); ++ni)
        new_vertices += g[*ni] ;

      // If a vertex has already been scheduled it can't be scheduled again,
      // so remove visited vertices
      new_vertices = new_vertices - visited_vertices    ;
      // We only schedule vertices that are also in the only_vertices set
      working = new_vertices & only_vertices ;
      new_vertices = EMPTY ;
      // Find any vertex from this working set that has had all vertices leading
      // to it scheduled
      for(ni=working.begin();ni != working.end(); ++ni)
        if((gt[*ni] & visited_vertices) == gt[*ni])
          new_vertices += *ni ;
      working = new_vertices ;
      // and these new vertices to the schedule
      if(new_vertices != EMPTY)
        schedule.push_back(new_vertices) ;
      // update visited vertices set to include scheduled vertices
      visited_vertices += new_vertices ;
    }
    return schedule ;
  }

  /*The existential information is required to generate an execution
    schedule . This routine returns a set of entities such that the
    rule can be applied over those entities. */
  void existential_rule_analysis(rule r, fact_db &facts, sched_db &scheds) {
    FATAL(r.type() == rule::INTERNAL) ;
    entitySet sources = ~EMPTY ;
    entitySet constraints = ~EMPTY ;
    entitySet my_entities = ~EMPTY ;
    const rule_impl::info &rinfo = r.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    /*The function vmap_source_exist takes into consideration the maps
      in the body of the rule . By looping over each of the sources in
      the rule and also the constraints we make sure that the
      attribute specified by the target is implied by the satisfaction
      of the attributes in the body of the rule. */
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      sources &= vmap_source_exist(*si,facts, scheds) ;
    }

    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
      constraints &= vmap_source_exist(*si,facts, scheds) ;

    entitySet comp_sources, comp_constraints;
    if(duplicate_work) {
      comp_sources = sources;
      comp_constraints = constraints;
    }
    if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      sources &= d->my_entities ;
      constraints &= d->my_entities ;
      my_entities = d->my_entities ;
    }
    if(rinfo.constraints.begin() != rinfo.constraints.end())
      if((sources & constraints) != constraints) {
        if(MPI_processes == 1) {
          cerr << "Warning, rule " << r <<
            " cannot supply all entities of constraint" << endl ;
          cerr << "constraints = " << constraints << endl ;
          cerr << "sources & constraints = " << (sources & constraints) << endl ;
        } else {
          debugout << "Warning, rule " << r <<
            " cannot supply all entities of constraint" << endl ;
          debugout << "constraints = " << constraints << endl ;
          debugout << "sources & constraints = " << (sources & constraints) << endl ;
        }
	scheds.set_error() ;

        for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
          entitySet sources = vmap_source_exist(*si,facts, scheds) ;
          sources &= my_entities ;
          if((sources & constraints) != constraints) {
            if(MPI_processes == 1)
              cerr << "sources & constraints != constraints for input"
                   << endl
                   << sources  << " -- " << *si << endl ;
            else
              debugout << "sources & constraints != constraints for input"
                       << endl
                       << sources  << " -- " << *si << endl ;

            if(si->mapping.size() > 0) {
              entitySet working = constraints ;
              for(size_t i=0;i<si->mapping.size();++i) {
                entitySet images ;
                variableSet::const_iterator vi ;
                for(vi=si->mapping[i].begin();vi!=si->mapping[i].end();++vi)
                  images |= scheds.image(*vi,working) ;
                working = images ;
              }
              variableSet::const_iterator vi ;
              for(vi=si->var.begin();vi!=si->var.end();++vi) {
                entitySet exist = scheds.variable_existence(*vi) ;
                entitySet fails = working & ~exist ;
                if(fails != EMPTY) {
                  if(MPI_processes == 1)
                    cerr  << "expecting to find variable " << *vi << " at entities " << fails << endl << *vi << " exists at entities " << exist << endl ;
                  else
                    debugout  << "expecting to find variable " << *vi << " at entities " << fails << endl << *vi << " exists at entities " << exist << endl ;
                }
              }
            }
          }
        }
      }
    sources &= constraints ;
    //The context over which the rule is applied is given by the intersection
    // of the existential information of the sources with that of the
    //  constraints.
    entitySet context = sources & constraints ;
    entitySet comp_context;
    if(duplicate_work) {
      comp_sources &= comp_constraints;
      comp_context= comp_sources;
    }
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      entitySet targets = vmap_target_exist(*si,facts,context, scheds) ;
      entitySet comp_targets;
      if(duplicate_work)
	comp_targets = vmap_target_exist(*si, facts, comp_context, scheds);
      const variableSet &tvars = si->var ;
      variableSet::const_iterator vi ;
      for(vi=tvars.begin();vi!=tvars.end();++vi) {
	scheds.set_existential_info(*vi,r,targets) ;

	//Add information for duplication computation of this rule
	if(duplicate_work) {
	  scheds.set_my_proc_able_entities(*vi, r, targets);
	  scheds.set_proc_able_entities(*vi, r, comp_targets);

	  if(r.get_info().rule_impl->thread_rule()
	     && r.targets().begin()->get_info().name != "OUTPUT") {
	    if(r.get_info().rule_impl->get_rule_class() == rule_impl::POINTWISE) {
	      if(pointwise_duplication) {
		if(!use_duplicate_model)
		  scheds.add_policy(*vi, sched_db::ALWAYS);
		else
		  scheds.add_policy(*vi, sched_db::MODEL_BASED);
	      }
	      else
		scheds.add_policy(*vi, sched_db::NEVER);
	    }
	    else if(r.get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
	      if(reduction_duplication && !isPARAMETER(facts.get_variable(*vi))) {
		if(!use_duplicate_model)
		  scheds.add_policy(*vi, sched_db::ALWAYS);
		else
		  scheds.add_policy(*vi, sched_db::MODEL_BASED);
	      }
	      else
		scheds.add_policy(*vi, sched_db::NEVER);
	    }
	    else
	      scheds.add_policy(*vi, sched_db::NEVER);
	  }
	  else
	    scheds.add_policy(*vi, sched_db::NEVER);
	}
#ifdef VERBOSE
	debugout << "rule " << r << " generating variable " << *vi
		 << " for entities " << targets << endl << endl << endl ;
#endif
      }
    }
    /*Since the unit rules are used for the existential deduction for
      the case of reduction rules these need to be treated
      separately. The information need to be communicated at this
      stage because the unit rules initializes the entities. */
    if(facts.isDistributed()) {
      if(r.get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
        WARN(r.targets().size() != 1) ;
        variable v = *r.targets().begin() ;
        entitySet exist = scheds.get_existential_info(v, r) ;
        exist += fill_entitySet(exist,facts) ;
        scheds.set_existential_info(v,r,exist) ;
	if(duplicate_work) {
	  scheds.set_proc_able_entities(v, r,exist);
	  scheds.set_my_proc_able_entities(v, r,exist);
	}
      }
    }
  }

  //is_request_modification_allowed flag defines if the function is allowed
  //to add more requests for a variable.
  //Flag should be false if we are just interested in finding context for a rule
  //for given requests without modifying requests in the sched_db.
  //Mainly, this flag is useful with work duplication with a model
  //to precalculate the execution time to make decision of variable duplication.
  entitySet vmap_target_requests(const vmap_info &vmi, const vdefmap &tvarmap,
                                 fact_db &facts, sched_db &scheds,
				 bool is_request_modification_allowed) {
    // Here we will compute the context implied by a particular target
    // mapping
    variableSet::const_iterator vi ;
    entitySet targets ;
    // First we get the queries for all the variables that the mapping
    // is applied to.
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi) {
      // The variable should be in tvarmap, but we will report an error
      // if something fishy happens
      FATAL(tvarmap.find(*vi) == tvarmap.end()) ;
      // Get the requests for variable *vi and union it with the target set.
      targets |= tvarmap.find(*vi)->second ;
    }
    // Here we do a hack to make sure that if there are multiple
    // variable that a map applies to, then they will request their
    // union.  e.g. for a->(b,c) we make sure that b and c both have
    // the same requests.
    if(is_request_modification_allowed){
      for(vi=vmi.var.begin();vi!=vmi.var.end();++vi) {
	scheds.variable_request(*vi,targets) ;
      }
    }

    // Now we are applying the mapping that is applied to the target
    // variables. We do this by finding the preimage of each map.
    // We use the union preimage since any value that is touched
    // will need to be computed.  The union preimage is in the second
    // element of the pair that preimage returns.
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      // working is the entityset that becomes the union of all preimages
      // on this level
      entitySet working = EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!scheds.is_a_Map(*vi)) ;
        working |= scheds.preimage(*vi,targets).second ;
      }
      // Now we have evaluated this map, we move targets to this level
      targets = working ;
    }
    // When we are finished, we have followed the maps back from the targets to
    // their root.  We now have the set of entities that will be in the context
    // of the rule that will be used to satisfy this set of requests.
    return targets ;
  }

  entitySet vmap_source_requests(const vmap_info &vmi, fact_db &facts,
				 entitySet context, sched_db &scheds) {
    // this routine computes the set of entities that a source mapping will
    // imply.  It does this by following the images of the mapping.
    // The resulting entitySet contains all entities that will be accessed
    // when a loop over context is executed.
    entitySet compute = context ;
    vector<variableSet>::const_iterator mi ;
    variableSet::const_iterator vi ;
    for(mi=vmi.mapping.begin();mi!=vmi.mapping.end();++mi) {
      entitySet working ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!scheds.is_a_Map(*vi)) ;
	working |= scheds.image(*vi,compute) ;
      }
      compute = working ;
    }
    return compute ;
  }

  entitySet process_rule_requests(rule r, fact_db &facts, sched_db &scheds) {
    // Internal rules should be handling the appropriate rule requests via
    // their associated compiler.
    FATAL(r.type() == rule::INTERNAL) ;

    // First we get the target variables of this rule ;
    variableSet targets = r.targets() ;
    // We will be iterating over the target variables so we need an iterator
    variableSet::const_iterator vi ;
    // The vdefmap data structure is a map from variables to entitySets.
    // We use the tvarmap to record the requests for target variables
    // Here we are filling in the requests.
    vdefmap tvarmap ;

    // Loop over target variables and get requests from fact database
    // Here we compute the context of the rule.  This is the union of all of
    // the requests for the variables that this rule produces
    set<vmap_info>::const_iterator si ;
    entitySet context,isect = ~EMPTY ;

    entitySet filter = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      filter = d->my_entities ;
      isect = d->my_entities ;
    }
      
    for(vi=targets.begin();vi!=targets.end();++vi) {
      // This is a hack for the special case of a rule with OUTPUT
      // as a target.  In that case we will request OUTPUT for
      // all entities that exist.  So we add a request for OUTPUT
      // to the fact database

      if(vi->get_info().name == string("OUTPUT"))
	scheds.variable_request(*vi,scheds.variable_existence(*vi)&filter) ;

      // Now fill tvarmap with the requested values for variable *vi
      tvarmap[*vi] = scheds.get_variable_request(r,*vi) ;
      if(r.get_info().rule_impl->get_rule_class() == rule_impl::UNIT)
	tvarmap[*vi] += scheds.get_extra_unit_request(*vi);
      if(duplicate_work) {
	if(scheds.is_duplicate_variable(*vi))
	  tvarmap[*vi] = scheds.get_proc_able_entities(*vi, r) & tvarmap[*vi];
	else
	  tvarmap[*vi] = scheds.get_my_proc_able_entities(*vi, r) & tvarmap[*vi];
      }
    }

    const rule_impl::info &rinfo = r.get_info().desc ;
    bool mapping_in_output = false ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      // Check to see if there is mapping in the output
      if(si->mapping.size() != 0)
        mapping_in_output = true ;
      // Transform the variable requests using the mapping constructs
      // in *si
      entitySet tmp = vmap_target_requests(*si,tvarmap,facts, scheds) ;
      //The context is the union
      context |= tmp ;
      isect &= tmp ;
    }

    if(mapping_in_output) {
      entitySet sources = ~EMPTY ;
      entitySet constraints = ~EMPTY ;
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;
      /*The function vmap_source_exist takes into consideration the maps
        in the body of the rule . By looping over each of the sources in
        the rule and also the constraints we make sure that the
        attribute specified by the target is implied by the satisfaction
        of the attributes in the body of the rule. */
      for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
        sources &= vmap_source_exist(*si,facts, scheds) ;
      }
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints &= vmap_source_exist(*si,facts, scheds) ;

      sources &= constraints ;
      context &= sources ;
      isect &= sources ;
    }

    // Unit rules need to apply in the clone region as well, so
    // here we make an exception for unit rules.  (this is because
    // we will be reducing to the clone region and then communicating
    // partial results.
    if(!duplicate_work) {
      if(r.get_info().rule_impl->get_rule_class() != rule_impl::UNIT) {
	context &= filter ;
      }
    }

    // If the interstection and the union are not equal, then we are in
    // danger of not properly allocating variables for computations.  It is
    // an optimization to check this. For the distributed memory version it
    // may be useful to always do this if there is a mapping in the
    // targets of the rule.
    if(isect != context) {
      entitySet working = context ;
      vector<variableSet>::const_reverse_iterator mi ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
	for(mi=si->mapping.rbegin();mi!=si->mapping.rend();++mi) {
	  entitySet tmp ;
	  for(vi=mi->begin();vi!=mi->end();++vi)
	    tmp |= scheds.image(*vi,working) ;
	  working = tmp ;
	}
	for(vi=si->var.begin();vi!=si->var.end();++vi)
	  scheds.variable_request(*vi,working) ;
      }
    }

    // Loop over all sources for this rule and pass on the requests.
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      // First map the context through source mappings
      entitySet requests;
      requests = vmap_source_requests(*si,facts,context, scheds) ;
      entitySet var ;

      // Now we have the actual requests we are making of other rules
      // so we can tell the fact database that we are now requesting
      // these values.
#ifdef VERBOSE
      debugout << "local pruning : rule " << r << " requesting variables "
	       << si->var << " for entities " << requests << endl << endl << endl  ;
#endif
      for(vi=si->var.begin();vi!=si->var.end();++vi)
	scheds.variable_request(*vi,requests) ;

      // We also need to pass the requests on to any conditional variables
      // this rule may have.
      for(vi=rinfo.conditionals.begin();vi!=rinfo.conditionals.end();++vi)
	scheds.variable_request(*vi,context) ;
    }
#ifdef VERBOSE
    debugout << "rule " << r << " computes over " << context << endl ;
#endif
    return context ;
  }

  ////////////////////////////////////////////////////////////////////
  // function version of apply_compiler::set_var_existence
  // used inside the chomp compiler
  ////////////////////////////////////////////////////////////////////
  void existential_applyrule_analysis(rule apply, fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {

      // Compute the shadow entities produced by using this apply rules.
      // Any shadow entities that we don't own we will need to exchange
      // the partial results with other processors.
      WARN(apply.targets().size() != 1) ;
      variable reduce_var = *apply.targets().begin() ;


      const rule_impl::info &rinfo = apply.get_info().desc ;

      bool outputmap = false ;
      set<vmap_info>::const_iterator si ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
	if(si->mapping.size() != 0)
	  outputmap = true ;
      }

      // If there is no mapping in the output, then there will be no
      // shadow cast from this rule application.
      if(!outputmap)
	// If we are duplicating computations, we need to collect information.
	// Therefore, even if there is no shadow cast, we need to continue
	if(!duplicate_work)
	  return ;

      entitySet sources = ~EMPTY;
      entitySet constraints = ~EMPTY;

      for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
        sources &= vmap_source_exist_apply(*si,facts,reduce_var, scheds) ;
      }
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints &= vmap_source_exist(*si,facts, scheds) ;

      entitySet comp_sources, comp_constraints;
      if(duplicate_work) {
	comp_sources = sources;
	comp_constraints = constraints;
      }

      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      sources &= d->my_entities;
      constraints &= d->my_entities;

      sources &= constraints ;

      if(duplicate_work)
	comp_sources &= comp_constraints;

      entitySet context = sources & constraints ;

      entitySet comp_context;
      if(duplicate_work)
	comp_context = comp_sources & comp_constraints;

      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        entitySet targets = vmap_target_exist(*si,facts,context, scheds) ;
	entitySet comp_targets;
	if(duplicate_work)
	  comp_targets = vmap_target_exist(*si, facts, comp_context, scheds);

        const variableSet &tvars = si->var ;
        variableSet::const_iterator vi ;
	for(vi=tvars.begin();vi!=tvars.end();++vi) {
#ifdef VERBOSE
          debugout << "shadow is " << targets << endl ;
          debugout << "shadow not owned is "
                   << targets - d->my_entities << endl
                   << "variable is " << *vi << endl ;
#endif
	  if(outputmap)
	    scheds.variable_shadow(*vi,targets) ;

	  //Collect information regarding duplication of rule computation
	  if(duplicate_work) {
	    scheds.set_my_proc_able_entities(*vi, apply, targets);
	    scheds.set_proc_able_entities(*vi, apply, comp_targets);
	    scheds.set_existential_info(*vi, apply, EMPTY);

	    if(is_intensive_rule_output_mapping(apply, facts)) {
	      scheds.add_policy(*vi, sched_db::NEVER);
	    }
	  }
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////
  // function version of apply_compiler::process_var_requests
  // used inside the chomp compiler
  ////////////////////////////////////////////////////////////////////
  entitySet process_applyrule_requests(rule apply, rule unit_tag, bool &output_mapping, fact_db &facts, sched_db &scheds) {

#ifdef VERBOSE
    debugout << "in process_applyrule_requests" << endl ;
#endif
    vdefmap tvarmap ;

    variableSet targets = apply.targets() ;
    variableSet sources = apply.sources() ;
    FATAL(targets.size() != 1) ;
    variable tvar = *(targets.begin()) ;

    if(isPARAMETER(facts.get_variable(tvar)))
      tvarmap[tvar] = scheds.variable_existence(tvar) ;
    else
      tvarmap[tvar] = scheds.get_variable_request(unit_tag,tvar) ;

    entitySet filter = ~EMPTY;
    entitySet reduce_filter = ~EMPTY;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      filter = d->my_entities ;
      if(multilevel_duplication)
	reduce_filter = d->comp_entities;
      else
	reduce_filter = d->my_entities;
    }

    if(duplicate_work && scheds.is_duplicate_variable(tvar)) {
      //If mapping in output, we will only compute entities which can be
      //definitely computed successfully on a processor
      if(scheds.is_reduction_outputmap(tvar))
	tvarmap[tvar] &= reduce_filter;
      //If we have no mapping in output we will be able to compute more entities
      else
	tvarmap[tvar] &= (reduce_filter + scheds.get_reduce_proc_able_entities(tvar));
    }

    const rule_impl::info &rinfo = apply.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    entitySet compute ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si)
      compute |= vmap_target_requests(*si,tvarmap,facts, scheds) ;

    if(!duplicate_work || !scheds.is_duplicate_variable(tvar))
      compute &= filter;

    entitySet cnstrnts = ~EMPTY ;

    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
      cnstrnts &= vmap_source_exist(*si,facts, scheds) ;

    if(!duplicate_work || !scheds.is_duplicate_variable(tvar))
      cnstrnts &= filter;

    compute &= cnstrnts ;

    output_mapping = false ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end(); ++si) {
      variableSet::const_iterator vi ;
      entitySet comp = compute;
      vector<variableSet>::const_iterator mi ;
      for(mi=si->mapping.begin();mi!=si->mapping.end();++mi) {
        output_mapping = true ;
        entitySet working ;
        for(vi=mi->begin();vi!=mi->end();++vi) {
          FATAL(!scheds.is_a_Map(*vi)) ;
          working |= scheds.image(*vi,comp) ;
        }
        comp = working ;
      }
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        if((comp - scheds.variable_existence(*vi)) != EMPTY) {
          cerr << "ERROR: Apply rule " << apply <<  endl
               << " output mapping forces application to entities where unit does not exist." << endl ;
          cerr << "error occurs for entities " <<
            entitySet(comp-scheds.variable_existence(*vi)) << endl ;
          cerr << "error occurs when applying to variable " << *vi << endl;
          cerr << "error is not recoverable, terminating scheduling process"
               << endl ;
          scheds.set_error() ;
          Loci::Abort();
        }
        scheds.add_extra_unit_request(*vi,comp) ;
      }
    }

    entitySet srcs = ~EMPTY ;

    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si)
      srcs &= vmap_source_exist(*si,facts, scheds) ;

    if(!duplicate_work || !scheds.is_duplicate_variable(tvar))
      srcs &= filter;

    if(rinfo.constraints.begin() != rinfo.constraints.end())
      if((srcs & cnstrnts & filter) != (cnstrnts & filter)) {
        if(MPI_processes == 1) {
          cerr << "Warning, reduction rule:" << apply
               << "cannot supply all entities of constraint" << endl ;
          cerr << "constraints = " << (cnstrnts&filter) << endl ;
          entitySet sac = srcs & cnstrnts & filter ;
          cerr << "srcs & constraints = " << sac << endl ;
        } else {
          debugout << "Warning, reduction rule:" << apply
                   << "cannot supply all entities of constraint" << endl ;
          debugout << "constraints = " << (cnstrnts&filter) << endl ;
          entitySet sac = srcs & cnstrnts & filter; ;
          debugout << "srcs & constraints = " << sac << endl ;
        }
        scheds.set_error();

        for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
          entitySet sources = vmap_source_exist(*si,facts, scheds) ;
          sources &= filter;
          if((sources & cnstrnts & filter) != (cnstrnts&filter)) {
            if(MPI_processes == 1)
              cerr << "sources & constraints != constraints for input"
                   << endl
                   << (sources & filter) << " -- " << *si << endl ;
            else
              debugout << "sources & constraints != constraints for input"
                       << endl
                       << (sources & filter) << " -- " << *si << endl ;

            scheds.set_error() ;

            if(si->mapping.size() > 0) {
              entitySet working = cnstrnts & filter;
              for(size_t i=0;i<si->mapping.size();++i) {
                entitySet images ;
                variableSet::const_iterator vi ;
                for(vi=si->mapping[i].begin();vi!=si->mapping[i].end();++vi)
                  images |= scheds.image(*vi,working) ;
                working = images ;
              }
              variableSet::const_iterator vi ;
              for(vi=si->var.begin();vi!=si->var.end();++vi) {
                entitySet exist = scheds.variable_existence(*vi) ;
                entitySet fails = working & ~exist ;
                if(fails != EMPTY) {
                  if(MPI_processes == 1)
                    cerr << "expecting to find variable " << *vi << " at entities " << fails << endl << *vi << " exists at entities " << exist << endl ;
                  else
                    debugout << "expecting to find variable " << *vi << " at entities " << fails << endl << *vi << " exists at entities " << exist << endl ;
                  scheds.set_error();
                }
              }
            }
          }
        }
      }
    srcs &= cnstrnts ;

    // now trim compute to what can be computed.
    compute &= srcs ;

    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      entitySet requests;
      requests = vmap_source_requests(*si,facts,compute, scheds) ;
      variableSet::const_iterator vi ;
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        variable v = *vi ;
	if(v != tvar)
	  scheds.variable_request(v,requests) ;
	else
	  scheds.add_extra_unit_request(v, requests);

#ifdef VERBOSE
	debugout << "rule " << apply << " requesting variable "
		 << v << " for entities " << requests << endl ;
#endif
      }
    }
    
    return compute ;
  }

  // special existential analysis functions designed for
  // the blackbox rules
  void existential_blackboxrule_analysis
  (rule r, fact_db &facts, sched_db &scheds) {
    FATAL(r.type() == rule::INTERNAL) ;
    entitySet sources = ~EMPTY ;
    entitySet constraints = ~EMPTY ;
    entitySet my_entities = ~EMPTY ;
    const rule_impl::info &rinfo = r.get_info().desc ;
    set<vmap_info>::const_iterator si ;
    // get the sources
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      sources &= vmap_source_exist(*si,facts, scheds) ;
    }
    // and then constraints
    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
      constraints &= vmap_source_exist(*si,facts, scheds) ;

    entitySet all_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      sources &= d->my_entities ;
      constraints &= d->my_entities ;
      my_entities = d->my_entities ;
      all_entities = d->my_entities ;
    }
    // we only check constraints requirement if and only
    // if the constraints is not "~EMPTY"
    if(!rinfo.constraints.empty() && constraints != all_entities) {
      if((sources & constraints) != constraints) {
        if(MPI_processes == 1) {
          cerr << "Warning, rule " << r <<
            " cannot supply all entities of constraint" << endl ;
          cerr << "constraints = " << constraints << endl ;
          cerr << "sources & constraints = " << (sources & constraints) << endl ;
        } else {
          debugout << "Warning, rule " << r <<
            " cannot supply all entities of constraint" << endl ;
          debugout << "constraints = " << constraints << endl ;
          debugout << "sources & constraints = " << (sources & constraints) << endl ;
        }
	scheds.set_error() ;
      } 
      sources &= constraints ;
    }
    //The context over which the rule is applied is given by the intersection
    // of the existential information of the sources with that of the
    //  constraints.
    entitySet context = sources ;

    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      entitySet targets = vmap_target_exist(*si,facts,context, scheds) ;
      const variableSet &tvars = si->var ;
      variableSet::const_iterator vi ;
      for(vi=tvars.begin();vi!=tvars.end();++vi) {
	scheds.set_existential_info(*vi,r,targets) ;
#ifdef VERBOSE
	debugout << "rule " << r << " generating variable " << *vi
		 << " for entities " << targets << endl << endl << endl ;
#endif
      }
    }
    // there is no possible for a blackbox rule to be a unit rule
  }

  entitySet process_blackboxrule_requests
  (rule r, fact_db &facts, sched_db &scheds) {
    // Internal rules should be handling the appropriate rule requests via
    // their associated compiler.
    FATAL(r.type() == rule::INTERNAL) ;

    // First we get the target variables of this rule ;
    variableSet targets = r.targets() ;
    // We will be iterating over the target variables so we need an iterator
    variableSet::const_iterator vi ;
    // The vdefmap data structure is a map from variables to entitySets.
    // We use the tvarmap to record the requests for target variables
    // Here we are filling in the requests.
    vdefmap tvarmap ;

    // Loop over target variables and get requests from fact database
    // Here we compute the context of the rule.  This is the union of all of
    // the requests for the variables that this rule produces
    set<vmap_info>::const_iterator si ;
    entitySet context,isect = ~EMPTY ;

    entitySet filter = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      filter = d->my_entities ;
      isect = d->my_entities ;
    }

    for(vi=targets.begin();vi!=targets.end();++vi) {
      // we will request all entities exist for all
      // blackbox rule targets
      scheds.variable_request(*vi,scheds.variable_existence(*vi)) ;

      tvarmap[*vi] = scheds.get_variable_request(r,*vi) ;
    }

    const rule_impl::info &rinfo = r.get_info().desc ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      // Transform the variable requests using the mapping constructs
      // in *si
      entitySet tmp = vmap_target_requests(*si, tvarmap, facts, scheds) ;
      //The context is the union
      context |= tmp ;
      isect &= tmp ;
    }

    // Loop over all sources for this rule and pass on the requests.
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      // First map the context through source mappings
      entitySet requests =
        vmap_source_requests(*si,facts,context, scheds) ;

      for(vi=si->var.begin();vi!=si->var.end();++vi)
	scheds.variable_request(*vi,requests) ;

      // We also need to pass the requests on to any conditional variables
      // this rule may have.
      for(vi=rinfo.conditionals.begin();vi!=rinfo.conditionals.end();++vi)
	scheds.variable_request(*vi,context) ;
    }

    return context ;
  }

  /* This routine, in addition to sending the entities that are not
     owned by a particular processor,  information is stored for
     performing this communication during the execution of the
     schedule . We know the entities that a particular processor is
     supposed to send (send_entities)  . But we need to inform its neighbours
     that they are supposed to receive those entities. */
  std::list<comm_info>
  put_precomm_info(vector<pair<variable,entitySet> > send_entities,
                   fact_db &facts) {


    std::list<comm_info> plist ;
    if(send_entities.size() == 0)
      return plist ;

    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      const int sesz = send_entities.size() ;
      int **send_buffer = 0 ;
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;
      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;

        recv_buffer[0] = new int[2*d->xmit_total_size*sesz+sesz*d->xmit.size()] ;
        recv_size[0] = 2*d->xmit[0].size*sesz + sesz ;

        for(size_t i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+recv_size[i-1] ;
          recv_size[i] = 2*d->xmit[i].size*sesz+sesz ;
        }
      }

      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[2*d->copy_total_size*sesz+sesz*d->copy.size()] ;
        for(size_t i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+2*d->copy[i-1].size*sesz+sesz ;
      }

      Map l2g ;
      l2g = d->l2g.Rep() ;
      store<unsigned char> key_domain ;
      key_domain = d->key_domain.Rep() ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(size_t i=0;i<d->xmit.size();++i) {
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 2,
                  MPI_COMM_WORLD, &recv_request[i] ) ;
      }
      for(size_t i=0;i<d->copy.size();++i) {
        int j=sesz ;
        for(int k=0;k<sesz;++k) {
          entitySet temp = send_entities[k].second & d->copy[i].entities ;
          send_buffer[i][k] = temp.size() ;
	  
          for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei) {
	    send_buffer[i][j++] = key_domain[*ei] ;
            send_buffer[i][j++] = l2g[*ei] ;
	  }

          if(temp != EMPTY) {
            comm_info ci ;
            ci.v = send_entities[k].first ;
            ci.processor = d->copy[i].proc ;
            ci.send_set = temp ;
            plist.push_back(ci) ;
	  }
	}
        int send_size = j ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 2,MPI_COMM_WORLD) ;
      }

      if(d->xmit.size() > 0) {
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->xmit.size(), recv_request, status) ;
        FATAL(err != MPI_SUCCESS) ;
      }


      for(size_t i=0;i<d->xmit.size();++i) {
#ifdef DEBUG
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
#endif
        int j=sesz ;
        for(int k=0;k<sesz;++k) {
          sequence seq ;
          for(int l=0;l<recv_buffer[i][k];++l) {
	    int kd = recv_buffer[i][j++] ;
            seq += d->g2lv[kd][recv_buffer[i][j++]] ;
	  }
          if(seq != EMPTY) {
            comm_info ci ;
            ci.v = send_entities[k].first ;
            ci.processor = d->xmit[i].proc ;
            ci.recv_set = seq ;
            plist.push_back(ci) ;
	  }
        }
        WARN(j!=recieved) ;
      }


      if(d->xmit.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->copy.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] recv_request ;
      delete [] status ;

    }
    return plist ;
  }

  set<vector<variableSet> >
  get_rule_output_mappings(rule my_rule, const fact_db &facts) {
    set<vector<variableSet> > return_maps;
    set<vmap_info>::const_iterator vmsi ;
    for(vmsi = my_rule.get_info().desc.targets.begin();
	vmsi != my_rule.get_info().desc.targets.end();
	++vmsi) {
      if(vmsi->mapping.size() != 0) {
	vector<variableSet> vvs ;
	for(size_t i = 0; i < vmsi->mapping.size(); ++i) {
	  variableSet v ;
	  for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
	      vi != vmsi->mapping[i].end();
	      ++vi) {
	    v += variable(*vi,time_ident()) ;
	  }
	  vvs.push_back(v) ;
	}
	return_maps.insert(vvs) ;
      }
    }
    return return_maps;
  }

  bool is_same_mappings(const vector<variableSet> &map1, const vector<variableSet> &map2) {
    if(map1.size() != map2.size())
      return false;
    for(unsigned int i = 0; i < map1.size(); i++) {
      if(map1[i] != map2[i])
	return false;
    }

    return true;

  }

  bool is_intensive_rule_output_mapping(rule my_rule, const fact_db &facts) {
    if(!rule_has_mapping_in_output(my_rule))
      return false;
    set<vector<variableSet> > my_output_mappings = get_rule_output_mappings(my_rule, facts);
    std::set<std::vector<variableSet> >::const_iterator smi ;
    for(smi = my_output_mappings.begin(); smi != my_output_mappings.end() ; ++smi) {
      std::set<std::vector<variableSet> >::const_iterator fsmi ;
      for(fsmi = facts.intensive_output_maps.begin(); fsmi != facts.intensive_output_maps.end(); ++fsmi) {
	if(is_same_mappings(*smi, *fsmi))
	  return true;
      }
    }

    return false;
  }

  bool rule_has_mapping_in_output(rule r) {
    std::set<vmap_info>::const_iterator si;
    const rule_impl::info &rinfo = r.get_info().desc ;
    for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
      if(si->mapping.size() != 0)
        return true ;
    }
    return false ;
  }

  variableSet input_variables_with_mapping(rule r) {
    variableSet vars;
    std::set<vmap_info>::const_iterator si;
    const rule_impl::info &rinfo = r.get_info().desc ;
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      if(si->mapping.size() != 0)
        vars += si->var;
    }
    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si) {
      if(si->mapping.size() != 0)
        vars += si->var;
    }
    return vars;
  }

  variableSet input_variables(rule r) {
    variableSet vars;
    std::set<vmap_info>::const_iterator si;
    const rule_impl::info &rinfo = r.get_info().desc ;
    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
      vars += si->var;
    }
    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si) {
      vars += si->var;
    }
    return vars;
  }

  /* In the case with mapping in the output we might end up computing
     values for some of the entities in the clone region. In that case
     we need to send these values to the processor that actually owns
     them. The information as to what entities are to be sent for a
     particular variable is returned by the barrier_existential_rule_analysis 
     routine. */
  /*! vlst: input, variables that need synchronization
    scheds: input and output, first existential_info of vlst is obtained from scheds.
    after send_entitySet and fill_entitySet, the existential_info of scheds is updated.
  */
  vector<pair<variable,entitySet> >
  barrier_existential_rule_analysis(variableSet vlst,
                                    fact_db &facts, sched_db &scheds) {
    vector<pair<variable,entitySet> > send_entities ;
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    vector<entitySet> exinfo ;
    vector<ruleSet> rules ;
    vector<variable> vars ;

    vector<int> exent ;
    vector<variable> send_vars ;
    vector<rule> send_rule ;
    //exinfo is originally  obtained from scheds.
    int ent = 0 ;
    for(variableSet::const_iterator vi=vlst.begin();vi!=vlst.end();++vi) {
      variable v = *vi ;
      ruleSet r = scheds.get_existential_rules(v) ;

      vars.push_back(v) ;
      rules.push_back(r) ;
    }

    for(size_t i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;
      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
        if(rule_has_mapping_in_output(*rsi)) {
          exent.push_back(ent) ;
          send_vars.push_back(v) ;
          send_rule.push_back(*rsi) ;
        }
        exinfo.push_back(scheds.get_existential_info(v, *rsi)) ;
        ent++ ;
      }
    }

    //seinfo  is the exinfo that located in my clone region
    vector<entitySet> seinfo ;
    
    map<variable,entitySet> vmap ;
    for(size_t i=0;i<send_vars.size();++i) {
      variable v = send_vars[i] ;
      entitySet send_ents = exinfo[exent[i]] - d->my_entities ;
#ifdef VERBOSE
      debugout << "v=" << v << ", (map on output) send_ents = " << send_ents << endl ;
#endif
      seinfo.push_back(send_ents) ;
      vmap[v] += send_ents ;
    }

    for(map<variable,entitySet>::const_iterator mi = vmap.begin(); mi != vmap.end(); mi++)
      send_entities.push_back(make_pair(mi->first,mi->second));
    //after I send seinfo in my clone region, I will receive more entities
    //add these entities to my exinfo
    if(seinfo.size() != 0) {
      vector<entitySet> send_sets = send_entitySet(seinfo,facts) ;
      for(size_t i=0;i<seinfo.size();++i) {
        exinfo[exent[i]] += send_sets[i] ;
        exinfo[exent[i]] &= d->my_entities ;
      }
    }
    int j = 0 ;
#ifdef VERBOSE
    for(size_t i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;
      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
        debugout << "v=" << v << ",rule ="<<*rsi
		 <<"   exinfo="<<exinfo[j++] << endl ;
      }
    }
#endif
    /*!send the exinfo in my xmit region,
      I receive more entities that fill in my clone region
      add these entities to exinfo
      finally set_existential_info of scheds with exinfo*/ 
    vector<entitySet> fill_sets = fill_entitySet(exinfo,facts) ;
    j=0;
    for(size_t i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;

      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
	exinfo[j] += fill_sets[j] ;
#ifdef VERBOSE
	debugout << "rule " << *rsi << ", fill_sets=" << fill_sets[j] << endl ;
#endif
        variableSet tvars = rsi->targets() ;
        variable rv = v ;
	// Find corresponding name from rule targets
	for(variableSet::const_iterator vi=tvars.begin();vi!=tvars.end();++vi) {
	  
	  if(vi->get_info().name == v.get_info().name &&
	     vi->get_info().namespac == v.get_info().namespac &&
	     vi->get_info().v_ids == v.get_info().v_ids) {
	    rv = *vi ;
	  }
	}

#ifdef VERBOSE
	debugout << "set_existential_info(" << rv << ",rule=" << *rsi <<
	  exinfo[j] << endl ;
#endif
	scheds.set_existential_info(rv,*rsi,exinfo[j]) ;
	++j ;
      }
    }
    return send_entities ;
  }


  /*In this routine we fill in the communication data structure needed
    for filling in the clone region . From the "copy" data structure we know
    what all we have to receive from the neighbouring processors(clone
    region information) . But we need to inform the neighbouring
    processors to send us those information also. The
    clist data structure is set up so that it stores what information
    need to be send or received from a particular processor so that
    the clone region is filled up . */

  entitySet send_requests(const entitySet& e, variable v, fact_db &facts,
                          list<comm_info> &clist) {
    entitySet re ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      int **send_buffer = 0 ;
      int **recv_buffer = 0;
      int *recv_size = 0 ;

      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;

        recv_buffer[0] = new int[2*d->xmit_total_size] ;
        recv_size[0] = 2*d->xmit[0].size ;

        for(size_t i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+2*d->xmit[i-1].size ;
          recv_size[i] = 2*d->xmit[i].size ;
        }
      }

      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[2*d->copy_total_size] ;
        for(size_t i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+2*d->copy[i-1].size ;
      }
      Map l2g ;
      l2g = d->l2g.Rep() ;
      store<unsigned char> key_domain ;
      key_domain = d->key_domain.Rep() ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(size_t i=0;i<d->xmit.size();++i) {
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 3,
                  MPI_COMM_WORLD, &recv_request[i] ) ;
      }

      for(size_t i=0;i<d->copy.size();++i) {
        entitySet temp = e & d->copy[i].entities ;

        comm_info ci ;
        ci.v = v ;
        ci.processor = d->copy[i].proc ;
        ci.recv_set = temp ;
        clist.push_back(ci) ;

        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei) {
	  send_buffer[i][j++] = key_domain[*ei] ;
          send_buffer[i][j++] = l2g[*ei] ;
	}
        int send_size = 2*temp.size() ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 3,MPI_COMM_WORLD) ;
      }

      if(d->xmit.size() > 0) {
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->xmit.size(), recv_request, status) ;
        FATAL(err != MPI_SUCCESS) ;
      }

      for(size_t i=0;i<d->xmit.size();++i) {
        int recieved ;
        MPI_Get_count(&status[i], MPI_INT, &recieved) ;
        entitySet temp ;
        for(int j=0;j<recieved;++j) {
	  int kd = recv_buffer[i][j++] ;
          temp += d->g2lv[kd][recv_buffer[i][j]] ;
	}
        re += temp ;
        comm_info ci ;
        ci.v = v ;
        ci.processor = d->xmit[i].proc ;
        ci.send_set = temp ;
        clist.push_back(ci) ;
      }

      if(d->xmit.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->copy.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] recv_request ;
      delete [] status ;

    }
    return re ;
  }


  list<comm_info>
  barrier_process_rule_requests(variableSet vars, fact_db &facts, sched_db &scheds) {
    list<comm_info> clist ;
    entitySet reduce_filter = ~EMPTY ;
    fact_db::distribute_infoP d;
    if(facts.isDistributed())
      d = facts.get_distribute_info() ;

    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      variable v = *vi ;
      entitySet requests = scheds.get_variable_requests(v) ;
      if(duplicate_work) {
	//Find information for reduction variables
	ruleSet r = scheds.get_existential_rules(v);
	bool reduction = false;

	for(ruleSet::const_iterator ri = r.begin();
	    ri != r.end(); ri++)
	  if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT)
	    reduction = true;

	//Minimize apropriate requests
	if(scheds.is_duplicate_variable(v)) {
	  if(!reduction) {
	    for(ruleSet::const_iterator ri = r.begin();
		ri != r.end(); ri++) {
	      requests -= scheds.get_proc_able_entities(v, *ri);
	    }
	  }
	  else {
	    //We do not need to request entities to the owner if
	    //those entities can be definitely computed on this processor
	    requests -= scheds.get_reduce_proc_able_entities(v);
	  }
	}
      }

      entitySet recv_requests = send_requests(requests, v, facts, clist ) ;
      requests += recv_requests;

      //Since entities are guranteed to being computed on owner processor,
      //no need to request them on the other processors
      if(duplicate_work) {
	if(!scheds.is_duplicate_variable(v))
	  requests += fill_entitySet(requests, facts) ;
      }
      else {
	// check to see if there is mapping in the output rules ;
	ruleSet r = scheds.get_existential_rules(v);
	bool map_output = false;
	for(ruleSet::const_iterator ri = r.begin();
	    ri != r.end(); ri++)
	  if(rule_has_mapping_in_output(*ri)) 
	    map_output = true ;

	if(map_output) {
	  // If mapping in output send requests from other processors
	  requests += fill_entitySet(requests, facts) ;
	}
      }

      scheds.variable_request(v,requests) ;
    }
    return clist ;
  }

  int execute_comm2::tag_base = 1500 ;
  
  execute_comm2::execute_comm2(list<comm_info>& plist, fact_db &facts) {
    HASH_MAP(int,vector<send_unit>) send_data ;
    HASH_MAP(int,vector<recv_unit>) recv_data ;

    list<comm_info>::const_iterator cli ;
    intervalSet send_processes, recv_processes ;

    for(cli=plist.begin();cli!=plist.end();++cli) {
      int proc = cli->processor ;
      variable v = cli->v ;
      if(cli->send_set.size() > 0) {
        send_unit su ;
        su.v = v ;
        su.rep = facts.get_variable(v) ;
        su.send = cli->send_set ;

        send_data[proc].push_back(su) ;
        send_processes += proc ;
      }
      if(cli->recv_set.size() > 0) {
        recv_unit ru ;
        ru.v = v ;
        ru.rep = facts.get_variable(v) ;
        ru.recv = cli->recv_set ;
        
        recv_data[proc].push_back(ru) ;
        recv_processes += proc ;
      }
    }

    for(intervalSet::const_iterator ii=send_processes.begin();
        ii!=send_processes.end();++ii) {
      send_proc sp ;

      sp.proc = *ii ;
      sp.send_size = 0 ;
      sp.max_send_size = 0 ;
      sp.buf = 0 ;
      sp.units = send_data[*ii] ;

      send_info.push_back(sp) ;
    }
    for(intervalSet::const_iterator ii=recv_processes.begin();
        ii!=recv_processes.end();++ii) {
      recv_proc rp ;

      rp.proc = *ii ;
      // initial recv size is set to sizeof(int) because we
      // at least need to recv the size of the messeage.
      rp.recv_size = sizeof(int) ;
      rp.buf = 0 ;
      rp.units = recv_data[*ii] ;

      recv_info.push_back(rp) ;
    }
    tag1 = tag_base ;
    tag2 = tag_base+1 ;
    tag_base += 2 ;
  }

  namespace {
    // these are the shared global buffer for
    // send/recv in all execute_comm2 objects
    unsigned char* execute_comm2_send_buf = 0 ;
    int execute_comm2_send_buf_size = 0 ;
    unsigned char* execute_comm2_recv_buf = 0 ;
    int execute_comm2_recv_buf_size = 0 ;
    // small object used in execute_comm2 communication protocol
    struct recomm {
      Entity proc ;
      int idx ;
      recomm(Entity p,int i):proc(p),idx(i) {}
    } ;
  }

  void execute_comm2::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ; s.start() ;

    vector<MPI_Request> send_req(send_info.size()) ;
    vector<MPI_Request> recv_req(recv_info.size()) ;
    vector<MPI_Status> recv_status(recv_info.size()) ;
    // first we'll set up the recv size/buffer and post the request
    if(!recv_info.empty()) {
      int total_recv_size = 0 ;
      for(size_t i=0;i<recv_info.size();++i)
        total_recv_size += recv_info[i].recv_size ;

      // reallocate global buffer if its size is smaller
      if(execute_comm2_recv_buf_size < total_recv_size) {
        if(execute_comm2_recv_buf)
          delete[] execute_comm2_recv_buf ;
        execute_comm2_recv_buf = new unsigned char[total_recv_size] ;
        execute_comm2_recv_buf_size = total_recv_size ;
      }
      recv_info[0].buf = execute_comm2_recv_buf ;
      for(size_t i=1;i<recv_info.size();++i)
        recv_info[i].buf = recv_info[i-1].buf + recv_info[i-1].recv_size ;

      // post the actual recv requests
      for(size_t i=0;i<recv_info.size();++i)
        MPI_Irecv(recv_info[i].buf, recv_info[i].recv_size, MPI_PACKED,
                  recv_info[i].proc, tag1, MPI_COMM_WORLD, &recv_req[i]) ;
    }

    // now we need to setup the send size
    // we use a protocol to avoid sending two messages in most
    // cases. we will first compute the send size, if it is greater
    // than any previous size, we will send the size to the dest
    // process followed by the real message.
    vector<recomm> resend, rerecv ;
    vector<bool> resend_flag(send_info.size(),false) ;
    vector<bool> rerecv_flag(recv_info.size(),false) ;
    
    if(!send_info.empty()) {
      int total_send_size = 0 ;
      // compute the send size for each dest process
      for(size_t i=0;i<send_info.size();++i) {
        int& send_size = send_info[i].send_size ;
        send_size = 0 ;
        vector<send_unit>& units = send_info[i].units ;
        for(size_t k=0;k<units.size();++k) {
          send_unit& su = units[k] ;
          send_size += su.rep->pack_size(su.send) ;
        }
        // current size is larger than all previous ones
        // we need to inform the receiving process the size
        // we first send a message whose size is an integer
        // to that particular process
        if(send_size > send_info[i].max_send_size ||
        // this condition is needed because if the real msg size
        // is equal to the integer size, the receiving process will
        // think that it needs to re-receive the msg.
           send_size == sizeof(int)) {
          if(send_size > send_info[i].max_send_size)
            send_info[i].max_send_size = send_size ;
          send_size = sizeof(int) ;
          resend.push_back(recomm(send_info[i].proc,i)) ;
          resend_flag[i] = true ;
        }
        total_send_size += send_size ;
      }
      // check to see if the global send buffer needs to be reallocated
      if(execute_comm2_send_buf_size < total_send_size) {
        if(execute_comm2_send_buf)
          delete[] execute_comm2_send_buf ;
        execute_comm2_send_buf = new unsigned char[total_send_size] ;
        execute_comm2_send_buf_size = total_send_size ;
      }
      // do the actual pack/send
      send_info[0].buf = execute_comm2_send_buf ;
      for(size_t i=1;i<send_info.size();++i)
        send_info[i].buf = send_info[i-1].buf + send_info[i-1].send_size ;

      for(size_t i=0;i<send_info.size();++i) {
        int pack_offset = 0 ;
        vector<send_unit>& units = send_info[i].units ;
        if(resend_flag[i])
          // pack the message size
          MPI_Pack(&send_info[i].max_send_size,1,MPI_INT,send_info[i].buf,
                   send_info[i].send_size,&pack_offset,MPI_COMM_WORLD) ;
        else
          for(size_t k=0;k<units.size();++k) { // pack the actual message
            send_unit& su = units[k] ;
            su.rep->pack(send_info[i].buf,pack_offset,
                         send_info[i].send_size,su.send) ;
          }
      }
      for(size_t i=0;i<send_info.size();++i)
        MPI_Isend(send_info[i].buf,send_info[i].send_size,MPI_PACKED,
                  send_info[i].proc,tag1,MPI_COMM_WORLD,&send_req[i]) ;
    }
    // wait for all messages to complete
    if(!send_info.empty())
      MPI_Waitall(send_info.size(), &send_req[0], MPI_STATUSES_IGNORE) ;
    if(!recv_info.empty())
      MPI_Waitall(recv_info.size(), &recv_req[0], &recv_status[0]) ;

    // extract the received message
    if(!recv_info.empty()) {
      for(size_t i=0;i<recv_info.size();++i) {
        int rs ;
        MPI_Get_count(&recv_status[i], MPI_BYTE, &rs) ;
        if(rs == sizeof(int)) {
          rerecv_flag[i] = true ;
          rerecv.push_back(recomm(recv_info[i].proc,i)) ;
        }
      }
      for(size_t i=0;i<recv_info.size();++i) {
        int unpack_offset = 0 ;
        vector<recv_unit>& units = recv_info[i].units ;
        if(rerecv_flag[i]) {
          // update the recv size if necessary
          int rs ;
          MPI_Unpack(recv_info[i].buf,recv_info[i].recv_size,
                     &unpack_offset,&rs,1,MPI_INT,MPI_COMM_WORLD) ;
          if(rs > recv_info[i].recv_size)
            recv_info[i].recv_size = rs ;
        } else
          for(size_t k=0;k<units.size();++k) {
            recv_unit& ru = units[k] ;
            ru.rep->unpack(recv_info[i].buf,unpack_offset,
                           recv_info[i].recv_size,ru.recv) ;
          }
      }
    }

    // post all the rerecv request
    send_req.resize(resend.size()) ;
    recv_req.resize(rerecv.size()) ;
    recv_status.resize(rerecv.size()) ;

    for(size_t i=0;i<rerecv.size();++i) {
      int proc = rerecv[i].proc ;
      int idx = rerecv[i].idx ;
      recv_info[idx].buf = new unsigned char[recv_info[idx].recv_size] ;
      MPI_Irecv(recv_info[idx].buf,recv_info[idx].recv_size,
                MPI_PACKED,proc,tag2,MPI_COMM_WORLD,&recv_req[i]) ;
    }

    // do the repack/resend
    for(size_t i=0;i<resend.size();++i) {
      int pack_offset = 0 ;
      int idx = resend[i].idx ;
      send_info[idx].send_size = send_info[idx].max_send_size ;
      send_info[idx].buf = new unsigned char[send_info[idx].send_size] ;
      for(size_t k=0;k<send_info[idx].units.size();++k) {
        send_unit& su = send_info[idx].units[k] ;
        su.rep->pack(send_info[idx].buf,pack_offset,
                     send_info[idx].send_size,su.send) ;
      } 
    }
    for(size_t i=0;i<resend.size();++i) {
      int proc = resend[i].proc ;
      int idx = resend[i].idx ;
      MPI_Isend(send_info[idx].buf,send_info[idx].send_size,
                MPI_PACKED,proc,tag2,MPI_COMM_WORLD,&send_req[i]) ;
    }

    // wait for all send/recv to complete
    if(!resend.empty())
      MPI_Waitall(resend.size(), &send_req[0], MPI_STATUSES_IGNORE) ;
    if(!rerecv.empty())
      MPI_Waitall(rerecv.size(), &recv_req[0], &recv_status[0]) ;
    // delete all resend buffer
    for(size_t i=0;i<resend.size();++i)
      delete[] send_info[resend[i].idx].buf ;

    // extract the real message
    for(size_t i=0;i<rerecv.size();++i) {
      int idx = rerecv[i].idx ;
      vector<recv_unit>& units = recv_info[idx].units ;
      int unpack_offset = 0 ;
      for(size_t k=0;k<units.size();++k) {
        recv_unit& ru = units[k] ;
        ru.rep->unpack(recv_info[idx].buf,unpack_offset,
                       recv_info[idx].recv_size,ru.recv) ;
      }
      delete[] recv_info[idx].buf ;
    }

    double tt = s.stop() ;
    timer.addTime(tt,1) ;
  }

  void
  execute_comm2::Print(ostream& s) const {
    int sz = 0 ;
    if(send_info.size()+recv_info.size() > 0) {
      printIndent(s) ;
      s << "communication block {" << endl ;
      if(send_info.size() > 0) {
        printIndent(s) ;
        s << "Send:" << endl ;
        printIndent(s) ;
        for(size_t i=0;i<send_info.size();++i) {
          for(size_t j=0;j<send_info[i].units.size();++j) {
            s << send_info[i].units[j].v << ' ' ;
	    sz += (send_info[i].units[j].send).size() ;
	  }
	  s << " to " << send_info[i].proc << endl ;
          printIndent(s) ;
        }
	s << " Total entities sent = " << sz << endl ;
      }
      if(recv_info.size() > 0) {
        printIndent(s) ;
        s << "Recv:" << endl ;
        printIndent(s) ;
        for(size_t i=0;i<recv_info.size();++i) {
          for(size_t j=0;j<recv_info[i].units.size();++j)
            s << recv_info[i].units[j].v << ' '  ;
          s << " from " << recv_info[i].proc << endl ;
          printIndent(s) ;
        }
      }
      s << "}" << endl ;
    }
  }

  void
  execute_comm2::dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "comm: " ;

    variableSet vars  ;
    for(size_t i=0;i<send_info.size();++i)
      for(size_t j=0;j<send_info[i].units.size();++j) 
        vars += send_info[i].units[j].v ;

    for(size_t i=0;i<recv_info.size();++i) 
      for(size_t j=0;j<recv_info[i].units.size();++j)
        vars += recv_info[i].units[j].v ;

    oss << vars ;

    data_collector.accumulateTime(timer,EXEC_COMMUNICATION,oss.str()) ;
  }

  execute_comm::execute_comm(list<comm_info> &plist, fact_db &facts) {
    HASH_MAP(int,vector<send_var_info>) send_data ;
    HASH_MAP(int,vector<recv_var_info>) recv_data ;
    list<comm_info>::const_iterator cli ;
    intervalSet send_procs, recv_procs ;
    for(cli=plist.begin();cli!=plist.end();++cli) {
      variable v = cli->v ;
      if(cli->send_set.size() > 0) {
        int send_proc = cli->processor ;
        send_procs += send_proc ;
        entitySet send_set = cli->send_set ;
        send_data[send_proc].push_back(send_var_info(v,send_set)) ;
      }
      if(cli->recv_set.size() > 0) {
        int recv_proc = cli->processor ;
        sequence recv_seq = cli->recv_set ;
        recv_procs += recv_proc ;
        recv_data[recv_proc].push_back(recv_var_info(v,recv_seq)) ;
      }
    }

    for(intervalSet::const_iterator ii=send_procs.begin();
        ii!=send_procs.end();
        ++ii) {
      send_info.push_back(make_pair(*ii,send_data[*ii])) ;
      send_vars.push_back(std::vector<storeRepP>()) ;
      for(size_t i=0;i<send_data[*ii].size();++i)
        send_vars.back().push_back(facts.get_variable(send_data[*ii][i].v)) ;
    }
    for(intervalSet::const_iterator ii=recv_procs.begin();
        ii!=recv_procs.end();
        ++ii) {
      recv_info.push_back(make_pair(*ii,recv_data[*ii])) ;
      recv_vars.push_back(std::vector<storeRepP>()) ;
      for(size_t i=0;i<recv_data[*ii].size();++i)
        recv_vars.back().push_back(facts.get_variable(recv_data[*ii][i].v)) ;
    }

    /* This part sets up the memory allocation needed for the execute routine.
       Instead of allocating and deallocating memory each time for
       receiving a message from a processor we allocate a fixed
       message size for the receive buffer. Initially the maximum
       receive size and the maximum send size is set to be the size of
       an integer. This approach also reduces the cost incurred in
       sending the sizes in advance before the actual message is
       sent.
    */

    int nsend = send_info.size() ;
    int nrecv = recv_info.size() ;
    r_size = new int[nrecv] ;
    recv_sizes = new int[nrecv] ;
    maxr_size = new int[nrecv] ;
    maxs_size = new int[nsend] ;
    s_size = new int[nsend] ;
    for(int i = 0; i < nrecv; ++i) {
      r_size[i] = 0 ;
      maxr_size[i] = sizeof(int) ;
    }
    for(int i = 0; i < nsend; ++i) {
      maxs_size[i] = sizeof(int) ;
      s_size[i] = 0 ;
    }

    recv_ptr = new unsigned char*[max(nrecv,1)] ;
    send_ptr = new unsigned char*[max(nsend,1)] ;
    request =  new MPI_Request[nrecv] ;
    status =  new MPI_Status[nrecv] ;
  }

  execute_comm::~execute_comm() {
    delete [] maxr_size ;
    delete [] maxs_size ;
    delete [] r_size ;
    delete [] recv_sizes ;
    delete [] s_size ;
    delete [] recv_ptr ;
    delete [] send_ptr ;
    delete [] request ;
    delete [] status ;

  }

  static unsigned char *recv_ptr_buf = 0;
  static int recv_ptr_buf_size = 0;
  static unsigned char *send_ptr_buf = 0 ;
  static int send_ptr_buf_size = 0 ;
  void execute_comm::execute(fact_db  &facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;
    const int nrecv = recv_info.size() ;
    int resend_size = 0, rerecv_size = 0 ;
    std::vector<int> send_index ;
    std::vector<int> recv_index ;
    int total_size = 0 ;
    MPI_Request *re_request = 0 ;
    MPI_Status *re_status = 0 ;
    for(int i=0;i<nrecv;++i) {
      r_size[i] = maxr_size[i] ;
      total_size += maxr_size[i] ;
    }
    /*
      #ifdef DEBUG
      entitySet rem = entitySet((recv_info[i].second[j].seq)) - sp->domain() ;
      if(rem != EMPTY)
      debugout << "variable " << recv_info[i].second[j].v << " not allocated, but recving entities " << rem << endl ;
      #endif
    */
    if(recv_ptr_buf_size < total_size) {
      if(recv_ptr_buf)
        delete[] recv_ptr_buf ;
      recv_ptr_buf = new unsigned char[total_size] ;
      recv_ptr_buf_size = total_size ;
    }
    recv_ptr[0] = recv_ptr_buf ;
    for(int i=1;i<nrecv;++i)
      recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;

    for(int i=0;i<nrecv;++i) {
      int proc = recv_info[i].first ;
      MPI_Irecv(recv_ptr[i], r_size[i], MPI_PACKED, proc, 1,
                MPI_COMM_WORLD, &request[i]) ;
    }

    /*First we find out the size of the message we are trying to
      receive using the pack_size method associated with that
      container. For static containers pack_size returns the correct
      size for the messages to be received. But for containers like
      multiStore and storeVec whose sizes gets set at run time
      pack_size returns a erroneous value. To avoid allocating the
      wrong buffer size we send the sizes first if the size returned by
      pack_size is the size of an integer or if it is greater than the
      maximum send_size(to that particular processor) . In this
      approach - in the worst case we might end up sending two messages
      always .

    */
    total_size = 0 ;
    const int nsend = send_info.size() ;
    entitySet resend_procs, rerecv_procs ;
    for(int i=0;i<nsend;++i) {
      s_size[i] = 0 ;
      for(size_t j=0;j<send_info[i].second.size();++j) {
	storeRepP sp = send_vars[i][j] ; //facts.get_variable(send_info[i].second[j].v) ;
        s_size[i] += sp->pack_size(send_info[i].second[j].set) ;
      }
      if((s_size[i] > maxs_size[i]) || (s_size[i] == sizeof(int))) {
	if(s_size[i] > maxs_size[i])
	  maxs_size[i] = s_size[i] ;
	int proc = send_info[i].first ;
	s_size[i] = sizeof(int) ;
	resend_procs += proc ;
	send_index.push_back(i) ;
      }
      total_size += maxs_size[i] ;
    }
    if(send_ptr_buf_size < total_size) {
      if(send_ptr_buf)
        delete[] send_ptr_buf ;
      send_ptr_buf = new unsigned char[total_size] ;
      send_ptr_buf_size = total_size ;
    }
    send_ptr[0] = send_ptr_buf ;
    for(int i = 1; i < nsend; i++)
      send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
    // Pack the buffer for sending
    for(int i=0;i<nsend;++i) {
      int loc_pack = 0 ;
      if(!resend_procs.inSet(send_info[i].first)) {
	for(size_t j=0;j<send_info[i].second.size();++j) {
	  storeRepP sp = send_vars[i][j] ; //facts.get_variable(send_info[i].second[j].v) ;
	  sp->pack(send_ptr[i], loc_pack,s_size[i],send_info[i].second[j].set);
	}
      }
      else
	MPI_Pack(&maxs_size[i], sizeof(int), MPI_BYTE, send_ptr[i], s_size[i], &loc_pack, MPI_COMM_WORLD) ;
    }
    // Send Buffer
    for(int i=0;i<nsend;++i) {
      int proc = send_info[i].first ;
      MPI_Send(send_ptr[i],s_size[i],MPI_PACKED,proc,1,MPI_COMM_WORLD) ;
    }
    /* We receive a message from all the processors in the
       neighbourhood. Whether the message needs to be received a second
       time is determined from the size of the message received. If the
       size of the message is equal to the size of an integer it is added
       to the list to be received a second time(even if the sent value is
       a store value or the size of the message to be received. */
    if(nrecv > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(nrecv, request, status) ;
      FATAL(err != MPI_SUCCESS) ;
      for(int i = 0 ; i < nrecv; i++) {
	MPI_Get_count(&status[i], MPI_BYTE, &recv_sizes[i]) ;
	if(recv_sizes[i] == sizeof(int)) {
	  rerecv_procs += recv_info[i].first ;
	  recv_index.push_back(i) ;
	}
      }
    }
    for(int i=0;i<nrecv;++i) {
      int loc_unpack = 0;
      if(rerecv_procs.inSet(recv_info[i].first)) {
	int temp ;
	/*If the size of the message received is that of an integer
	  then we need to check whether it is greater than the maximum
	  size received so far from that processor. If it is not then
	  the maximum size is set to that value. */
	MPI_Unpack(recv_ptr[i], r_size[i], &loc_unpack, &temp, sizeof(int), MPI_BYTE, MPI_COMM_WORLD) ;
	if(temp > maxr_size[i])
	  maxr_size[i] = temp ;
      }
      else
	for(size_t j=0;j<recv_info[i].second.size();++j) {
	  storeRepP sp = recv_vars[i][j] ; // facts.get_variable(recv_info[i].second[j].v) ;
	  sp->unpack(recv_ptr[i], loc_unpack, r_size[i],
		     recv_info[i].second[j].seq) ;
	}
    }
    rerecv_size = rerecv_procs.size() ;
    resend_size = resend_procs.size() ;
    if(rerecv_size > 0) {
      re_request =  new MPI_Request[rerecv_size] ;
      re_status =  new MPI_Status[rerecv_size] ;
    }
    for(int i = 0; i < rerecv_size; i++) {
      int proc = recv_info[recv_index[i]].first ;
      recv_ptr[recv_index[i]] = new unsigned char[maxr_size[recv_index[i]]] ;
      MPI_Irecv(recv_ptr[recv_index[i]], maxr_size[recv_index[i]], MPI_PACKED, proc, 2, MPI_COMM_WORLD, &re_request[i]) ;
    }

    for(int i=0;i<resend_size;++i) {
      int loc_pack = 0 ;
      send_ptr[send_index[i]] = new unsigned char[maxs_size[send_index[i]]] ;
      for(size_t j=0;j<send_info[send_index[i]].second.size();++j) {
	storeRepP sp = send_vars[send_index[i]][j] ; //facts.get_variable(send_info[send_index[i]].second[j].v) ;
	sp->pack(send_ptr[send_index[i]], loc_pack,maxs_size[send_index[i]],send_info[send_index[i]].second[j].set);
      }
    }
    // Send Buffer
    for(int i=0;i<resend_size;++i) {
      int proc = send_info[send_index[i]].first ;
      MPI_Send(send_ptr[send_index[i]],maxs_size[send_index[i]],MPI_PACKED,proc,2,MPI_COMM_WORLD) ;
      delete [] send_ptr[send_index[i]] ;
    }
    if(rerecv_size > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(rerecv_size, re_request, re_status) ;
      FATAL(err != MPI_SUCCESS) ;
    }
    for(int i=0;i<rerecv_size;++i) {
      int loc_unpack = 0;
      for(size_t j=0;j<recv_info[recv_index[i]].second.size();++j) {
	//vset += recv_info[recv_index[i]].second[j].v ;
	storeRepP sp = recv_vars[recv_index[i]][j] ; //facts.get_variable(recv_info[recv_index[i]].second[j].v) ;
	sp->unpack(recv_ptr[recv_index[i]], loc_unpack, maxr_size[recv_index[i]],
		   recv_info[recv_index[i]].second[j].seq) ;
      }
      delete [] recv_ptr[recv_index[i]] ;
    }
    if(rerecv_size > 0) {
      delete [] re_status ;
      delete [] re_request ;
    }

    double tt = s.stop() ;
    timer.addTime(tt,1) ;
  }

  void execute_comm::Print(ostream &s) const {
    int sz = 0 ;
    if(send_info.size()+recv_info.size() > 0) {
      printIndent(s) ;
      s << "communication block {" << endl ;
      if(send_info.size() > 0) {
        printIndent(s) ;
        s << "Send:" << endl ;
        printIndent(s) ;
        for(size_t i=0;i<send_info.size();++i) {
          for(size_t j=0;j<send_info[i].second.size();++j) {
            s << send_info[i].second[j].v << ' ' ;
	    sz += (send_info[i].second[j].set).size() ;
	  }
	  s << " to " << send_info[i].first << endl ;
          printIndent(s) ;
        }
	s << " Total entities sent = " << sz << endl ;
      }
      if(recv_info.size() > 0) {
        printIndent(s) ;
        s << "Recv:" << endl ;
        printIndent(s) ;
        for(size_t i=0;i<recv_info.size();++i) {
          for(size_t j=0;j<recv_info[i].second.size();++j)
            s << recv_info[i].second[j].v << ' '  ;
          s << " from " << recv_info[i].first << endl ;
          printIndent(s) ;
        }
      }
      s << "}" << endl ;
    }
  }

  void execute_comm::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "comm: " ;

    variableSet vars  ;
    for(size_t i=0;i<send_info.size();++i)
      for(size_t j=0;j<send_info[i].second.size();++j) 
        vars += send_info[i].second[j].v ;

    for(size_t i=0;i<recv_info.size();++i) 
      for(size_t j=0;j<recv_info[i].second.size();++j)
        vars += recv_info[i].second[j].v ;

    oss << vars ;

    data_collector.accumulateTime(timer,EXEC_COMMUNICATION,oss.str()) ;
  }


  // Sort the communication list so that the receive sequence is in the
  // order corresponding to the sending entitySet
  list<comm_info> sort_comm(list<comm_info> slist, fact_db &facts) {
    vector<pair<int,vector<send_var_info> > > send_info ;
    vector<pair<int,vector<recv_var_info> > > recv_info ;

    // First collect information from slist

    HASH_MAP(int,vector<send_var_info>) send_data ;
    HASH_MAP(int,vector<recv_var_info>) recv_data ;
    list<comm_info>::const_iterator cli ;
    intervalSet send_procs, recv_procs ;
    for(cli=slist.begin();cli!=slist.end();++cli) {
      variable v = cli->v ;
      if(cli->send_set.size() > 0) {
        int send_proc = cli->processor ;
        send_procs += send_proc ;
        entitySet send_set = cli->send_set ;
	send_data[send_proc].push_back(send_var_info(v,send_set)) ;
      }
      if(cli->recv_set.size() > 0) {
        int recv_proc = cli->processor ;
        sequence recv_seq = cli->recv_set ;
        recv_procs += recv_proc ;
	recv_data[recv_proc].push_back(recv_var_info(v,recv_seq)) ;
      }
    }

    for(intervalSet::const_iterator ii=send_procs.begin();
        ii!=send_procs.end();
        ++ii) {
      send_info.push_back(make_pair(*ii,send_data[*ii])) ;
    }
    for(intervalSet::const_iterator ii=recv_procs.begin();
        ii!=recv_procs.end();
        ++ii) {
      recv_info.push_back(make_pair(*ii,recv_data[*ii])) ;
    }


    // Now build sorted comm list

    list<comm_info> clist ;

    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    Map l2g ;
    l2g = d->l2g.Rep() ;
    store<unsigned char> key_domain ;
    key_domain = d->key_domain.Rep() ;

    const int nrecv = recv_info.size() ;
    int *r_size = new int[nrecv] ;
    int total_size = 0 ;
    for(int i=0;i<nrecv;++i) {
      r_size[i] = 0 ;
      for(size_t j=0;j<recv_info[i].second.size();++j) {
        r_size[i] += 2*recv_info[i].second[j].seq.size() ;
      }
      total_size += r_size[i] ;
    }
    int  **recv_ptr = new int*[max(nrecv,1)] ;
    recv_ptr[0] = new int[total_size] ;
    for(int i=1;i<nrecv;++i)
      recv_ptr[i] = recv_ptr[i-1]+r_size[i-1] ;

    const int nsend = send_info.size() ;
    int *s_size = new int[nsend] ;
    total_size = 0 ;
    for(int i=0;i<nsend;++i) {
      s_size[i] = 0 ;
      for(size_t j=0;j<send_info[i].second.size();++j) {
        s_size[i] += 2*send_info[i].second[j].set.size() ;
      }
      total_size += s_size[i] ;
    }
    int **send_ptr = new int*[max(nsend,1)] ;
    send_ptr[0] = new int[total_size] ;
    for(int i=1;i<nsend;++i)
      send_ptr[i] = send_ptr[i-1]+s_size[i-1] ;

    MPI_Request *request =  new MPI_Request[nrecv] ;
    MPI_Status *status =  new MPI_Status[nrecv] ;

    for(int i=0;i<nrecv;++i) {
      int proc = recv_info[i].first ;
      MPI_Irecv(recv_ptr[i], r_size[i], MPI_INT, proc, 1,
                MPI_COMM_WORLD, &request[i]) ;

    }

    // Pack the buffer for sending
    for(int i=0;i<nsend;++i) {
      int loc_pack = 0 ;
      for(size_t j=0;j<send_info[i].second.size();++j) {
        comm_info ci ;
        ci.v = send_info[i].second[j].v ;
        ci.processor = send_info[i].first ;
        ci.send_set = send_info[i].second[j].set ;
        clist.push_back(ci) ;
        for(entitySet::const_iterator ei=ci.send_set.begin();
            ei!=ci.send_set.end();
            ++ei) {
	  send_ptr[i][loc_pack++] = key_domain[*ei] ;
          send_ptr[i][loc_pack++] = l2g[*ei] ;
	}
      }
      WARN(loc_pack != s_size[i]) ;
    }

    // Send Buffer
    for(int i=0;i<nsend;++i) {
      int proc = send_info[i].first ;
      MPI_Send(send_ptr[i],s_size[i],MPI_INT,proc,1,MPI_COMM_WORLD) ;
    }

    if(nrecv > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(nrecv, request, status) ;
      FATAL(err != MPI_SUCCESS) ;
    }


    for(int i=0;i<nrecv;++i) {
      int loc_unpack = 0;
      for(size_t j=0;j<recv_info[i].second.size();++j) {
        sequence seq ;
        for(size_t k=0;k<recv_info[i].second[j].seq.size();++k) {
	  int kd = recv_ptr[i][loc_unpack++] ;
	  seq += d->g2lv[kd][recv_ptr[i][loc_unpack++]] ;
	}
	comm_info ci ;
        ci.v = recv_info[i].second[j].v ;
        ci.processor = recv_info[i].first ;
        ci.recv_set = seq ;
	clist.push_back(ci) ;
      }
      WARN(loc_unpack != r_size[i]) ;
    }

    delete [] status ;
    delete [] request ;
    delete [] send_ptr[0] ;
    delete [] send_ptr ;
    delete [] s_size ;
    delete [] recv_ptr[0] ;
    delete [] recv_ptr ;
    delete [] r_size ;

    return clist ;
  }

   
  void barrier_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()){
      std::vector<std::pair<variable,entitySet> > send_entities = barrier_existential_rule_analysis(barrier_vars, facts, scheds) ;
      scheds.update_send_entities(send_entities, sched_db::BARRIER);
    }
  }

  void barrier_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
   
    
#ifdef VERBOSE
    Loci::debugout << "entering barrier process requests, " << barrier_vars
                   << endl ;
#endif
    if(facts.isDistributed()) {
      std::list<comm_info> clist;
      std::list<comm_info> plist;
      std::vector<std::pair<variable,entitySet> > send_entities =
        scheds.get_send_entities(barrier_vars, sched_db::BARRIER);
    
      
      vector<pair<variable,entitySet> >::const_iterator vi ;
      vector<pair<variable,entitySet> > send_requested ;
      if(duplicate_work) {
	//Find out which variables are duplicate variables
	set_duplication_of_variables(barrier_vars, scheds, facts);
      }

      list<comm_info> request_comm ;
      /* The list<comm_info> returned by the
	 barrier_process_rule_requests contains the communication
	 information to send and receive the entities in the clone region*/
      request_comm = barrier_process_rule_requests(barrier_vars, facts, scheds) ;

      clist = request_comm ;
      clist = sort_comm(request_comm,facts) ;
      
     
      
      //Find out which entities are to be sent on owner processor
      //based on duplication policies
      if(!duplicate_work) {
	for(vi=send_entities.begin();vi!=send_entities.end();++vi) {
	  variable v = vi->first ;
	  entitySet send_set = vi->second ;
	  send_requested.push_back(make_pair(v,send_set &
					     scheds.get_variable_requests(v))) ;
	}
      }
      else {
	send_requested = send_ent_for_plist(barrier_vars, facts, scheds);
      }

      /*The put_precomm_info is used in case there is a mapping in the
	output for any of the rules. */
      plist = put_precomm_info(send_requested, facts) ;
      
      scheds.update_comm_info_list(clist, sched_db::BARRIER_CLIST);
      scheds.update_comm_info_list(plist, sched_db::BARRIER_PLIST);
    }
#ifdef VERBOSE
    Loci::debugout << "exiting barrier process requests, " << barrier_vars
                   << endl ;
#endif
  }

  executeP barrier_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {

    if(facts.isDistributed()) {
      std::list<comm_info> clist = scheds.get_comm_info_list(barrier_vars, facts, sched_db::BARRIER_CLIST);
      std::list<comm_info> plist = scheds.get_comm_info_list(barrier_vars, facts, sched_db::BARRIER_PLIST);
      
      CPTR<execute_list> el = new execute_list ;

      bool isGPUSync = false ;
      for(variableSet::const_iterator vi=barrier_vars.begin();
	  vi != barrier_vars.end();++vi) {
	if(vi->name.substr(0,7) == string("__GPU__"))
	  isGPUSync = true ;
      }

      // Add gpu syncronization point
      if(isGPUSync) {
	executeP tmp = new execute_gpuSync(barrier_vars) ;
	el->append_list(tmp) ;
      }
      execute_comm2::inc_comm_step() ;
      if(!plist.empty()) {
        executeP tmp2 = new execute_comm2(plist, facts) ;
        el->append_list(tmp2) ;
      }
      execute_comm2::inc_comm_step() ;
      if(!clist.empty()) {
        executeP tmp2 = new execute_comm2(clist, facts) ;
        el->append_list(tmp2) ;
      }

      if(MPI_processes==1 && !isGPUSync) {
	ostringstream oss ;
	oss << "Sync: " << barrier_vars << endl ;
	executeP exec_thrd_sync = new execute_msg(oss.str()) ;
	el->append_list(exec_thrd_sync) ;
      }
      return executeP(el) ;
    }
    ostringstream oss ;
      oss << "Sync: " << barrier_vars << endl ;
    executeP exec_thrd_sync = new execute_msg(oss.str()) ;
    return exec_thrd_sync;
  }

  void execute_msg::execute(fact_db &facts, sched_db& scheds) {  }

  void execute_msg::Print(std::ostream &s) const {
    printIndent(s) ;
    s << msg << endl ;
  }

  void singleton_var_compiler::set_var_existence(fact_db &facts, sched_db &scheds)  {
    if(facts.isDistributed())
      barrier_existential_rule_analysis(barrier_vars, facts, scheds) ;
    if(duplicate_work) {
      for(variableSet::const_iterator vi = barrier_vars.begin();
	  vi != barrier_vars.end(); vi++)
	scheds.add_policy(*vi, sched_db::NEVER);
    }
  }

  void singleton_var_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      barrier_process_rule_requests(barrier_vars, facts, scheds) ;
    }
  }

  executeP singleton_var_compiler::create_execution_schedule(fact_db &facts,
                                                             sched_db &scheds){
    if(verbose) {
      variableSet vars ;
      vars = barrier_vars ;
      ostringstream oss ;
      oss << "singleton param " << vars ;
      executeP execute = executeP(new execute_msg(oss.str())) ;
      return execute;
    }
    return executeP(0) ;
  }

  /////////////////////////////////////////////////////////////////////////
  //////// allocate and deallocate compiler code //////////////////////////
  /////////////////////////////////////////////////////////////////////////
  class execute_allocate_var : public execute_modules {
    variableSet allocate_vars ;
    map<variable,entitySet> v_requests ;
    map<variable,variable> v2alias ;
    map<variable,double> v_max_sizes ;
    timeAccumulator timer ;
  public:
    execute_allocate_var(const variableSet& vars,
                         const map<variable,entitySet> &vr,
                         const map<variable,variable> &va)
      : allocate_vars(vars), v_requests(vr),v2alias(va) {
      for(variableSet::const_iterator vi=allocate_vars.begin();
          vi!=allocate_vars.end();++vi)
        v_max_sizes[*vi] = 0 ;
    }
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() { return "execute_allocate_var";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  void execute_allocate_var::execute(fact_db &facts, sched_db &scheds) {

    stopWatch s ;
    s.start() ;
    for(variableSet::const_iterator vi=allocate_vars.begin();
        vi!=allocate_vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      entitySet alloc_dom = v_requests[*vi] + srp->domain() ;
      
      if(srp->domain() == EMPTY) {
#ifdef VERBOSE
	debugout << "allocate " << *vi << ", alloc_dom =" << alloc_dom << endl ;
#endif
	srp->allocate(alloc_dom) ;
      }else {
#ifdef VERBOSE
	debugout << "reallocate " << *vi << ", alloc_dom =" 
		 << alloc_dom << endl ;
#endif
        if(profile_memory_usage || collect_memory_info) {
          // this variable is reallocated, we take
          // the space off from the counter, since
          // it will be recounted in profiling
          int packsize = srp->pack_size(srp->domain()) ;
          v_max_sizes[*vi] = max(v_max_sizes[*vi],double(packsize)) ;
          LociAppPMTemp -= packsize ;
          LociAppAllocRequestBeanCounting -= packsize ;
          LociAppFreeRequestBeanCounting -= packsize ;
        }
	if(isSTORE(srp)) {

	  entitySet tmp = interval(alloc_dom.Min(), alloc_dom.Max()) ;
	  if(verbose && tmp.size() >= 2*srp->domain().size() )
	    Loci::debugout << "Variable = " << *vi << "  more than twice the space allocated :  allocated over " << alloc_dom << " size = " << tmp.size()  << "  while domain is only  " << srp->domain() << " size = " << srp->domain().size() << endl ;
	  if(alloc_dom != srp->domain()) {
	    if(verbose)
	      Loci::debugout << "reallocating " << *vi << "  over  " 
			     << alloc_dom << " initially it was over  " 
			     << srp->domain() << endl ;
	    srp->allocate(alloc_dom) ;
	  }
	}
      }
    }
    timer.addTime(s.stop(),1) ;
  }

  void execute_allocate_var::Print(std::ostream &s) const {
    if(verbose && allocate_vars != EMPTY) {
      printIndent(s) ;
      s << "allocating variables " << allocate_vars << endl ;
    }
  }

  void execute_allocate_var::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "allocate: "<<allocate_vars ;

    data_collector.accumulateTime(timer,EXEC_CONTROL,oss.str()) ;
    for(variableSet::const_iterator vi=allocate_vars.begin();
        vi!=allocate_vars.end();++vi) {
      ostringstream vname ;
      map<variable,variable>::const_iterator mi=v2alias.find(*vi) ;
      if(mi != v2alias.end())
        vname << mi->second ;
      else
        vname << *vi ;
      string vn = vname.str() ;
      map<variable,double>::const_iterator mi2 = v_max_sizes.find(*vi) ;
      data_collector.accumulateMemory(vn,ALLOC_CREATE,0,mi2->second) ;
    }
  }

  class execute_free_var : public execute_modules {
    variableSet free_vars ;
    map<variable,variable> v2alias ;
    timeAccumulator timer ;
    map<variable,double> v_max_sizes ;

    double max_memory ;
  public:
    execute_free_var(const variableSet& vars,
                     const map<variable,variable> &va)
      :free_vars(vars),v2alias(va) {
      max_memory=0;
#ifdef VERBOSE
      debugout << "free vars: " << free_vars << endl ;
#endif
      for(variableSet::const_iterator vi=free_vars.begin();
          vi!=free_vars.end();++vi) {
        v_max_sizes[*vi] = 0 ;
      }
    }
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() { return "execute_free_var";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  void execute_free_var::execute(fact_db &facts, sched_db &scheds) {
    stopWatch s ;
    s.start() ;

    if(profile_memory_usage || collect_memory_info) {
      for(variableSet::const_iterator vi=free_vars.begin();
          vi!=free_vars.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        if(srp != 0) {
          int packsize = srp->pack_size(srp->domain()) ;
          v_max_sizes[*vi] = max(v_max_sizes[*vi],double(packsize)) ;
        }
      }
      max_memory=max(max_memory,currentMem()) ;
    }
    
    for(variableSet::const_iterator vi=free_vars.begin();
        vi!=free_vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      if(srp != 0) {
	srp->allocate(EMPTY) ;
#ifdef VERBOSE
	debugout << "deallocating " << *vi << endl ;
#endif
      }
    }
    timer.addTime(s.stop(),1) ;
  }

  void execute_free_var::Print(std::ostream &s) const {
    if(verbose && free_vars != EMPTY) {
      printIndent(s) ;
      s << "deallocating variables " << free_vars << endl ;
    }
  }

  void execute_free_var::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss <<"freevar: "<<free_vars ;

    string s = oss.str() ;
    data_collector.accumulateTime(timer,EXEC_CONTROL,s) ;
    for(variableSet::const_iterator vi=free_vars.begin();
        vi!=free_vars.end();++vi) {
      ostringstream vname ;
      map<variable,variable>::const_iterator mi=v2alias.find(*vi) ;
      if(mi != v2alias.end())
        vname << mi->second ;
      else
        vname << *vi ;
      string vn = vname.str() ;
      map<variable,double>:: const_iterator mi2 = v_max_sizes.find(*vi) ;
      data_collector.accumulateMemory(vn,ALLOC_DELETE,max_memory,
                                      mi2->second) ;
    }
  }

  void allocate_var_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void allocate_var_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP allocate_var_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    variableSet::const_iterator vi,vii ;

    map<variable,entitySet> v_requests ;
    map<variable,variable> v2alias ;
    for(vi=allocate_vars.begin();vi!=allocate_vars.end();++vi) {
      variableSet aliases = variableSet(scheds.get_aliases(*vi)+
                                        scheds.get_antialiases(*vi)+
                                        scheds.get_synonyms(*vi)+
                                        scheds.get_rotations(*vi)) ;
      if(aliases.size() > 1) {
        // If it looks like there are aliases, then collect all of the
        // information about name aliasing
        // Note: Not sure if there is still a problem due to the one-way
        // nature of aliasing....
        variableSet work ;
        for(vii=aliases.begin();vii!=aliases.end();++vii) {
          work += scheds.get_aliases(*vii) ;
          work += scheds.get_antialiases(*vii) ;
          work += scheds.get_synonyms(*vii) ;
          work += scheds.get_rotations(*vii) ;
        }
        work -= aliases ;
        aliases += work ;

        while(work != EMPTY) {
          variableSet new_vars ;

          for(vii=work.begin();vii!=work.end();++vii) {
            new_vars += scheds.get_aliases(*vii) ;
            new_vars += scheds.get_antialiases(*vii) ;
            new_vars += scheds.get_synonyms(*vii) ;
            new_vars += scheds.get_rotations(*vii) ;
          }
          work=new_vars ;
          work-=aliases ;
          aliases += work ; ;

        }
      }
#ifdef VERBOSE
      debugout << "allocating v=" << *vi << ", aliases = " << aliases << endl ;
#endif
      if(aliases.size() <=1)
        v2alias[*vi] = *vi ;
      else
        v2alias[*vi] = *aliases.begin() ;
      entitySet requests ;
      for(vii=aliases.begin();vii!=aliases.end();++vii) {
	requests += scheds.get_variable_requests(*vii) ;
      }
      v_requests[*vi] = requests ;
    }
    executeP execute = executeP(new execute_allocate_var(allocate_vars,
                                                         v_requests,
                                                         v2alias)) ;
    return execute;
  }

  void free_var_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }


  void free_var_compiler::process_var_requests(fact_db &facts, sched_db &scheds) { }

  executeP free_var_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    variableSet::const_iterator vi,vii ;
    map<variable,variable> v2alias ;
    for(vi=free_vars.begin();vi!=free_vars.end();++vi) {
      variableSet aliases = variableSet(scheds.get_aliases(*vi)+
                                        scheds.get_antialiases(*vi)+
                                        scheds.get_synonyms(*vi)+
                                        scheds.get_rotations(*vi)) ;
      if(aliases.size() > 1) {
        // If it looks like there are aliases, then collect all of the
        // information about name aliasing
        // Note: Not sure if there is still a problem due to the one-way
        // nature of aliasing....
        variableSet work ;
        for(vii=aliases.begin();vii!=aliases.end();++vii) {
          work += scheds.get_aliases(*vii) ;
          work += scheds.get_antialiases(*vii) ;
          work += scheds.get_synonyms(*vii) ;
          work += scheds.get_rotations(*vii) ;
        }
        work -= aliases ;
        aliases += work ;

        while(work != EMPTY) {
          variableSet new_vars ;

          for(vii=work.begin();vii!=work.end();++vii) {
            new_vars += scheds.get_aliases(*vii) ;
            new_vars += scheds.get_antialiases(*vii) ;
            new_vars += scheds.get_synonyms(*vii) ;
            new_vars += scheds.get_rotations(*vii) ;
          }
          work=new_vars ;
          work-=aliases ;
          aliases += work ; ;

        }
      }
      if(aliases.size() <=1)
        v2alias[*vi] = *vi ;
      else
        v2alias[*vi] = *aliases.begin() ;
    }

    executeP execute = executeP(new execute_free_var(free_vars,v2alias)) ;
    return execute;
  }

  /////////////////////////////////////////////////////////////////////////
  //////////////////    memory profiling compiler code       //////////////
  /////////////////////////////////////////////////////////////////////////
  void execute_memProfileAlloc::Print(std::ostream &s) const {
    if(vars != EMPTY) {
      printIndent(s) ;
      s << "memory profiling check point (allocate: " << vars
        << ")" << endl ;
    }
  }

  void execute_memProfileAlloc::dataCollate(collectData &data_collector) const {
  }

  void execute_memProfileAlloc::execute(fact_db& facts, sched_db &scheds) {
    for(variableSet::const_iterator vi=vars.begin();
        vi!=vars.end();++vi) {
      //cerr<<"memory profiling (allocation) on: "<<*vi<<endl;
      storeRepP srp = facts.get_variable(*vi) ;
      entitySet alloc_dom = srp->domain() ;
      LociAppPeakMemory = max(LociAppPeakMemory,currentMem()) ;
      //      Loci::debugout << "LociAppPeakMemory="<< LociAppPeakMemory
      //                     << ", mem=" << mem
      //                     << endl ;

      int packsize = srp->pack_size(alloc_dom) ;
      LociAppAllocRequestBeanCounting += packsize ;
      LociAppPMTemp += packsize ;

      if(LociAppPMTemp > LociAppPeakMemoryBeanCounting) {
        LociAppPeakMemoryBeanCounting = LociAppPMTemp ;
      }
      if(packsize > LociAppLargestAlloc) {
        LociAppLargestAlloc = packsize ;
        LociAppLargestAllocVar = *vi ;
      }
    }
  }

  void execute_memProfileFree::Print(std::ostream &s) const {
    if(vars != EMPTY) {
      printIndent(s) ;
      s << "memory profiling check point (free: " << vars
        << ")" << endl ;
    }
  }

  void execute_memProfileFree::dataCollate(collectData &data_collector) const {
  }

  void execute_memProfileFree::execute(fact_db& facts, sched_db &scheds) {
    for(variableSet::const_iterator vi=vars.begin();
        vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      entitySet alloc_dom = srp->domain() ;

      int packsize = srp->pack_size(alloc_dom) ;
      LociAppFreeRequestBeanCounting += packsize ;
      LociAppPMTemp -= packsize ;
      if(LociAppPMTemp < 0)
        std::cout << "MEMORY PROFILING WARNING: negative memory size"
                  << endl ;
      if(packsize > LociAppLargestFree) {
        LociAppLargestFree = packsize ;
        LociAppLargestFreeVar = *vi ;
      }
    }
  }

  executeP memProfileAlloc_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    executeP execute = executeP(new execute_memProfileAlloc(vars)) ;
    return execute;
  }

  executeP memProfileFree_compiler::create_execution_schedule (fact_db &facts, sched_db &scheds) {
    executeP execute = executeP(new execute_memProfileFree(vars));
    return execute;
  }

  //Finds the plist that contains information of entities need to be sent on the
  //owner processor based on the duplication policies.
  //Also defines which variables are duplicate variables.
  std::vector<std::pair<variable,entitySet> >
  send_ent_for_plist(variableSet vlst, fact_db &facts, sched_db &scheds) {
    vector<pair<variable,entitySet> > send_entities ;
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    vector<entitySet> exinfo ;
    vector<variable> vars ;
    vector<ruleSet> rules;
    vector<variable> send_vars ;

    for(variableSet::const_iterator vi=vlst.begin();vi!=vlst.end();++vi) {
      ruleSet r = scheds.get_existential_rules(*vi) ;
      vars.push_back(*vi);
      rules.push_back(r);
    }

    for(size_t i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;
      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
        if(rule_has_mapping_in_output(*rsi)) {
	  send_vars.push_back(v) ;
	  exinfo.push_back(scheds.get_my_proc_able_entities(v, *rsi)) ;
	}
      }
    }

    map<variable,entitySet> vmap ;
    for(size_t i=0;i<send_vars.size();++i) {
      variable v = send_vars[i] ;
      entitySet send_ents;
      if(!scheds.is_duplicate_variable(v))
	send_ents = exinfo[i] - d->my_entities;

      vmap[v] += send_ents ;
    }

    for(map<variable,entitySet>::const_iterator mi = vmap.begin();
	mi != vmap.end(); mi++)
      send_entities.push_back(make_pair(mi->first,mi->second & scheds.get_variable_requests(mi->first)));
    return send_entities;
  }

  entitySet sending_comm_processors(entitySet sendSet, fact_db &facts) {
    entitySet send_procs;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info();
      for(unsigned int i = 0; i < d->copy.size(); i++) {
	if((d->copy[i].entities & sendSet) != EMPTY)
	  send_procs += d->copy[i].proc;
      }
    }
    return send_procs;
  }

  //Based on the policy selected for a variable,
  //it sets the duplication
  bool process_policy_duplication(variable v, sched_db &scheds, fact_db &facts) {
    if(!scheds.is_policy(v, sched_db::NEVER)) {
      if(scheds.is_policy(v, sched_db::ALWAYS)) {
	//scheds.set_duplicate_variable(v, true);
	return true;
      }
      else if(scheds.is_policy(v, sched_db::MODEL_BASED)){
	double original_comm_time = 0, duplication_comm_time = 0;
	double original_comp_time = 0, duplication_comp_time = 0;

	fact_db::distribute_infoP d = facts.get_distribute_info();
	ruleSet r = scheds.get_existential_rules(v);
	bool reduction = false;
	for(ruleSet::const_iterator ri = r.begin();
	    ri != r.end(); ri++)
	  if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT)
	    reduction = true;

	entitySet initial_requests = scheds.get_variable_requests(v);
	entitySet proc_able;
	entitySet my_proc_able;
	entitySet shadow;
	if(!reduction) {
	  for(ruleSet::const_iterator ri = r.begin();
	      ri != r.end(); ri++) {
	    proc_able += scheds.get_proc_able_entities(v, *ri);
	    my_proc_able += scheds.get_my_proc_able_entities(v, *ri);
	  }
	}
	else {
	  proc_able += scheds.get_reduce_proc_able_entities(v);
	  shadow += scheds.get_variable_shadow(v);
	}

	entitySet original_clist_entities = initial_requests - d->my_entities;
	entitySet duplication_clist_entities = (initial_requests - d->my_entities) - proc_able;
	entitySet original_clist_procs = sending_comm_processors(original_clist_entities, facts);

	entitySet duplication_clist_procs = sending_comm_processors(duplication_clist_entities, facts);
	vector<entitySet> temp = send_entitySetv(original_clist_entities, facts);
	for(int i = 0; i < Loci::MPI_processes; i++) {
	  if(temp[i].size() > 0) {
	    original_clist_procs += i;
	    original_clist_entities += temp[i];
	  }
	}

	temp = send_entitySetv(duplication_clist_entities, facts);
	for(int i = 0; i < Loci::MPI_processes; i++) {
	  if(temp[i].size() > 0) {
	    duplication_clist_procs += i;
	    duplication_clist_entities += temp[i];
	  }
	}

	storeRepP vRep = facts.get_variable(v);
	double original_clist_comm = 0, duplication_clist_comm = 0;
	if(original_clist_entities != EMPTY) {
 	  int size = vRep->pack_size(original_clist_entities);
 	  original_clist_comm = (scheds.get_comm_model().ts*original_clist_procs.size()) +
	    scheds.get_comm_model().tw*size;
 	}

	if(duplication_clist_entities != EMPTY) {
 	  int size = vRep->pack_size(duplication_clist_entities);
 	  duplication_clist_comm = (scheds.get_comm_model().ts*duplication_clist_procs.size()) +
	    scheds.get_comm_model().tw*size;
 	}

	entitySet recv_requests = send_entitySet(initial_requests , facts);
	entitySet requests;
	requests = initial_requests + recv_requests;
	requests += fill_entitySet(requests, facts);

	entitySet original_plist_entities;
	if(!reduction)
	  original_plist_entities = (requests - d->my_entities) & my_proc_able;
	else
	  original_plist_entities = (requests - d->my_entities) & shadow;

	entitySet original_plist_procs = sending_comm_processors(original_plist_entities, facts);

	temp = send_entitySetv(original_plist_entities, facts);
	for(int i = 0; i < Loci::MPI_processes; i++) {
	  if(temp[i].size() > 0) {
	    original_plist_procs += i;
	    original_plist_entities += temp[i];
	  }
	}

	double original_plist_comm = 0;
	if(original_plist_entities != EMPTY) {
	  int size = vRep->pack_size(original_plist_entities);
	  original_plist_comm = scheds.get_comm_model().ts*original_plist_procs.size() + scheds.get_comm_model().tw*size;
	}

	if(original_clist_comm > 0)
	  original_comm_time += original_clist_comm;
	if(duplication_clist_comm > 0)
	  duplication_comm_time += duplication_clist_comm;
	if(original_plist_comm > 0)
	  original_comm_time += original_plist_comm;

 	for(ruleSet::const_iterator ri = r.begin();
 	    ri != r.end(); ri++) {
	  entitySet original_context;
	  entitySet comp_context;
	  vdefmap original_tvarmap;
	  vdefmap comp_tvarmap;

	  if(reduction){
	    FATAL(ri->targets().size() > 1);
	    if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
	      original_tvarmap[v] = (requests &
				     scheds.get_my_proc_able_entities(v, *ri));
	      comp_tvarmap[v] = (initial_requests + duplication_clist_entities)
		& scheds.get_proc_able_entities(v, *ri);
	    }
	    else {
	      original_tvarmap[v] = requests;
	      comp_tvarmap[v] = initial_requests + duplication_clist_entities;
	      entitySet reduce_filter;
	      if(multilevel_duplication)
		reduce_filter = d->comp_entities;
	      else
		reduce_filter = d->my_entities;
	      if(scheds.is_reduction_outputmap(v))
		comp_tvarmap[v] &= reduce_filter;
	      else
		comp_tvarmap[v] &= (reduce_filter + scheds.get_reduce_proc_able_entities(v));
	    }
	  }
	  else {
	    original_tvarmap[v] = requests & scheds.get_my_proc_able_entities(v, *ri);
	    comp_tvarmap[v] = (initial_requests + duplication_clist_entities)
	      & scheds.get_proc_able_entities(v, *ri);

	  }

	  const rule_impl::info &rinfo = ri->get_info().desc ;
	  bool mapping_in_output = false ;
	  set<vmap_info>::const_iterator si ;
	  for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
	    // Transform the variable requests using the mapping constructs
	    // in *si
	    //The context is the union
	    if((si->var).inSet(v)) {
	      // Check to see if there is mapping in the output
	      if(si->mapping.size() != 0)
		mapping_in_output = true ;

	      if(si->var.size() > 1) {
		variableSet::const_iterator vi ;
		for(vi = si->var.begin(); vi != si->var.end(); vi++) {
		  if(*vi != v) {
		    comp_tvarmap[*vi] = comp_tvarmap[v];
		    original_tvarmap[*vi] = original_tvarmap[v];
		  }
		}
	      }
	      original_context |= vmap_target_requests(*si, original_tvarmap,facts, scheds, false) ;
	      comp_context |= vmap_target_requests(*si,comp_tvarmap,facts, scheds, false) ;
	    }
	  }

	  if(reduction)
	    if(ri->get_info().rule_impl->get_rule_class() != rule_impl::UNIT)
	      original_context &= d->my_entities;

	  if(mapping_in_output || (reduction && ri->get_info().rule_impl->get_rule_class() != rule_impl::UNIT)) {
	    entitySet sources = ~EMPTY ;
	    entitySet constraints = ~EMPTY ;
	    set<vmap_info>::const_iterator si ;
	    for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si)
	      sources &= vmap_source_exist(*si,facts, scheds) ;

	    for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
	      constraints &= vmap_source_exist(*si,facts, scheds) ;

	    sources &= constraints ;
	    original_context &= sources ;
	    comp_context &= sources;
	  }

	  double ts, tw;
	  scheds.get_comp_model(*ri).get_parameters(ts, tw);
	  if(scheds.get_comp_model(*ri).is_valid_val(ts)) {
	    original_comp_time += ts;
	    duplication_comp_time += ts;
	    if(scheds.get_comp_model(*ri).is_valid_val(tw)) {
	      if(original_context.size() > 0)
		original_comp_time += tw*original_context.size();

	      if(comp_context.size() > 0)
		duplication_comp_time += tw*comp_context.size();
	    }
	    else if(comp_context.size() > 0 || original_context.size() > 0) {
	      cerr << "Error: rule " << " has model problems." << endl;
	      cerr << "Value of tw is invalid for a rule." << endl;
	      cerr << "It may be because no linear regression is performed using data of a previous run" << endl;
	    }
	  }
	  else {
	    cerr << "Error: rule " << r << " is not in the model."  << endl;
	  }
	}

	if(duplication_comp_time < 0)
	  duplication_comp_time = 0;

	if(duplication_comp_time < 0)
	  original_comp_time = 0;

	scheds.add_original_communication_time(v, original_comm_time);
	scheds.add_original_computation_time(v, original_comp_time);
	scheds.add_duplication_communication_time(v, duplication_comm_time);
	scheds.add_duplication_computation_time(v, duplication_comp_time);

	double time[4];
	double max_time[4];
	time[0] = original_comm_time;
	time[1] = original_comp_time;
	time[2] = duplication_comm_time;
	time[3] = duplication_comp_time;

	MPI_Allreduce(time, max_time, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if((max_time[0] + max_time[1]) - (max_time[2] + max_time[3])
	   > -0.000000000000000001)
	  //scheds.set_duplicate_variable(v, true);
	  return true;
	else
	  return false;
      }
      else
	return false;
    }
    else
      return false;
  }

  //It considers all variables which are associated with rules that
  //compute tvars, and figures out if they are duplicate variables
  void set_duplication_of_variables(variableSet tvars, sched_db &scheds,
				    fact_db &facts) {
    vector<variable> vars ;
    vector<ruleSet> rules;
    for(variableSet::const_iterator vi=tvars.begin();vi!=tvars.end();++vi) {
      variable v = *vi ;
      ruleSet r = scheds.get_existential_rules(v) ;
      vars.push_back(v) ;
      rules.push_back(r);
    }

    variableSet current_possible_duplicate_vars, future_possible_duplicate_vars;
    current_possible_duplicate_vars += tvars & scheds.get_possible_duplicate_vars();
    //First add variables that have mapping in input or output
    for(size_t i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;
      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
	if(rule_has_mapping_in_output(*rsi)) {
          current_possible_duplicate_vars += v;
        }
        current_possible_duplicate_vars += tvars & (input_variables_with_mapping(*rsi));
	future_possible_duplicate_vars += input_variables_with_mapping(*rsi) - tvars;
      }
    }

    //Figure out duplication of variables which are subset of tvars
    for(variableSet::const_iterator vi = current_possible_duplicate_vars.begin();
	vi != current_possible_duplicate_vars.end(); vi++)
      scheds.set_duplicate_variable(*vi, process_policy_duplication(*vi, scheds, facts));

    //Now if target variable is duplicate variable, add the variables
    //those are in input of the rule which compute that variable
    for(size_t i=0;i<vars.size();++i) {
      variable v = vars[i] ;
      ruleSet &rs = rules[i] ;
      if(scheds.is_duplicate_variable(v))
	for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi)
	  future_possible_duplicate_vars += input_variables(*rsi);
    }

    scheds.add_possible_duplicate_vars(future_possible_duplicate_vars);
  }

} // end of namespace Loci

