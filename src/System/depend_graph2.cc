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
#include <depend_graph.h>
#include "dist_tools.h"
#include "loci_globs.h"
#include <unistd.h>
#include <map>
using std::map ;
#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <string>
using std::string ;
#include <sstream>
using std::ostringstream ;
#include <ostream>

#include <iostream>
using std::endl ;
using std::cerr ;
using std::cout ;
#include <algorithm>
using std::set_intersection ;
#include <iterator>
using std::back_inserter ;

//#define VERBOSE

namespace Loci {

  namespace {
    rule create_rule(variable sv, variable tv, string qualifier) {
      ostringstream oss ;
      oss << "source(" << sv << ')' ;
      oss << ",target(" << tv << ')' ;
      oss << ",qualifier(" << qualifier << ')' ;
      string sig = oss.str() ;
      rule r(sig) ;
      return r ;
    }
    
    rule create_rule(variableSet source, variableSet target,
                     string qualifier) {
      ostringstream oss ;
      oss << "source(" << source << ')' ;
      oss << ",target(" << target << ')' ;
      oss << ",qualifier(" << qualifier << ")" ;
      string sig = oss.str() ;
      rule r(sig) ;
      return r ;
    }
    
    inline void invoke_rule(rule f, digraph &gr) {
      gr.add_edges(f.sources(),f.ident()) ;
      gr.add_edges(f.ident(),f.targets()) ;
    }

    void invoke_rule_wp(rule f, digraph &gr) {
      variableSet targets = f.targets() ;
      variableSet sources = f.sources() ;
      invoke_rule(f,gr) ;
      for(variableSet::const_iterator vi=targets.begin();
          vi!=targets.end();
          ++vi) {
        if(vi->get_info().priority.size() != 0) {
          variable v = *vi ;
          variable vold = v ;
          while(v.get_info().priority.size() != 0) {
            v = v.drop_priority() ;
            rule priority_rule = create_rule(vold,v,"priority") ;
            invoke_rule(priority_rule,gr) ;
            vold = v ;
          }
        }
      }
    }
    
    variableSet convert_stationary(const variableSet &v) {
      variableSet result ;
      variableSet::const_iterator vi ;
      time_ident stationary_time ;
      for(vi=v.begin();vi!=v.end();++vi) {
        if(vi->time() != stationary_time)
          result += variable(*vi,time_ident()) ;
        else
          result += *vi ;
      }
      return result ;
    }

    variableSet convert_time(const variableSet &v, time_ident t) {
      variableSet result ;
      variableSet::const_iterator vi ;
      for(vi=v.begin();vi!=v.end();++vi) {
        if(vi->time() != t)
          result += variable(*vi,t) ;
        else
          result += *vi ;
      }
      return result ;
    }
    
    inline bool time_before(time_ident t1, time_ident t2) {
      return t1.before(t2) ;
    }
    
    inline bool time_equal(time_ident t1, time_ident t2) {
      return (!time_before(t1,t2) && !time_before(t2,t1)) ;
    }
    
    inline variable drop_all_priorities(variable v) {
      while(v.get_info().priority.size() != 0)
        v = v.drop_priority() ;
      return v ;
    }

    inline variableSet drop_all_priorities(variableSet vset) {
      variableSet vnew ;
      for(variableSet::const_iterator vi=vset.begin();vi!=vset.end();++vi)
        vnew += drop_all_priorities(*vi) ;
      return vnew ;
    }
    
    // forward declaration
    struct iteration_info ;
    
    struct iteration {
      bool active ;
      digraph iteration_graph ;
      time_ident iteration_time ;
      rule iteration_rule ;
      ruleSet build, advance, collapse, time_specific ;
      variableSet changing_vars ; // variables that changing in this level
      variableSet dont_promote ; // variables don't need to be promoted
      variableSet requests ; // pass to the parent level
      variableSet search_requests ; // used to start building the graph
      void init(iteration_info&,const digraph&,const digraph&) ;
      iteration() {active=false ;}
    } ;

    struct iteration_info {
      map<time_ident,iteration> iteration_rules ;
      map<rule,time_ident> iteration_time_ident ;
      ruleSet iteration_set ;
    } ;

    // function that classifies all rules in the rule database
    // all iteration related rules are put in an iteration_info
    // object, all other rules are put in a ruleSet
    void classify_ruledb(const ruleSet& all_rules,
                         iteration_info& iter,
                         ruleSet& working_rules) {
      ruleSet::const_iterator ri ;
      ruleSet time_specific ;
      for(ri=all_rules.begin();ri!=all_rules.end();++ri) {
        const rule::rule_type rtype = ri->type() ;
        if(rtype == rule::BUILD) {
          if(ri->target_time().parent()  != ri->source_time()) {
            cerr << "ERROR: malformed build rule, time levels should increment only one level" << endl
                 << *ri << endl ;
          }          
          iter.iteration_rules[ri->target_time()].build += *ri ;
        } else if(rtype == rule::COLLAPSE) {
          variableSet cond = ri->get_info().desc.conditionals ;
          if(cond.size() != 1) {
            cerr << "ERROR: malformed collapse rule with no conditionals:"
                 << endl
                 << *ri << endl ;
          } else if(ri->target_time()  != ri->source_time().parent()) {
            cerr << "ERROR: malformed collapse rule, time levels should increment only one level" << endl
                 << *ri << endl ;
          }          
            iter.iteration_rules[ri->source_time()].collapse += *ri ;
        } else if(!ri->time_advance) {
          working_rules += *ri ;
        }
        if(ri->time_advance) {
          if(rtype != rule::COLLAPSE)
            iter.iteration_rules[ri->target_time()].advance += *ri ;
        }
        if(ri->target_time() != time_ident()) {
          if(rtype != rule::BUILD && rtype != rule::COLLAPSE &&
             !ri->time_advance) {
            time_specific += *ri ;
          }
        }
      }

      for(ri=time_specific.begin();ri!=time_specific.end();++ri) {
        map<time_ident,iteration>::iterator ii ;
        ii = iter.iteration_rules.find(ri->target_time()) ;
        if(ii != iter.iteration_rules.end()) {
          ii->second.time_specific += *ri ;
#ifdef VERBOSE
          debugout << "adding time specific rule " << *ri << endl ;
#endif
        }
      }
    }

    // function that creates a representative rule for
    // each iteration, these representatives are also
    // put into the working_rules set
    void create_iteration_rep(iteration_info& iter,ruleSet& working_rules) {
      map<time_ident,iteration>::iterator mi ;
      for(mi=iter.iteration_rules.begin();mi!=iter.iteration_rules.end();
          ++mi) {
        // set up the time level for this iteration
        mi->second.iteration_time = mi->first ;
        ruleSet build = mi->second.build ;
        ruleSet collapse = mi->second.collapse ;
        variableSet sources, targets ;
        // get the source and target variables for this loop
        ruleSet::const_iterator ri ;
        for(ri=build.begin();ri!=build.end();++ri)
          sources += ri->sources() ;
        for(ri=collapse.begin();ri!=collapse.end();++ri)
          targets += ri->targets() ;
        rule i_rule = create_rule(sources,targets,"iterating_rule") ;
        
        // set up the record ;
        iter.iteration_time_ident[i_rule] = mi->first ;
        mi->second.iteration_rule = i_rule ;
        iter.iteration_set += i_rule ;
        if(i_rule.get_info().time_advance) {
          iter.iteration_rules[i_rule.target_time()].advance += i_rule ;
        }
        else {
          working_rules += i_rule ;
        }
      }
    }

    // function that checks the validity of all the iterations
    // return true if passed, false if failed
    bool check_iteration(const iteration_info& iter) {
      map<time_ident,iteration>::const_iterator mi ;
      bool total_check = true ;
      for(mi=iter.iteration_rules.begin();mi!=iter.iteration_rules.end();
          ++mi) {
        ruleSet build = mi->second.build ;
        ruleSet collapse = mi->second.collapse ;
        ruleSet advance = mi->second.advance ;
        bool check = true ;
        if(build == EMPTY) {
          cerr << "Malformed iteration: " << mi->first
               << ", no build rules exist." << endl ;
          check = false ;
        }
        if(advance == EMPTY) {
          cerr << "Malformed iteration: " << mi->first
               << ", no advance rules exist." << endl ;
          check = false ;
        }
        if(collapse == EMPTY) {
          cerr << "Malformed iteration: " << mi->first
               << ", no collapse rules exist." << endl ;
          check = false ;
        }
        if(check == false) {
          cerr << "iteration rules remaining" << endl
               << "build: " << build  << endl
               << "advance: " << advance << endl
               << "collapse:"<< collapse << endl ;
          total_check = false ;
        }
      }
      // return status
      return total_check ;
    }

    // function prototype for mutual recursion
    rule promote_iterating_rule(const rule&,const time_ident&,
                                iteration_info&) ;
    
    // promote the entire iteration object
    // return the promoted iteration_rule
    rule promote_iteration(const iteration& i,
                           const time_ident& tl,
                           iteration_info& iter) {
      // a new iteration object
      iteration pi ;
      // first generate a new time_ident
      pi.iteration_time = prepend_time(tl,i.iteration_time) ;
      // check if the promotion is already in the record
      map<time_ident,iteration>::const_iterator mi =
        iter.iteration_rules.find(pi.iteration_time) ;
      if(mi!=iter.iteration_rules.end())
        return  pi.iteration_rule ;
      // then promote all the rules
      pi.iteration_rule = prepend_rule(i.iteration_rule,tl) ;
      ruleSet new_build, new_advance, new_collapse, new_time_specific ;
      for(ruleSet::const_iterator ri=i.build.begin();
          ri!=i.build.end();++ri) {
        new_build += prepend_rule(*ri,tl) ;
      }
      for(ruleSet::const_iterator ri=i.time_specific.begin();
          ri!=i.time_specific.end();++ri) {
        new_time_specific += prepend_rule(*ri,tl) ;
      }

      // in advance rule we'll need to check if there
      // are any iteration rules, if any, we also need
      // to promote them to the correct time level
      for(ruleSet::const_iterator ri=i.advance.begin();
          ri!=i.advance.end();++ri)
        if( (ri->type() == rule::INTERNAL) &&
            (ri->qualifier() == "iterating_rule")
            ) {
          new_advance += promote_iterating_rule(*ri,tl,iter) ;
        }else {
          new_advance += prepend_rule(*ri,tl) ;
        }
      
      for(ruleSet::const_iterator ri=i.collapse.begin();
          ri!=i.collapse.end();++ri) {
        new_collapse += prepend_rule(*ri,tl) ;
      }

      pi.build = new_build ;
      pi.advance = new_advance ;
      pi.collapse = new_collapse ;
      pi.time_specific = new_time_specific ;

      // then create a new record in the iteration info object
      iter.iteration_rules[pi.iteration_time] = pi ;
      iter.iteration_time_ident[pi.iteration_rule] = pi.iteration_time ;

      return pi.iteration_rule ;
    }

    // promote iteration rule, return the promoted iteration rule
    rule promote_iterating_rule(const rule& r,
                                const time_ident& tl,
                                iteration_info& iter) {
      map<rule,time_ident>::const_iterator tp =
        iter.iteration_time_ident.find(r) ;
      if(tp==iter.iteration_time_ident.end()) {
        cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
        cerr << "ERROR: iteration rule " << r << " not seen before." << endl ;
        Abort() ;
      }
      map<time_ident,iteration>::iterator ip =
        iter.iteration_rules.find(tp->second) ;
      if(ip==iter.iteration_rules.end()) {
        cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
        cerr << "ERROR: iteration rule record does not exist for rule "
             << tp->second
             << endl ;
        Abort() ;
      }
      return promote_iteration( (ip->second),tl,iter) ;
    }

    inline time_ident
    get_iterating_rule_time(const rule& r,
                            const iteration_info& iter) {
      map<rule,time_ident>::const_iterator tp =
        iter.iteration_time_ident.find(r) ;
      if(tp==iter.iteration_time_ident.end()) {
        cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
        cerr << "ERROR: iteration rule " << r << " not seen before." << endl ;
        Abort() ;
      }
      return tp->second ;
    }

    // return true if t1 contains t2 as a sub-level, false otherwise
    // e.g., t1={n,it}, t2={n} return true;
    //       t1={n,it}, t2={igs} return false ;
    inline bool time_contain(const time_ident& t1,
                             const time_ident& t2) {
      set<string> t1_level ;
      set<string> t2_level ;
      time_ident stationary ;
      time_ident parent_t1 = t1 ;
      while(parent_t1 != stationary) {
        t1_level.insert(parent_t1.level_name()) ;
        parent_t1 = parent_t1.parent() ;
      }
      time_ident parent_t2 = t2 ;
      while(parent_t2 != stationary) {
        t2_level.insert(parent_t2.level_name()) ;
        parent_t2 = parent_t2.parent() ;
      }
      vector<string> diff ;
      set_intersection(t1_level.begin(),t1_level.end(),
                       t2_level.begin(),t2_level.end(),
                       back_inserter(diff)) ;
      if(diff.empty())
        return false ;
      else
        return true ;
    }
    
    // function prototype
    variableSet create_graph(const digraph&,const digraph&,const variableSet&,
                             iteration_info&,digraph&,time_ident) ;
    
    variableSet instantiate_iteration(time_ident,iteration_info&,
                                      const digraph&,const digraph&) ;

    // member function of struct iteration, which
    // initialize necessary iteration information
    void iteration::init(iteration_info& iter,
                         const digraph& rule_graph,
                         const digraph& rule_graph_transpose) {
      active = true ;

      // build generalize rules first
      map<variable,intervalSet> build_offset ;
      for(ruleSet::const_iterator ri=build.begin();
          ri!=build.end();++ri) {
        variableSet build_target = ri->targets() ;
        for(variableSet::const_iterator vi=build_target.begin();
            vi!=build_target.end();++vi) {
          if(!vi->assign) {
            cerr << "incorrect specification of build rule " << *ri
                 << endl
                 << "A correct build rule should have an assign in the target"
                 << " variable.  For example, instead of " << *vi
                 << " use " << vi->name << "{"<< vi->time_id<<"=0}"
                 << endl ;
          }
          variable vt = vi->drop_assign() ;
          variable vbase = vt.new_offset(0) ;
          build_offset[vbase] += vi->offset ;
          // add generalize rule
          rule generalize_rule = create_rule(*vi,vt,"generalize") ;
          invoke_rule(generalize_rule,iteration_graph) ;
        }
      }

      variableSet iteration_input, iteration_output ;
      for(map<variable,intervalSet>::iterator mvi=build_offset.begin();
          mvi!=build_offset.end();++mvi) {
        variable vbase = mvi->first ;
        for(intervalSet::const_iterator ii=mvi->second.begin();
            ii!=mvi->second.end();++ii) {
          variable newv = vbase.new_offset(*ii) ;
          
          iteration_output += newv ;
          changing_vars += drop_all_priorities(newv) ;
        }
        // this is the variable that rotates
        variable newv = vbase.new_offset(mvi->second.Max()+1) ;

        iteration_input += newv ;
        changing_vars += drop_all_priorities(newv) ;
      }
      for(ruleSet::const_iterator ri=build.begin();
          ri!=build.end();++ri) {
        variableSet build_target = ri->targets() ;
        for(variableSet::const_iterator vi=build_target.begin();
            vi!=build_target.end();++vi) {
          changing_vars += drop_all_priorities(*vi) ;
        }
      }

      for(ruleSet::const_iterator ri=time_specific.begin();
          ri!=time_specific.end();++ri) {
        variableSet computed_vars = ri->targets() ;
#ifdef VERBOSE
        debugout << "time specific rule " << *ri << endl ;
#endif
        for(variableSet::const_iterator vi=computed_vars.begin();
            vi!=computed_vars.end();++vi) {
          changing_vars += drop_all_priorities(*vi) ;
        }
      }

      
      // add an output variable
      variable ov("OUTPUT") ;
      variable output(ov,iteration_time) ;
      iteration_input += output ;
      iteration_output += output ;

      // create the time variable
      variable tvar(iteration_time) ;
      iteration_output += tvar ;
      changing_vars += tvar ;

      // create a looping rule that hook up the iteration rule set
      ostringstream oss ;
      oss << "source(" << iteration_input << "),target("
          << iteration_output << "),qualifier(looping)" ;
      rule looping(oss.str()) ;

      // invoke known rules first
      // first the build rules
      for(ruleSet::const_iterator ri=build.begin();
          ri!=build.end();++ri) {
        invoke_rule_wp(*ri,iteration_graph) ;
      }

      // then the collapse rules
      variableSet collapse_inputs ;
      for(ruleSet::const_iterator ri=collapse.begin();
          ri!=collapse.end();++ri) {
        collapse_inputs += ri->sources() ;
        invoke_rule_wp(*ri,iteration_graph) ;
      }
      // finally the advance rules
      // but they need to be treated specially
      // because they themselves could be iterations
      // therefore, we need to instantiate those iteration, if any
      variableSet advance_inputs ;
      for(ruleSet::const_iterator ri=advance.begin();
          ri!=advance.end();++ri) {
        changing_vars += drop_all_priorities(ri->targets()) ;

        if( (ri->type() == rule::INTERNAL) &&
            (ri->qualifier() == "iterating_rule")
            ) {
          // in this case, we must instantiate the iteration
          // we first need to look for the iteration time
          map<rule,time_ident>::const_iterator tp =
            iter.iteration_time_ident.find(*ri) ;
          if(tp==iter.iteration_time_ident.end()) {
            cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
            cerr << "ERROR: iteration rule " << *ri << " not seen before." << endl ;
            Abort() ;
          }
          // instantiate it
          variableSet iteration_requests =
            instantiate_iteration(tp->second,iter,rule_graph,
                                  rule_graph_transpose) ;
          advance_inputs += iteration_requests ;
        } else {
          advance_inputs += ri->sources() ;
          invoke_rule_wp(*ri,iteration_graph) ;
        }
      }
      // finally invoke the looping rule
      invoke_rule(looping,iteration_graph) ;

      search_requests += collapse_inputs ;
      search_requests += advance_inputs ;
      search_requests += iteration_input ;
      search_requests += iteration_output ;

      // search for additional changing variables
      variableSet add_changing ;
      variableSet working_vars = changing_vars ;
      // convert those necessary changing variables to
      // stationary time to searching
      working_vars += convert_stationary(changing_vars) ;
      
      variableSet visited_vars ;
      ruleSet visited_rules ;
      while(working_vars != EMPTY) {
        visited_vars += working_vars ;
        variableSet next ;
        for(variableSet::const_iterator vi=working_vars.begin();
            vi!=working_vars.end();++vi) {
          ruleSet rules = extract_rules(rule_graph[vi->ident()]) ;
          rules -= visited_rules ;
          for(ruleSet::const_iterator ri=rules.begin();
              ri!=rules.end();++ri) {
            variableSet varsPLUStime =
              convert_time(ri->targets(),iteration_time) ;
            variableSet varsATstationary =
              convert_stationary(ri->targets()) ;
            
            add_changing += varsPLUStime ;
            
            next += varsPLUStime ;
            next += varsATstationary ;
          }
          visited_rules += rules ;
        }
        next -= visited_vars ;
        working_vars = next ;
      }
      changing_vars += drop_all_priorities(add_changing) ;

      // for the OUTPUT variable, we need
      // to see if it is in the changing_vars
      // set. If it is, that means rule promotions
      // are allowed, otherwise
      // it is actually not used in the iteration.
      // we disable its promotion
      if(!changing_vars.inSet(output)) {
        dont_promote += output ;
        search_requests -= output ;
      }

      requests += iteration_rule.sources() ;
      // end of the function
#ifdef VERBOSE
      debugout << "Iteration: " << iteration_time << endl ;
      debugout << "changing_vars = " << changing_vars << endl ;
      debugout << "dont_promote = " << dont_promote << endl;
      debugout << "requests = " << requests << endl ;
      debugout << "search_requests = " << search_requests << endl ;
      debugout << "iteration_rule = " << iteration_rule << endl ;
#endif
               
    }

    
    // function that instantiates an iteration,
    // it actually builds the iteration graph.
    // it returns the requests from this time level to
    // its parent level
    variableSet instantiate_iteration(time_ident tlevel,
                                      iteration_info& iter,
                                      const digraph& rule_graph,
                                      const digraph& rule_graph_transpose) {
      // first find out the iteration object
      map<time_ident,iteration>::iterator ip =
        iter.iteration_rules.find(tlevel) ;
      if(ip==iter.iteration_rules.end()) {
        cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
        cerr << "ERROR: iteration rule record does not exist for " << tlevel
             << endl ;
        Abort() ;
      }
      // initialize all the required variables in the object
      iteration& io = ip->second ;
      io.init(iter,rule_graph,rule_graph_transpose) ;
      // then we build the graph
      io.requests += 
      create_graph(rule_graph,rule_graph_transpose,
                   io.search_requests,iter,io.iteration_graph,
                   tlevel) ;

      return io.requests ;
    }

    // function that creates the whole dependency graph
    // it is a recursive process. it instantiates necessary
    // iterations on the fly (see more detail in the
    // dependency graph documentation.).
    // return requests from the passed
    // in time level to its parent level
    variableSet create_graph(const digraph& rule_graph,
                             const digraph& rule_graph_transpose,
                             const variableSet& search_requests,
                             iteration_info& iter, digraph& gr,
                             time_ident tlevel) {
      // setting up useful variables
      variableSet dont_promote, changing_vars ;
      variableSet requests ;
      if(tlevel != time_ident()) {
        map<time_ident,iteration>::iterator ip =
          iter.iteration_rules.find(tlevel) ;
        if(ip==iter.iteration_rules.end()) {
          cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
          cerr << "ERROR: iteration rule record does not exist for " << tlevel
               << endl ;
          Abort() ;
        }
        dont_promote = ip->second.dont_promote ;
        changing_vars = ip->second.changing_vars ;
      }
      
      variableSet visited_vars ;
      ruleSet visited_rules ;
      variableSet working_vars = search_requests ;
      variableSet computed_vars ;
      while(working_vars != EMPTY) {
        visited_vars += working_vars ;
        variableSet next ;
        for(variableSet::const_iterator vi=working_vars.begin();
            vi!=working_vars.end();++vi) {
          ruleSet pre_rules ;
          // if we see any parent level variable, we stop looking
          // and we'll add this to the requests of this iteration
          if(time_equal(vi->time(),tlevel)) {
            pre_rules = extract_rules(rule_graph_transpose[vi->ident()]) ;
          } else if(time_before(vi->time(),tlevel)) {
            requests += *vi ;
            continue ;
          } else {
            // this is an error if (vi->time() after tlevel)
            cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
            cerr << "ERROR: variable time higher than iteration time."
                 << "variable: " << *vi << " iteration time: {"
                 << tlevel  << "}" << endl ;
            //            Abort() ;
          }
          pre_rules -= visited_rules ;
          for(ruleSet::const_iterator ri=pre_rules.begin();
              ri!=pre_rules.end();++ri) {
            if( (ri->type() == rule::INTERNAL) &&
                (ri->qualifier() == "iterating_rule")
                ) {
              // we see an iterating_rule, we need to instantiate it
              
              // we first need to look for the iteration time
              map<rule,time_ident>::const_iterator tp =
                iter.iteration_time_ident.find(*ri) ;
              if(tp==iter.iteration_time_ident.end()) {
                cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
                cerr << "ERROR: iteration rule " << *ri
                     << " not seen before." << endl ;
                Abort() ;
              }
              // instantiate it
              variableSet iteration_requests =
                instantiate_iteration(tp->second,iter,rule_graph,
                                      rule_graph_transpose) ;
              next += iteration_requests ;
            } else {
#ifdef VERBOSE
              debugout << "invoking rule " << *ri << endl ;
#endif
              invoke_rule_wp(*ri,gr) ;
              next += extract_vars(rule_graph_transpose[ri->ident()]) ;
              computed_vars += ri->targets() ;
            }
          }
          // if pre_rules are empty and we are in the stationary time
          // then this *vi is a top level requests because no one
          // at this time can generate *vi
          if(pre_rules == EMPTY) {
            if(tlevel == time_ident()) {
              requests += *vi ;
              continue ;
            }
          }
          visited_rules += pre_rules ;
          if(tlevel != time_ident()) {
            // If a priority variable, then definitely don't promote
            //            if(vi->get_info().priority.size() != 0)
            //              dont_promote += *vi ;
            // then we decide to do rule promotion or variable promotion
            if(!dont_promote.inSet(*vi)) {
              if(changing_vars.inSet(drop_all_priorities(*vi))) {
                // then we do rule promotions, if any
                variable stationary_var(*vi,time_ident()) ;
                ruleSet promote_rules =
                  extract_rules
                  (rule_graph_transpose[stationary_var.ident()]) ;
                // take off any time specific and internal rules
                ruleSet takeoff ;
                for(ruleSet::const_iterator ri=promote_rules.begin();
                    ri!=promote_rules.end();++ri) {
                  if(ri->type() == rule::TIME_SPECIFIC) {
                    takeoff += *ri ;
                  }
                  else if( (ri->type() == rule::INTERNAL) &&
                           (ri->qualifier() == "iterating_rule")) {
                    time_ident itime = get_iterating_rule_time(*ri,iter) ;
                    if(time_contain(tlevel,itime))
                      takeoff += *ri ;
                  }
                }
                promote_rules -= takeoff ;
                for(ruleSet::const_iterator ri=promote_rules.begin();
                    ri!=promote_rules.end();++ri) {
                  if( (ri->type() == rule::INTERNAL) &&
                      (ri->qualifier() == "iterating_rule")) {
                    //cerr<<"promote iteration: "<<*ri<<endl ;
                    rule new_irule = promote_iterating_rule(*ri,tlevel,iter) ;
                    // If successful in promoting, instantiate
                    if(new_irule != rule())
                      next += instantiate_iteration(get_iterating_rule_time
                                                    (new_irule,iter),
                                                    iter,rule_graph,
                                                    rule_graph_transpose) ;
                  }else {
                    rule pr = promote_rule(*ri,tlevel) ;
                    if(!visited_rules.inSet(pr)) {
                      visited_rules += pr ;
                      invoke_rule(pr,gr) ;
                      next += pr.sources() ;
                    }
                  }
                } // end of for(promote)
                
              } else {
                if(!computed_vars.inSet(*vi)) {
                  // we only do variable promotions
                  time_ident parent = tlevel.parent() ;
                  variable pv(*vi,parent) ;
                  invoke_rule(create_rule(pv,*vi,"promote"),gr) ;
#ifdef VERBOSE
                  debugout << "promote variable: " <<  *vi << endl ;
                  debugout << "computed_vars = " << computed_vars << endl ;
#endif
                  requests += pv ;
                }
              }
            } // end of if(!dont_promote.inSet)
          } // end of if(tlevel != time_ident())
        }
        next -= visited_vars ;
        working_vars = next ;
      }

      return requests ;
    }
    
    // function that adds dependency to rename rules
    void add_rename_dependencies(digraph &gr) {
      variableSet all_vars = extract_vars(gr.get_all_vertices()) ;
      ruleSet     all_rules = extract_rules(gr.get_all_vertices()) ;
      
      // extract the qualified rules, these are rules that are
      // automatically generated by the system.  Since these rules
      // cannot provide type information directly, they are
      // singled out so that they can be handled as a separate case.
      ruleSet qualified_rules,rename_rules ;
      for(ruleSet::const_iterator ri=all_rules.begin();
          ri!=all_rules.end();
          ++ri) {
        if(ri->type() == rule::INTERNAL)
          qualified_rules += *ri ;
        set<vmap_info>::const_iterator vmsi ;
        for(vmsi=ri->get_info().desc.targets.begin();
            vmsi!=ri->get_info().desc.targets.end(); ++vmsi)
          if(vmsi->assign.size() != 0) 
            rename_rules += *ri ;
      }
      // We need the transpose of the graph in order to find the rules that
      // generate a particular variable
      
      ruleSet check_rules = all_rules ;
      check_rules -= qualified_rules ;
      for(ruleSet::const_iterator ri=check_rules.begin();
          ri != check_rules.end();
          ++ri) {
        set<vmap_info>::const_iterator vmsi ;
        for(vmsi=ri->get_info().desc.targets.begin();
            vmsi!=ri->get_info().desc.targets.end(); ++vmsi)
          if(vmsi->assign.size() != 0) 
            for(size_t i=0;i<vmsi->assign.size();++i) {
              variable orig_name = vmsi->assign[i].second ;
              //              digraph grt = gr.transpose() ;
              ruleSet depend_rules = extract_rules(gr[orig_name.ident()]) ;
              depend_rules -= rename_rules ;
              // We ignore renaming dependencies in collapse rules, however
              // in reality we may need to follow the dependencies back to
              // the build rules.  However, since we know that the collapse
              // must follow the build, there is no need to add the rename
              // dependencies to collapse rules.
              if(depend_rules.size() != 0 && ri->type() != rule::COLLAPSE) {
                gr.add_edges(depend_rules,ri->ident()) ;
              }
            }
      }
      
    }
  } // end of unnamed namespace
    
  // function that clean the dependency graph at last
  void clean_graph(digraph &gr, const variableSet& given,
                                      const variableSet& target) {

    bool debugging = MPI_processes == 1 || verbose ;
    
    if(verbose) {
      debugout << "given = " << given << endl ;
      debugout << "targets = " << target << endl ;
    }

    // Adjustments for super rule
    ruleSet super_rules ;
    { 
      variable UNIVERSE("UNIVERSE") ;
      ruleSet cmp = extract_rules(gr.get_all_vertices()) ;
      for(ruleSet::const_iterator ri=cmp.begin();ri!=cmp.end();++ri) {
	if(ri->type() != rule::INTERNAL 
	   && ri->get_rule_implP()->get_rule_class() == rule_impl::SUPER_RULE) {
	  debugout << "super_rule " << *ri << endl ;
	  super_rules += *ri ;
	  gr.add_edge(UNIVERSE.ident(),ri->ident()) ;
	}
      }
    }


    // testing...
    //given -= variable("EMPTY") ;
    bool cleaned_rules = false ;
    ruleSet all_cleaned_rules ;

    // First remove any OUTPUT rules that are using values computed in
    // their iteration level.

    {
      digraph::vertexSet allvertices = gr.get_all_vertices() ;
      variableSet allvars = extract_vars(allvertices) ;
      variableSet outputs ;
      digraph grt = gr.transpose() ;
      for(variableSet::const_iterator vi=allvars.begin();
          vi!=allvars.end();++vi) {
        variable ov = *vi ;
        if(ov.get_info().name == "OUTPUT")
          outputs += *vi ;
      }
      ruleSet remove_rules ;
      for(variableSet::const_iterator vi=outputs.begin();
          vi!=outputs.end();++vi) {
        ruleSet outrules = extract_rules(grt[vi->ident()]) ;
        ruleSet outloop = extract_rules(gr[vi->ident()]) ;

        outrules -= outloop ;

        ruleSet::const_iterator ri ;
        for(ri=outrules.begin();ri!=outrules.end();++ri) {
          digraph::vertexSet working, breadth ;
          working += ri->ident() ;
          
          while(working != EMPTY) {
            digraph::vertexSet visit ;
            for(digraph::vertexSet::const_iterator dvi=working.begin();
                dvi != working.end();++dvi)
              visit += grt[*dvi] ;
            ruleSet vrules = extract_rules(visit) ;
            ruleSet::const_iterator vri ;
            for(vri=vrules.begin();vri!=vrules.end();++vri) {
              // Don't follow rules out of the iteration
              if( (vri->type() == rule::INTERNAL) &&
                  ((vri->qualifier() == "promote") ||
                  (vri->qualifier() == "generalize"))) {
                visit -= vri->ident() ;
              }
            }
            breadth += working ;
            visit -= breadth ;
            working = visit ;
          }

          if((extract_rules(breadth) & outloop) == EMPTY) {
            if(debugging)
              debugout << "removing non-participating output rule " << *ri << endl ;
            remove_rules += *ri ;
          }
        }
          
      }

      for(ruleSet::const_iterator ri=remove_rules.begin();
          ri!=remove_rules.end(); ++ri)
        gr.remove_vertex(ri->ident()) ;

      all_cleaned_rules += remove_rules ;

    }
    
    do { // Keep cleaning until nothing left to clean!
        
      // Remove unnecessary vertexes from graph.
      int virtual_vertex = gr.max_vertex() + 1 ;
      digraph::vertexSet allvertices = gr.get_all_vertices() ;
      variableSet allvars = extract_vars(allvertices) ;
        
      //  target += variable(expression::create("OUTPUT")) ;
        
      gr.add_edges(virtual_vertex, given) ;
      gr.add_edges(target,virtual_vertex) ;
        
      const vector<digraph::vertexSet> components =
        component_sort(gr).get_components() ;
        
      digraph::vertexSet subset = EMPTY ;
        
      for(size_t i=0;i<components.size();++i)
        if(components[i].inSet(virtual_vertex)) {
          subset = components[i] ;
          break ;
        }
        
      subset -= virtual_vertex ;
      ruleSet rules = extract_rules(subset) ;
      ruleSet::const_iterator fi ;
      for(fi=rules.begin();fi!=rules.end();++fi)
        subset += fi->targets() ;

      if(verbose) {
        // some debugout info
        ruleSet rulesNOTinComponent =
          ruleSet(extract_rules(gr.get_all_vertices()) - rules) ;
        debugout << "The following rules are cleaned out because they"
                 << " are not in the targets->sources component: {{{{{"
                 << endl ;
        debugout << rulesNOTinComponent ;
        debugout << "}}}}}" << endl << endl ;
      }
        
      // Check for looping rules here, don't clean looping rule if it is
      // in the subset.
      digraph grt = gr.transpose() ;
      digraph::vertexSet  cleanout ;
      ruleSet tmp = rules ;
      tmp -= super_rules ;
      for(fi=tmp.begin();fi!=tmp.end();++fi) {
        if(fi->get_info().qualifier() != "looping")
          if((subset & fi->sources()) != fi->sources()) {
            if(debugging) {
              debugout << "cleanout " << *fi << endl ;
              debugout << "because of variables "
                       << extract_vars(fi->sources()-subset)
                       << endl ;
            }
            cleanout += fi->ident() ;
          }
      }

      ruleSet cleanoutrules = extract_rules(cleanout) ;
      subset -= cleanout ;
        
      variableSet touched_variables = given ;
      ruleSet working_rules = extract_rules(subset) ;
      for(ruleSet::const_iterator ri = working_rules.begin();
          ri!=working_rules.end();
          ++ri) {
        touched_variables += ri->targets() ;
      }
      
      ruleSet looping_rules ;
        
      digraph::vertexSet cleanout2 ;
      // don't clean out super rules
      tmp = working_rules ;
      tmp -= super_rules ;
      for(ruleSet::const_iterator ri = tmp.begin(); ri!=tmp.end();++ri) {
        if(ri->get_info().qualifier() != "looping") {
          if((ri->sources() - touched_variables)!=EMPTY) {
            cleanout2 += ri->ident() ;
            if(debugging) {
              debugout << "cleanout " << *ri << endl ;
              debugout << "because of variables "
                       << extract_vars(ri->sources()-touched_variables)
                       << endl ;
            }
          }
        } else
          looping_rules += *ri ;
      }

      subset -= cleanout2 ;
        
      cleanoutrules += extract_rules(cleanout2) ;
        
      WARN(subset == EMPTY) ;

      gr = gr.subgraph(subset) ;
        
      if(looping_rules != EMPTY) {
        digraph grt = gr.transpose() ;
        for(ruleSet::const_iterator ri = looping_rules.begin();
            ri!=looping_rules.end();
            ++ri) {
          variableSet sources = ri->sources() ;
          variableSet::const_iterator vi ;
          variableSet unused_vars ;
          for(vi=sources.begin();vi!=sources.end();++vi) {
            ruleSet rs = extract_rules(grt.get_edges(vi->ident())) ;
            if(rs == EMPTY) {
              unused_vars += *vi ;
            }
          }
          variableSet targets = ri->targets() ;
          for(vi=targets.begin();vi!=targets.end();++vi) {
            ruleSet rs = extract_rules(gr.get_edges(vi->ident())) ;
            if(rs == EMPTY) {
              unused_vars += *vi ;
            }
          }
          variableSet shared_vars = variableSet(sources & targets) ;
          unused_vars -= shared_vars ;
          variableSet newtargets = variableSet(targets-unused_vars) ;
          variableSet newsources = variableSet(sources-unused_vars) ;
          for(vi=newtargets.begin();vi!=newtargets.end();++vi) {
            variable tvar = *vi ;
            if(tvar.get_info().name != "OUTPUT" && !tvar.get_info().tvar) {
              // If a variable isn't being advanced in time, then
              // it has no business in the time loop 
              while(newtargets.inSet(tvar)) {
                tvar = tvar.new_offset(tvar.get_info().offset + 1) ;
              }
              if(!newsources.inSet(tvar)) {
                unused_vars += *vi ;
              }
            }
              
          }
            
          if(unused_vars != EMPTY) {
            variableSet looping_input = variableSet(sources-unused_vars) ;
            variableSet looping_output = variableSet(targets-unused_vars) ;
            ostringstream oss ;
            oss << "source("<< looping_input
                << "),target(" << looping_output
                << "),qualifier(looping)" ;
            rule floop(oss.str()) ;
              
            invoke_rule(floop,gr) ;
            gr.remove_vertex(ri->ident()) ;
            if(debugging) {
              debugout << "restructure iteration: " << floop << endl ;
              debugout << "originally was: " << *ri << endl ;
            }
          }
        }
      }

      gr.remove_dangling_vertices() ;

      all_cleaned_rules += cleanoutrules ;
      cleaned_rules = (cleanoutrules != EMPTY) ;

    } while (cleaned_rules) ;
  }
    
  // return only the rules that can be used with the given facts
  ruleSet active_rules(ruleSet rin,variableSet given) {
    ruleSet outrules ;
    digraph gr ;
    for(ruleSet::const_iterator ri=rin.begin();ri!=rin.end();++ri) {
      variableSet inputs = ri->constraints() ;
      if(inputs == EMPTY)
        inputs = ri->sources() ;
      variableSet outputs = ri->targets() ;
      int rid = ri->ident() ;
      inputs -= outputs ;
      variableSet::const_iterator vi ;
      variableSet vinSet,voutSet ;
      for(vi=inputs.begin();vi!=inputs.end();++vi) {
        variable v = *vi ;
        variable::info vinfo = v.get_info() ;
        vinfo.assign = false ;
        vinfo.time_id = time_ident() ;
        vinfo.offset = 0 ;
        vinfo.priority = vector<std::string>() ;
        variable vn(vinfo) ;
        if(vinfo.tvar)
          given += vn ;
        vinSet += vn ;
      }
      for(vi=outputs.begin();vi!=outputs.end();++vi) {
        variable v = *vi ;
        variable::info vinfo = v.get_info() ;
        vinfo.assign = false ;
        vinfo.time_id = time_ident() ;
        vinfo.offset = 0 ;
        vinfo.priority = vector<std::string>() ;
        variable vn(vinfo) ;
        voutSet += vn ;
      }
      vinSet -= voutSet ;
      if(vinSet == EMPTY) {
        given += voutSet ;
        outrules += *ri ;
      }
      for(vi=vinSet.begin();vi!=vinSet.end();++vi) {
        variable vn = *vi ;
        gr.add_edge(vn.ident(),rid) ;
      }
      for(vi=voutSet.begin();vi!=voutSet.end();++vi) {
        variable vn = *vi ;
        gr.add_edge(rid,vn.ident()) ;
      }

    }
    variableSet working = given ;
    variableSet visited_vars = working ;
    ruleSet visited_rules ;

    ruleSet cmp = extract_rules(gr.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=cmp.begin();ri!=cmp.end();++ri) {
        int id = ri->ident() ;
	
	if(ri->type() != rule::INTERNAL 
	   && ri->get_rule_implP()->get_rule_class() == rule_impl::SUPER_RULE) {
	  debugout << "super_rule " << *ri << endl ;
	  visited_rules += *ri ;
          visited_vars += extract_vars(gr[id]) ;
	}
    }
    digraph gt = gr.transpose() ;
    while(working != EMPTY) {
      // While we have vertexes to work on, compute additional vertexes that
      // can be scheduled
      ruleSet rule_consider ;
      variableSet::const_iterator ni ;
      // loop over working set and create a list of candidate vertexes
      for(ni=working.begin();ni != working.end(); ++ni) 
        rule_consider += extract_rules(gr[ni->ident()]) ;
        
      rule_consider -= visited_rules ;
      ruleSet::const_iterator ri ;
      ruleSet new_rules ;
      variableSet new_vars ;
      for(ri=rule_consider.begin();ri!=rule_consider.end();++ri) {
        int id = ri->ident() ;
	
	if(ri->type() != rule::INTERNAL 
	   && ri->get_rule_implP()->get_rule_class() == rule_impl::SUPER_RULE) {
	  new_rules += *ri ;
          new_vars += extract_vars(gr[id]) ;
	} else if((entitySet(extract_vars(gt[id])) & entitySet(visited_vars))
           == entitySet(gt[id])) {
          new_rules += *ri ;
          new_vars += extract_vars(gr[id]) ;
        }
      }
      new_vars -= visited_vars ;
      visited_vars += new_vars ;
      working = new_vars ;
      visited_rules += new_rules ;
    }
    

    outrules += visited_rules ;

    if(MPI_processes == 1 || verbose) {
      rin -= outrules ;
      if(rin != EMPTY) {
        debugout << "precleaning rules that cannot be scheduled based on given facts:" << endl ;
        for(ruleSet::const_iterator ri=rin.begin();ri!=rin.end();++ri) {
          variableSet x = extract_vars(gt[ri->ident()]) ;
          x -= visited_vars ;
          debugout << "eliminating " << *ri  
                   << " due to " << x << endl ;
        }
      }
    }
    return outrules ;
  }

  std::ostream &output_graph(const digraph &dg,std::ostream &s) {
    s << "graph = {" << endl ;
    s << "nodes = {" << endl ;
    variableSet vars = extract_vars(dg.get_all_vertices()) ;
    ruleSet rules = extract_rules(dg.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      s << ri->ident() << ": " << '"' ;
      if(ri->type() != rule::INTERNAL) {
        rule_implP rp = ri->get_rule_implP() ;
        s<< '[' << rp->get_name() << ']';
      }
      s << *ri <<'"' << endl ;
    }
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      s << vi->ident() << ": " << '"' << *vi << '"' << endl ;
    }
    s << "}" << endl ;
    s << "connectivity = {"  << endl ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      int id = ri->ident() ;
      digraph::vertexSet conn = dg[id] ;
      digraph::vertexSet::const_iterator ii ;
      s << id <<':' ;
      for(ii=conn.begin();ii!=conn.end();++ii)
        s << ' ' << *ii ;
      s << " ;" << endl ;
    }
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      int id = vi->ident() ;
      digraph::vertexSet conn = dg[id] ;
      digraph::vertexSet::const_iterator ii ;
      s << id <<':' ;
      for(ii=conn.begin();ii!=conn.end();++ii)
        s << ' ' << *ii ;
      s << " ;" << endl ;
    }
    s << "}" << endl ;
    s << "}" << endl ;
    return s ;
  }
  

  variable flattenVariable(variable v) {
    variable::info vinfo = v.get_info() ;
    vinfo.assign = false ;
    vinfo.time_id = time_ident() ;
    vinfo.offset = 0 ;
    vinfo.priority = vector<std::string>() ;
    variable vn(vinfo) ;
    return vn ;
  }
  
  void makeFlatGraph(digraph &gr,const rule_db &rdb) {
    ruleSet all_rules = rdb.all_rules() ;

    ruleSet::const_iterator ri ;
    variableSet all_vars ;
    for(ri=all_rules.begin();ri!=all_rules.end();++ri) {
      int rid = ri->ident() ;
      
      variableSet vin = ri->sources() ;

      variableSet::const_iterator vi ;
      for(vi=vin.begin();vi!=vin.end();++vi) {
        variable v = flattenVariable(*vi) ;
        gr.add_edge(v.ident(),rid) ;
      }
      variableSet vout = ri->targets() ;

      for(vi=vout.begin();vi!=vout.end();++vi) {
        variable v = flattenVariable(*vi) ;
        gr.add_edge(rid,v.ident()) ;
      }
    }
    
  }

  //#define CLEANREPORT
  //#define REPORT
  dependency_graph2::dependency_graph2(const rule_db& rdb,
                                       const variableSet& given,
                                       const variableSet& target) {
    // at first, we get all the rules in the rule database
    //ruleSet all_rules = rdb.all_rules() ;
    ruleSet all_rules = active_rules(rdb.all_rules(),given) ;

#ifdef REPORT
    char *p ;
    char buf[2048] ;
    
    if(getcwd(buf,sizeof(buf))==0) 
      p = "PWD" ;
    else
      p=buf ;
    string prefix = "/usr/tmp/";
    char *t ;
    if((t = rindex(p,'/')) == 0)
      t = p ;
    else {
      *t = '_' ;
      char *t2 ;
      if((t2 = rindex(p,'/')) !=0) {
        t = t2 ;
        *t = '_' ;
      }
      if((t2 = rindex(p, '/')) !=0)
        t = t2 ;
      t++ ;
    }
    prefix += string(t) ;
    {
      string filename = prefix+string("rules.dat") ;
      
      std::ofstream file(filename.c_str(),std::ios::out) ;
      file << "rules = {" << endl ;
      for(ruleSet::const_iterator ri=all_rules.begin();
          ri!=all_rules.end();
          ++ri) {
        file << '"' ;
        if(ri->type() != rule::INTERNAL) {
          rule_implP rp = ri->get_rule_implP() ;
          file << '[' << rp->get_name() << ']';
        }
        file << *ri <<'"' << endl ;
      }
      file << "}" << endl ;
    }
    {
      string filename = prefix+string("globalgraph.dat") ;
      std::ofstream file(filename.c_str(),std::ios::out) ;
      digraph grg ;
      makeFlatGraph(grg,rdb) ;
      output_graph(grg,file) ;
    }
#endif

    // then we classify all the iterations in the rule database,
    // while also pick out non iteration rules
    // all the iteration information is stored in an
    // iteration_info object, other rules are in working_rules
    iteration_info iter ;
    ruleSet working_rules ;
    // then we call the function to classify rules
    classify_ruledb(all_rules,iter,working_rules) ;
    // next, we create a representative rule for
    // each iteration
    create_iteration_rep(iter,working_rules) ;
    // followed, we will need to check the iteration rules 
    // if check failed, we just return an empty graph
    if(!check_iteration(iter))
      return ; // because gr is empty now

    // finally, we will need to remove any rule that
    // generates any given variable. because given
    // variables are assumed to be extensional and
    // do not need to be computed again.
    {
      ruleSet rules_to_remove ;
      for(ruleSet::const_iterator ri=working_rules.begin();
          ri!=working_rules.end();++ri) {
        variableSet targets = ri->targets() ;
        // we only remove rule whose targets are COMPLETELY
        // included in the given variableSet
        if( (targets - given) == EMPTY)
          rules_to_remove += *ri ;
      }
      working_rules -= rules_to_remove ;
    }


    // then we build a digraph that has all the working
    // rules (stationary rules + time specific rules +
    // iterating rules) inside
    digraph rule_graph ;
    for(ruleSet::const_iterator ri = working_rules.begin();
        ri != working_rules.end();
        ++ri) {
      invoke_rule_wp(*ri,rule_graph) ;
    }

    //cerr<<"rule_graph size: "<<rule_graph.get_all_vertices().size()<<endl ;
    // Some checking to the graph
    {
      variableSet rg_allvars = extract_vars(rule_graph.get_all_vertices()) ;
      variableSet target_diff = variableSet(target-rg_allvars) ;
      if(target_diff != EMPTY) {
        cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
        cerr << "ERROR: insufficient rule database, "
             << " not all requested targets can be computed." << endl ;
        cerr << "\tVariable(s): " << target_diff
             << " cannot be inferred from the given rule database."
             << endl ;
        // gr is still empty, we return an empty graph
        // which means the query is not satisfied
	return ;
      }
    }
    // we are now ready to build the graph
    digraph rule_graph_transpose = rule_graph.transpose() ;
    // we collect the toplevel requests
    variableSet top_request = 
      create_graph(rule_graph,rule_graph_transpose,
                   target,iter,gr,time_ident()) ;
    // we need to compare the top_request with the given
    // variables to see if the schedule is possible
//     if( (top_request - given) != EMPTY) {
//       cerr << __FILE__ << ", Line " << __LINE__ << ": "  ;
//       cerr << "ERROR: insufficient fact database!" << endl ;
//       cerr << "These facts are required to compute the query: " << endl ;
//       cerr << "        " << variableSet(top_request-given) << endl ;
//       gr = digraph() ;
//       return ;
//     }
    
    // we now add these built iteration graphs
    for(map<time_ident,iteration>::iterator ip=iter.iteration_rules.begin();
        ip!=iter.iteration_rules.end();++ip)
      if(ip->second.active) {
        gr.add_graph(ip->second.iteration_graph) ;
      }

#ifndef CLEANREPORT
#ifdef REPORT
    {
      string filename = prefix+string("graph.dat") ;
      std::ofstream file(filename.c_str(),std::ios::out) ;

      file << "input = {" << given << "}"<<endl ;
       
      output_graph(gr,file) ;
    }
#endif
#endif
    //cerr<<"vertices size before cleaning: "
    //  <<gr.get_all_vertices().size()<<endl ;
    clean_graph(gr,given,target) ;

    //cerr<<"vertices size after cleaning: "
    //  <<gr.get_all_vertices().size()<<endl ;


    // Partition any iterations that remain that have multiple collapse
    // conditionals into independent iterations
    gr = partition_iteration(gr) ;

#ifdef CLEANREPORT
#ifdef REPORT
    {
      string filename = prefix+string("graph.dat") ;
      std::ofstream file(filename.c_str(),std::ios::out) ;
      gr.remove_dangling_vertices() ;

      file << "input = {" << given << "}"<<endl ;
       
      output_graph(gr,file) ;
    }
#endif
#endif
    // Add dependencies created by update-in-place rules
    add_rename_dependencies(gr) ;

    // Finish cleaning things up
    gr.remove_dangling_vertices() ;
  }

  //Search through a graph from a given set of vertexes to all connected
  // vertexes. Return the set of all vertexes related through edges in
  // the graph to the start set provided.
  digraph::vertexSet search_digraph(const digraph &gr,
                                    digraph::vertexSet start) {
    digraph::vertexSet working = start ;
    digraph::vertexSet delta ;
    do {
      digraph::vertexSet newset ;
      digraph::vertexSet::const_iterator vi ;
      for(vi=working.begin();vi!=working.end();++vi)
        newset += gr.get_edges(*vi) ;
      newset += working ;
      delta = newset - working ;
      working = newset ;
    } while(delta != EMPTY) ;
    return working ;
  }

  // This function partitions iterations that have collapse rules conditional
  // on different variables into completely different iterations.  This allows
  // us to load modules that share the same iteration tag but are iterating on
  // completely different data and have different termination conditions.
  digraph partition_iteration(digraph gr) {

    // First we search through the collapse rules to see if there are any
    // collapse rules that conflict.
    ruleSet all_rules = extract_rules(gr.get_all_vertices()) ;
    ruleSet looping_rules ;
    ruleSet::const_iterator ri ;
    map<time_ident,ruleSet> tmap ;
    for(ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->type() == rule::COLLAPSE) {
        tmap[ri->source_time()] += *ri ;
      }
    }

    vector<time_ident> conflicts ;
    map<time_ident,ruleSet>::const_iterator mi ;
    for(mi=tmap.begin();mi!=tmap.end();++mi) {
      if(mi->second.size() > 1) {
        ri = mi->second.begin() ;
        variableSet cond = ri->get_info().desc.conditionals ;
        ri++;
        for(;ri!=mi->second.end();++ri) {
          if(cond != ri->get_info().desc.conditionals) {
            conflicts.push_back(mi->first) ;
          }
        } 
      }
    }
    // If there are no conflicts then we are done
    if(conflicts.size() == 0)
      return gr ;

    // Now for each set of conflicting collapse rules, we need to create an
    // independent iteration
    for(size_t i=conflicts.size();i>0;--i) {
      time_ident iter = conflicts[i-1] ;
      ruleSet col = tmap[iter] ;
      // First we collect information about this iteration level
      ruleSet all_rules = extract_rules(gr.get_all_vertices()) ;
      variableSet all_vars = extract_vars(gr.get_all_vertices()) ;
      variableSet::const_iterator vi;
      ruleSet iter_rules ;
      variableSet iter_vars ;
      ruleSet looping_rule ;
      variableSet promote_vars ;
      for(ri=all_rules.begin();ri!=all_rules.end();++ri) {
        if(ri->target_time() == iter || (iter < ri->target_time())) {
          iter_rules += *ri ;
          if(ri->target_time() == iter && ri->get_info().qualifier() == "looping")
            looping_rule += *ri ;
          if(ri->get_info().qualifier() == "promote") {
            promote_vars += ri->targets() ;
            promote_vars += ri->sources() ;
          }
          iter_vars += ri->targets() ;
          iter_vars += ri->sources() ;
        }
      }
      for(vi=all_vars.begin();vi!=all_vars.end();++vi)
        if(vi->time() == iter || iter < vi->time())
          iter_vars += *vi ;

      if(looping_rule.size() == 0)
        continue ;
      
      // There should only be one looping rule at this stage, if not
      // something weird is going on.
      if(looping_rule.size() != 1) {
        cerr << "something confused in iteration conflict analysis" << endl;
        cerr << "looping_rule = " << looping_rule << endl ;
        return gr ;
      }

      iter_rules -= looping_rule ;

      // Now we need to determine what subset of the graph contributes to each
      // independent loop
      digraph::vertexSet ss = iter_rules ;
      ss += iter_vars ;
      digraph sg = gr.subgraph(ss) ;
      sg.remove_dangling_vertices();
      digraph sgt = sg.transpose() ;

      map<variableSet,ruleSet> collapse_groups ;
      for(ri=col.begin();ri!=col.end();++ri) {
        collapse_groups[ri->get_info().desc.conditionals] += *ri ;
      }
      map<variableSet,ruleSet>::const_iterator cmi ;
      vector<digraph::vertexSet> groups ;

      for(cmi=collapse_groups.begin();cmi!=collapse_groups.end();++cmi) {
        digraph::vertexSet v1 ;
        for(ri=cmi->second.begin();ri!=cmi->second.end();++ri) {
          v1 = search_digraph(sgt,ri->sources()) ;
          v1 += ri->ident() ;
        }
        variableSet v1var = extract_vars(v1) ;
        v1var &= looping_rule.begin()->targets() ;
        
        for(vi=v1var.begin();vi!=v1var.end();++vi) {
          variable::info vinfo = vi->get_info() ;
          if(vinfo.tvar || vinfo.name == "OUTPUT")
            continue ;
          vinfo.offset++ ;
          variable fv = variable(vinfo) ;
          v1 += fv.ident() ;
        }
        v1 = search_digraph(sgt,v1) ;
        groups.push_back(v1) ;
      }


      // Here we check to make sure the loops really are independent of each
      // other.

      digraph::vertexSet tot ;
      for(size_t i=0;i<groups.size();++i) {
        if((tot & groups[i]) != EMPTY) {
          variableSet check = extract_vars(tot&groups[i]) ;
          check -= promote_vars ;
          check -= variable(iter) ;
          if(check != EMPTY) {
            cerr << "Warning:" << endl 
                 << "  iteration will yield duplicate computations when collapse"
                 << endl
                 << "  are split to accommodate different collapse conditionals"
                 << endl
                 << " Duplicated Variables are " << check << endl
                 << " Collapse rules are " << col << endl ;
          }
        }
        tot += groups[i] ;
      }

      // Now we remove the old set of rules from the graph and replace them
      // with rules that have their variables renamed such that they are
      // independent loops
      digraph newgr = gr ;
      for(size_t i=0;i<groups.size();++i) {
        newgr.remove_vertices(extract_rules(groups[i])) ;
      }
      newgr.remove_vertices(looping_rule) ;
      
      for(size_t i=0;i<groups.size();++i) {
        string level_name = iter.level_name() ;
        char buf[512] ;
        int ii = i ;
        snprintf(buf,512,"_%s_%d_",level_name.c_str(),ii) ;
        string new_level = string(buf) ;
        time_ident ntime = time_ident(iter.parent(),new_level) ;
        map<variable,variable> rvm ;
        variable tv1 = variable(iter) ;
        variable tv2 = variable(ntime) ;
        rvm[tv1] = tv2 ;
        variableSet vars = extract_vars(groups[i]) ;
        for(vi=vars.begin();vi!=vars.end();++vi) {
          variable::info vinfo = vi->get_info() ;
          if(vinfo.time_id < iter) {
            rvm[*vi] = *vi ;
            continue ;
          }
          if(vinfo.time_id != iter) {
            cerr << "time_id = " << vinfo.time_id << endl ;
            cerr << "current code unable to decompose nested iterations!"
                 << endl ;
            return gr ;
          }
          if(vinfo.tvar)
            continue ;
          vinfo.time_id = ntime ;
          rvm[*vi] = variable(vinfo) ;
        }

        ruleSet grules = extract_rules(groups[i]) ;
        for(ri=grules.begin();ri!=grules.end();++ri) {
          rule nr = (*ri).rename_vars(rvm) ;
          invoke_rule(nr,newgr) ;
        }
        variableSet loopin = looping_rule.begin()->sources() ;
        variableSet loopout = looping_rule.begin()->targets() ;
        loopin &= extract_vars(groups[i]) ;
        loopout &= extract_vars(groups[i]) ;
        ostringstream oss ;
        oss << "source("<< loopin
            << "),target(" << loopout
            << "),qualifier(looping)" ;
        rule floop(oss.str()) ;
        rule nr = floop.rename_vars(rvm) ;
        invoke_rule(nr,newgr) ;

      }
      gr = newgr ;
    }
    return gr ;
  }

}// End namespace Loci
