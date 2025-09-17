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
#include <store.h>
#include "comp_tools.h"
#include "dist_tools.h"
#include "visitorabs.h"
#include "thread.h"
#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <list>
using std::list ;
#include <set>
using std::set ;

#include "loci_globs.h"
//#define VERBOSE

namespace Loci {
  extern bool in_internal_query;
  extern bool threading_recursion;
  extern int num_total_recursion;
  extern int num_threaded_recursion;

  void impl_recurse_compiler::accept(visitor& v) {
    v.visit(*this) ;
  }

  void impl_recurse_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    
#ifdef VERBOSE
    debugout << "set var existence for recursive impl rule " << impl << endl ;
#endif

    variableSet::const_iterator vi ;
    variableSet tvars ;
    tvars = impl.targets() ;
    fcontrol &fctrl = control_set ;
    fctrl.recursion_maps.clear();
    fastseq = EMPTY;
    const rule_impl::info &finfo = impl.get_info().desc ;
    warn(impl.type() == rule::INTERNAL) ;

    entitySet my_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      variableSet recurse_vars = variableSet(impl.sources() & impl.targets()) ;

      barrier_existential_rule_analysis(recurse_vars,facts, scheds) ;

      my_entities = d->my_entities ;
    }
    entitySet sources = ~EMPTY ;
    entitySet constraints = ~EMPTY ;

    set<vmap_info>::const_iterator si ;
    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
      if((si->var & tvars) == EMPTY)
        sources &= vmap_source_exist(*si,facts, scheds) ;
      else {
        const int num_maps = si->mapping.size() ;
        if(num_maps != 0) {
          for(vi=si->mapping[0].begin();vi!=si->mapping[0].end();++vi) {
            sources &= scheds.variable_existence(*vi) ;
          }
        }

        tvars += si->var ;
        vector<variableSet::const_iterator> miv(num_maps) ;
        for(int i=0;i<num_maps;++i)
          miv[i] = si->mapping[i].begin() ;

        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          int i = 0 ;
          do {
            fcontrol::mapping_info minfo ;
            for(int j=0;j!=num_maps;++j) {
              FATAL(!scheds.is_a_Map(*(miv[j]))) ;
              MapRepP m = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
              minfo.mapvec.push_back(m) ;
              minfo.mapvar.push_back(*(miv[j])) ;
            }
            minfo.v = *vi ;

            for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
              miv[i] = si->mapping[i].begin() ;
            }
            fctrl.recursion_maps.push_back(minfo) ;
          } while(i!=num_maps) ;
        }
      }
    }
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      int num_maps = si->mapping.size() ;
      vector<variableSet::const_iterator> miv(num_maps) ;
      for(int i=0;i<num_maps;++i)
        miv[i] = si->mapping[i].begin() ;

      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        int i = 0 ;
        do {
          fcontrol::mapping_info minfo ;
          for(int j=0;j!=num_maps;++j) {
            FATAL(!scheds.is_a_Map(*(miv[j]))) ;
            MapRepP m = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
            minfo.mapvec.push_back(m) ;
            minfo.mapvar.push_back(*(miv[j])) ;
          }
          minfo.v = *vi ;

          for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
            miv[i] = si->mapping[i].begin() ;
          }
          fctrl.target_maps.push_back(minfo) ;
        } while(i!=num_maps) ;
      }
    }
    fctrl.use_constraints = false ;
    for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si) {
      warn((si->var & tvars) != EMPTY) ;
      constraints &= vmap_source_exist(*si,facts, scheds) ;
      fctrl.use_constraints = true ;
    }
    entitySet comp_sources;
    if(duplicate_work)
      comp_sources = sources;

    sources += fill_entitySet(sources,facts) ;
    if(fctrl.use_constraints)
      constraints += fill_entitySet(constraints,facts) ;

    fctrl.nr_sources = sources ;

    fctrl.constraints = constraints ;
    warn(fctrl.recursion_maps.size() == 0) ;

    fatal(tvars.size() != 1 ||
          control_set.recursion_maps.size() != 1 ||
          control_set.target_maps.size() != 1 ||
          fctrl.use_constraints) ;

    fcontrol::mapping_info &rmap = fctrl.recursion_maps[0] ;
    fcontrol::mapping_info &tmap = fctrl.target_maps[0] ;
    vector<multiMap> read_maps,read_map_inv ;

    for(size_t j=0;j<rmap.mapvec.size();++j) {
      read_maps.push_back(rmap.mapvec[j]->get_map()) ;
    }
    for(size_t j=0;j<read_maps.size();++j) {
      read_map_inv.push_back(multiMap()) ;
    }
    variable rvar = *(tvars.begin()) ;


    entitySet sdelta = scheds.variable_existence(rvar) ;

    entitySet initial = sdelta ;
    entitySet domain = fctrl.nr_sources + sdelta ;

    if(facts.isDistributed()) {
      variable rename_var = rvar ;
      if(finfo.targets.begin()->assign.size() != 0) {
        rename_var = finfo.targets.begin()->assign[0].second ;
      }
      entitySet start=scheds.variable_existence(rename_var) ;
      domain += start - my_entities ;
      sdelta += start - my_entities ;
    }


    for(int j=read_maps.size()-1;j>=0;--j) {
      entitySet newdomain = scheds.preimage(rmap.mapvar[j],domain).first ;
#ifdef VERBOSE
      debugout << "j = " << j << ", domain = " << domain
               << ", newdomain = " << newdomain << endl ;
#endif
      if(domain == ~EMPTY) {
        if(MPI_processes == 1)
          cerr << "problem in recurse compiler for rule = "<< impl << endl ;
        else
          debugout << "problem in recurse compiler for rule = "<< impl << endl ;
        domain = read_maps[j].domain() ;
        scheds.set_error() ;
      }
      inverseMap(read_map_inv[j],read_maps[j],domain,newdomain) ;
      domain = newdomain ;
    }

    for(int j=rmap.mapvec.size()-1;j>=0;--j)
      sdelta = rmap.mapvec[j]->preimage(sdelta).first ;

    entitySet comp_generated;
    if(duplicate_work) {
      entitySet comp_sdelta = sdelta;
      comp_sdelta &= comp_sources;
      comp_sdelta &= my_entities;
      entitySet comp_tdelta = comp_sdelta;
      for(size_t j=0;j<tmap.mapvec.size();++j)
	comp_tdelta = tmap.mapvec[j]->image(comp_tdelta) ;
      comp_generated = comp_tdelta;
    }
    sdelta &= fctrl.nr_sources ;
    sdelta &= my_entities ;

    entitySet tdelta = sdelta ;
    for(size_t j=0;j<tmap.mapvec.size();++j)
      tdelta = tmap.mapvec[j]->image(tdelta) ;

#ifdef VERBOSE
    debugout << "sdelta_init = " << sdelta << ", tdelta = " << tdelta << endl ;
#endif

    fastseq += sequence(sdelta) ;

    entitySet generated = tdelta ;
    const entitySet nr_sources = fctrl.nr_sources ;
    entitySet exists_alloc = nr_sources ;
    if(facts.isDistributed()) {
      entitySet working = exists_alloc ;
      for(size_t j=0;j<read_maps.size();++j) {
        MapRepP mp = MapRepP(read_maps[j].Rep()) ;
        working = mp->image(working) ;
      }
      exists_alloc += working ;
      working = exists_alloc ;
      for(size_t j=0;j<tmap.mapvec.size();++j) {
        working = tmap.mapvec[j]->image(working) ;
      }
      exists_alloc += working ;
    }
    exists_alloc += tdelta ;
    exists_alloc += initial ;
    store<bool> exists ;
    exists.allocate(exists_alloc) ;
#ifdef VERBOSE
    debugout << "exists_alloc = " << exists_alloc
             << "nr_sources = " << nr_sources << endl ;
#endif
    for(entitySet::const_iterator
          ei=exists_alloc.begin();ei!=exists_alloc.end();++ei)
      exists[*ei] = false ;

    for(entitySet::const_iterator
          ei=initial.begin();ei!=initial.end();++ei)
      exists[*ei] = true ;

    exists_alloc -= my_entities ;
#ifdef VERBOSE
    debugout << "exists_alloc-my_entities = " << exists_alloc << endl ;
    debugout << "my_entities = " << my_entities << endl ;
#endif
    for(entitySet::const_iterator
          ei=exists_alloc.begin();ei!=exists_alloc.end();++ei)
      exists[*ei] = true ;

    for(entitySet::const_iterator
          ei=tdelta.begin();ei!=tdelta.end();++ei)
      exists[*ei] = true ;

    do {
      for(int j=read_map_inv.size()-1;j>=0;--j) {
        entitySet candidates ;
        const int *mi ;
        if((sdelta - read_map_inv[j].domain()) != EMPTY) {
          cerr << "problem in processing recursive rule " << impl
               << endl ;
          Loci::Abort() ;
          break ;
        }

        for(entitySet::const_iterator di=sdelta.begin();di!=sdelta.end();++di) {
#ifdef DEBUG
          if(!read_map_inv[j].domain().inSet(*di)) {
            cerr << " di = " << *di << ", j = " << j << endl ;
          } else
#endif
            //        fatal(!read_map_inv[j].domain().inSet(*di)) ;
            for(mi=read_map_inv[j].begin(*di); mi!=read_map_inv[j].end(*di);++mi)
              candidates += *mi ;
        }
        sdelta = candidates ;
      }
#ifdef VERBOSE
      debugout << "candidates = " << sdelta << endl ;
      debugout << "nr_sources = " << nr_sources << endl ;
#endif
      // we now have candidate sdeltas, check them to see if they
      // are satisfied.
      entitySet satisfied ;
      vector<const int *> mi(read_maps.size()) ;
      vector<const int *> me(read_maps.size()) ;
      for(entitySet::const_iterator di=sdelta.begin();di!=sdelta.end();++di) {
        int c = *di ;
        mi[0] = read_maps[0].begin(c) ;
        me[0] = read_maps[0].end(c) ;
        int j=0 ;
        const int je=read_maps.size();
        bool chk = nr_sources.inSet(c) ;
        // Check that all paths from candidate go to generated cells
        while(chk && j>=0) {
          c = *mi[j] ;
          j++ ;
          if(j==je) {
            chk = chk && exists[c] ; //&& nr_sources.inSet(c) ;
            do {
              j-- ;
            } while(j>=0 && (++mi[j]) == me[j]) ;
          } else {
            mi[j] = read_maps[j].begin(c) ;
            me[j] = read_maps[j].end(c) ;
          }

        }
        if(chk)
          satisfied += *di ;
        if(!chk) {
          int c = *di ;
          mi[0] = read_maps[0].begin(c) ;
          me[0] = read_maps[0].end(c) ;
          chk = nr_sources.inSet(c) ;
          int j=0 ;
          const int je=read_maps.size();
          while(chk && j>=0) {
            c = *mi[j] ;
            j++ ;
            if(j==je) {
              chk = chk && exists[c] ;// && nr_sources.inSet(c) ;
              do {
                j-- ;
              } while(j>=0 && (++mi[j]) == me[j]) ;
            } else {
              mi[j] = read_maps[j].begin(c) ;
              me[j] = read_maps[j].end(c) ;
            }
          }
        }
      }
      sdelta = satisfied ;
      sdelta &= my_entities ;
      entitySet tdelta = sdelta ;
      for(size_t j=0;j<tmap.mapvec.size();++j)
        tdelta = tmap.mapvec[j]->image(tdelta) ;

#ifdef VERBOSE
      debugout << "sdelta = " << sdelta << ", tdelta = " << tdelta << endl ;
#endif
      fastseq += sequence(sdelta) ;

      for(entitySet::const_iterator
            ei=tdelta.begin();ei!=tdelta.end();++ei)
        exists[*ei] = true ;

    } while(sdelta != EMPTY) ;

    for(entitySet::const_iterator
          ei=nr_sources.begin();ei!=nr_sources.end();++ei)
      if(exists[*ei])
        generated += *ei ;

    if(duplicate_work) {
      comp_sources &= my_entities;
      for(entitySet::const_iterator
	    ei=comp_sources.begin();ei!=comp_sources.end();++ei) {
	if(exists[*ei])
	  comp_generated += *ei ;
      }
    }

    fctrl.generated[rvar] = generated ;
#ifdef VERBOSE
    debugout << "recursive rule " << impl << " generating " << generated << endl ;
#endif

    for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
        mi!=fctrl.generated.end();++mi) {
      scheds.set_existential_info(mi->first,impl,mi->second) ;
    }

    entitySet create = scheds.get_existential_info(rvar,impl) ;
    //Add duplication information to recursive variable
    if(duplicate_work) {
      scheds.set_my_proc_able_entities(rvar, impl, comp_generated);
      scheds.add_policy(rvar, sched_db::NEVER);
    }
    create += send_entitySet(create,facts) ;
    create += fill_entitySet(create,facts) ;
    scheds.set_existential_info(rvar,impl,create) ;
   
  }

  void impl_recurse_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    entitySet my_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      my_entities = d->my_entities ;
    }
    process_rule_requests(impl,facts, scheds) ;
    if(facts.isDistributed()) {

      // For the relaxed recursion we need to adjust our variable request
      variableSet tvars = impl.targets() ;
      variable rvar = *(tvars.begin()) ;
      variable rename_var = rvar ;
      const rule_impl::info &finfo = impl.get_info().desc ;
      if(finfo.targets.begin()->assign.size() != 0) {
        rename_var = finfo.targets.begin()->assign[0].second ;
      }
      entitySet request = scheds.get_variable_requests(rvar) ;
      request -= my_entities ;
      scheds.variable_request(rename_var,request) ;

      //Find duplication of variables that are associtated with
      //rules that compute tvars
      if(duplicate_work)
	set_duplication_of_variables(tvars, scheds, facts);

      list<comm_info> request_comm ;
      variableSet recurse_vars = variableSet(impl.sources() & impl.targets()) ;
      request_comm = barrier_process_rule_requests(recurse_vars, facts, scheds) ;
      // list<comm_info> clist = sort_comm(request_comm,facts) ;
      // scheds.update_comm_info_list(clist, sched_db::RECURSE_CLIST);
    }
     
  }

  executeP impl_recurse_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {

    executeP exe_rule = new execute_rule(impl, fastseq, facts, scheds);
    return exe_rule;
  }

  void recurse_compiler::accept(visitor& v) {
    v.visit(*this) ;
  }

  void recurse_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    entitySet my_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      
      std::vector<std::pair<variable,entitySet> > pre_send_entities
        = barrier_existential_rule_analysis(recurse_vars,facts, scheds) ;
      scheds.update_send_entities(pre_send_entities, sched_db::RECURSE_PRE);
      my_entities = d->my_entities ;
      
    }

    control_set.clear() ;
    recurse_send_entities.clear();
    send_req_var.clear();
    ruleSet::const_iterator fi ;
    variableSet::const_iterator vi ;
    variableSet tvars ;
    const ruleSet &fset = recurse_rules ;
    for(fi=fset.begin();fi!=fset.end();++fi)
      tvars += fi->targets() ;
    for(fi=fset.begin();fi!=fset.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      const rule_impl::info &finfo = fi->get_info().desc ;
      warn(fi->type() == rule::INTERNAL) ;
      entitySet sources = my_entities ;
      entitySet constraints = my_entities ;
      set<vmap_info>::const_iterator si ;
      for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
        if((si->var & tvars) == EMPTY)
          sources &= vmap_source_exist(*si,facts, scheds) ;
        else {
          tvars += si->var ;
          int num_maps = si->mapping.size() ;
          vector<variableSet::const_iterator> miv(num_maps) ;
          for(int i=0;i<num_maps;++i)
            miv[i] = si->mapping[i].begin() ;

          for(vi=si->var.begin();vi!=si->var.end();++vi) {
            int i = 0 ;
            do {
              fcontrol::mapping_info minfo ;
              for(int j=0;j!=num_maps;++j) {
                FATAL(!scheds.is_a_Map(*(miv[j]))) ;
                MapRepP mp = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
                minfo.mapvec.push_back(mp) ;
              }
              minfo.v = *vi ;

              for(i=0;(i!=num_maps)&&(++(miv[i]) == si->mapping[i].end());++i)
                miv[i] = si->mapping[i].begin() ;

              fctrl.recursion_maps.push_back(minfo) ;
            } while(i!=num_maps) ;
          }
        }
      }
      for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
        int num_maps = si->mapping.size() ;
        vector<variableSet::const_iterator> miv(num_maps) ;
        for(int i=0;i<num_maps;++i)
          miv[i] = si->mapping[i].begin() ;

        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          int i = 0 ;
          do {
            fcontrol::mapping_info minfo ;
            for(int j=0;j!=num_maps;++j) {
              FATAL(!scheds.is_a_Map(*(miv[j]))) ;
              MapRepP mp = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
              minfo.mapvec.push_back(mp) ;
            }
            minfo.v = *vi ;

            for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
              miv[i] = si->mapping[i].begin() ;
            }
            fctrl.target_maps.push_back(minfo) ;
          } while(i!=num_maps) ;
        }
      }
      fctrl.use_constraints = false ;
      for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si) {
        warn((si->var & tvars) != EMPTY) ;
#ifdef DEBUG
	if((si->var & tvars) != EMPTY) {
	  cerr << "fset=" << fset << endl ;
	}
#endif
        constraints &= vmap_source_exist(*si,facts, scheds) ;
        fctrl.use_constraints = true ;
      }
      fctrl.nr_sources = sources ;
      fctrl.constraints = constraints ;
      warn(fctrl.recursion_maps.size() == 0) ;
    }


    map<variable,entitySet> tvar_computed,tvar_update ;
    for(vi=tvars.begin();vi!=tvars.end();++vi) {
      tvar_computed[*vi] = scheds.variable_existence(*vi) ;
      tvar_update[*vi] = tvar_computed[*vi] ;
    }
    map<variable,entitySet>::iterator mi ;
    for(mi=tvar_update.begin();mi!=tvar_update.end();++mi)
      mi->second = EMPTY ;

    bool finished = false;

    while(!finished) {
      for(fi=fset.begin();fi!=fset.end();++fi) {
        fcontrol &fctrl = control_set[*fi] ;
        // partial vmap_souce_exist
        entitySet srcs = fctrl.nr_sources ;
        for(size_t i=0;i<fctrl.recursion_maps.size();++i) {
          entitySet sdelta = tvar_computed[fctrl.recursion_maps[i].v] ;
          for(int j=int(fctrl.recursion_maps[i].mapvec.size())-1;j>=0;--j)
            sdelta = fctrl.recursion_maps[i].mapvec[j]->preimage(sdelta).first ;
          fctrl.recursion_maps[i].total += sdelta ;
          srcs &= fctrl.recursion_maps[i].total ;
        }
        if(fctrl.use_constraints) {
          if((srcs & fctrl.constraints) != fctrl.constraints) {
            if(MPI_processes == 1) {
              cerr << "recursive rule: " << *fi
                   << " cannot supply all entitites in constraint" << endl ;
              cerr << "constraints = " << fctrl.constraints << endl ;
              entitySet sac = srcs & fctrl.constraints ;
              cerr << "srcs & constraints = " << sac << endl ;
            } else {
              debugout << "recursive rule: " << *fi
                       << " cannot supply all entitites in constraint" << endl ;
              debugout << "constraints = " << fctrl.constraints << endl ;
              entitySet sac = srcs & fctrl.constraints ;
              debugout << "srcs & constraints = " << sac << endl ;
            }
            scheds.set_error() ;
            //exit(-1) ;

            const rule_impl::info &rinfo = fi->get_info().desc ;
            set<vmap_info>::const_iterator si ;



            if(MPI_processes == 1) {
              cerr << "detailed report:" << endl ;
              cerr<< "fctrl.nr_sources = " << fctrl.nr_sources << endl ;
            } else {
              debugout << "detailed report:" << endl ;
              debugout << "fctrl.nr_sources = " << fctrl.nr_sources << endl ;
            }
            entitySet constraints = fctrl.constraints ;
            for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
              entitySet sources = vmap_source_exist(*si,facts, scheds) ;
              sources &= my_entities ;
              if((sources & constraints) != constraints) {
                if(MPI_processes == 1) {
                  cerr << "sources & constraints != constraints for input"
                       << endl
                       << sources  << " -- " << *si << endl ;
                } else {
                  debugout << "sources & constraints != constraints for input"
                           << endl
                           << sources  << " -- " << *si << endl ;
                }

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
                      if(MPI_processes == 1) {
                        cerr << "expecting to find variable " << *vi
                             << " at entities " << fails << endl
                             << *vi << " exists at entities " << exist << endl ;
                      } else {
                        debugout << "expecting to find variable " << *vi
                                 << " at entities " << fails << endl
                                 << *vi << " exists at entities " << exist << endl ;
                      }
                    }
                  }
                }
              }
            }
          }
          srcs &= fctrl.constraints ;
        }

        fctrl.control_list.push_back(srcs) ;

        for(size_t i=0;i<fctrl.target_maps.size();++i) {
          entitySet trgts = srcs ;
          for(size_t j=0;j<fctrl.target_maps[i].mapvec.size();++j)
            trgts = fctrl.target_maps[i].mapvec[j]->image(trgts) ;
          tvar_update[fctrl.target_maps[i].v] += trgts ;
          fctrl.generated[fctrl.target_maps[i].v] += trgts ;
        }
#ifdef VERBOSE
        debugout << "recursive rule " << *fi << " generating "
                 << srcs << endl ;
#endif
      }

      recurse_send_entities.push_back(vector<pair<variable,entitySet> >()) ;

      int deltas = 0 ;
      finished = true ;
      for(vi=tvars.begin();vi!=tvars.end();++vi) {
        entitySet tvar_delta = tvar_update[*vi] - tvar_computed[*vi] ;

        entitySet tvar_other = tvar_delta - my_entities ;
        recurse_send_entities.back().push_back(make_pair(*vi,tvar_other)) ;

        entitySet recvset = send_entitySet(tvar_other,facts) ;
        tvar_delta += recvset ;
        tvar_delta += fill_entitySet(tvar_delta,facts) ;
        if(tvar_delta!=EMPTY)
          deltas++ ;
        tvar_update[*vi] += tvar_delta ;
        tvar_computed[*vi] = tvar_update[*vi] ;
      }
      int deltar = 0 ;
      if(!facts.isDistributed())
        deltar = deltas ;
      else
        MPI_Allreduce(&deltas,&deltar, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
      if(deltar != 0)
        finished = false ;
    }

    for(fi=fset.begin();fi!=fset.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
          mi!=fctrl.generated.end();++mi) {

        entitySet create = mi->second;
	if(duplicate_work) {
	  scheds.set_my_proc_able_entities(mi->first, *fi, create);
	  scheds.add_policy(mi->first, sched_db::NEVER);
	}
        create += send_entitySet(create,facts) ;
        create += fill_entitySet(create,facts) ;

        scheds.set_existential_info(mi->first,*fi,create) ;
      }
    }
  }

  void recurse_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    entitySet my_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      my_entities = d->my_entities ;
    }

    ruleSet::const_iterator fi ;
    if(facts.isDistributed()) {
      list<comm_info> request_comm ;
      map<variable, entitySet> orig_requests ;
      for(variableSet::const_iterator vi=recurse_vars.begin();
          vi!=recurse_vars.end();
          ++vi) {
        orig_requests[*vi] = scheds.get_variable_requests(*vi) ;
      }
      request_comm = barrier_process_rule_requests(recurse_vars, facts,  scheds) ;

      vector<pair<variable,entitySet> >::const_iterator vi ;
      vector<pair<variable,entitySet> > send_requested ;

      map<variable,entitySet> var_requests, recurse_entities,recurse_comm ;
      for(variableSet::const_iterator vi=recurse_vars.begin() ;
          vi!=recurse_vars.end();
          ++vi) {
        var_requests[*vi] = scheds.get_variable_requests(*vi) ;
        ruleSet::const_iterator ri ;
        for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
          recurse_entities[*vi] += scheds.get_existential_info(*vi,*ri) ;
        }
      }
      
      std::vector<std::pair<variable,entitySet> > pre_send_entities =
        scheds.get_send_entities(recurse_vars, sched_db::RECURSE_PRE);
      for(vi=pre_send_entities.begin();vi!=pre_send_entities.end();++vi) {
        variable v = vi->first ;
        entitySet send_set = vi->second - recurse_entities[v] ;
        send_requested.push_back(make_pair(v,send_set &
                                           scheds.get_variable_requests(v))) ;
      }
      list<comm_info> pre_plist = put_precomm_info(send_requested, facts) ;
      scheds.update_comm_info_list(pre_plist, sched_db::RECURSE_PRE_PLIST);
      for(fi=recurse_rules.begin();fi!=recurse_rules.end();++fi) {
        fcontrol &fctrl = control_set[*fi] ;
        entitySet control = process_rule_requests(*fi,facts, scheds) ;
        list<entitySet>::iterator ci ;
        entitySet total ;
        for(ci=fctrl.control_list.begin();ci!=fctrl.control_list.end();++ci) {
          *ci -= total ;
          total += *ci ;
        }
      }
      map<rule, list<entitySet>::reverse_iterator> rpos ;


      ruleSet::const_iterator ri ;
      for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri)
        rpos[*ri] = control_set[*ri].control_list.rbegin() ;

      map<variable,vector<entitySet> > recurse_send_req ;
      map<variable,entitySet> all_requests ;
      bool finished = false ;
      do {
        map<variable,entitySet> vreq_map ;
        for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
          entitySet &context = *rpos[*ri] ;

          fcontrol &fctrl = control_set[*ri] ;
          for(size_t i=0;i<fctrl.target_maps.size();++i) {
            entitySet ct = var_requests[fctrl.target_maps[i].v] ;
            for(int j=int(fctrl.target_maps[i].mapvec.size())-1;j>=0;--j)
              ct = fctrl.target_maps[i].mapvec[j]->preimage(ct).first ;
            context &= ct ;
          }
          for(size_t i=0;i<fctrl.recursion_maps.size();++i) {
            entitySet rq = context ;
            for(size_t j=0;j<fctrl.recursion_maps[i].mapvec.size();j++)
              rq = fctrl.recursion_maps[i].mapvec[j]->image(rq) ;
            vreq_map[fctrl.recursion_maps[i].v] += rq  ;
          }


          rpos[*ri]++ ;
          if(rpos[*ri] == control_set[*ri].control_list.rend())
            finished = true ;

          for(variableSet::const_iterator vi = recurse_vars.begin();
              vi!=recurse_vars.end();
              ++vi) {
            all_requests[*vi] += vreq_map[*vi] ;
            entitySet remain = vreq_map[*vi] & recurse_entities[*vi] ;

            remain -= my_entities ;
            recurse_send_req[*vi].push_back(remain) ;
          }
        }
      } while(!finished) ;

      for(variableSet::const_iterator vi = recurse_vars.begin();
          vi!=recurse_vars.end();
          ++vi) {
        std::reverse(recurse_send_req[*vi].begin(),recurse_send_req[*vi].end()) ;
        entitySet req_loc ;
        for(vector<entitySet>::iterator vei=recurse_send_req[*vi].begin();
            vei!=recurse_send_req[*vi].end();
            ++vei) {
          *vei -= req_loc ;
          list<comm_info> req_comm ;
          send_requests(*vei,*vi,facts,req_comm) ;
          send_req_var[*vi].push_back(sort_comm(req_comm,facts)) ;
          req_loc += *vei ;
        }
        recurse_comm[*vi] = req_loc ;
      }

      list<comm_info> pre_req_comm ;
      list<comm_info>::const_iterator li ;
      for(li=request_comm.begin();li!=request_comm.end();++li) {
        entitySet presend = li->send_set - recurse_entities[li->v] ;
        entitySet prerecv = entitySet(li->recv_set) - recurse_entities[li->v] ;
        if(presend != EMPTY || prerecv != EMPTY) {
          comm_info precomm = *li ;
          precomm.send_set = presend ;
          precomm.recv_set = sequence(prerecv) ;
          pre_req_comm.push_back(precomm) ;
        }
      }
      // Add communications for requests that come from recursive rules
      // for results from non-recursive rules.
      for(variableSet::const_iterator vi = recurse_vars.begin();
          vi!= recurse_vars.end();
          ++vi) {
        variable v = *vi ;
        entitySet pre_req = all_requests[v] & ~recurse_entities[v] & ~orig_requests[v]  & ~my_entities ;
        send_requests(pre_req,v,facts,pre_req_comm) ;
      }


      list<comm_info> post_req_comm ;
      for(variableSet::const_iterator vi = recurse_vars.begin();
          vi!= recurse_vars.end();
          ++vi) {
        variable v = *vi ;
        entitySet requests = orig_requests[v] ;
        requests &= recurse_entities[v] ;
        requests -= recurse_comm[*vi] ;

        send_requests(requests,v,facts,post_req_comm) ;
      }
      list<comm_info> pre_clist = sort_comm(pre_req_comm,facts) ;
      list<comm_info> post_clist = sort_comm(post_req_comm,facts) ;
      scheds.update_comm_info_list(pre_clist, sched_db::RECURSE_PRE_CLIST);
      scheds.update_comm_info_list(post_clist, sched_db::RECURSE_POST_CLIST);
    } else {
      for(fi=recurse_rules.begin();fi!=recurse_rules.end();++fi) {
        fcontrol &fctrl = control_set[*fi] ;
        entitySet control = process_rule_requests(*fi,facts, scheds) ;
        list<entitySet>::iterator ci ;
        entitySet total ;
        for(ci=fctrl.control_list.begin();ci!=fctrl.control_list.end();++ci) {
          *ci -= total ;
          total += *ci ;
        }
        do {
          if(fctrl.control_list.back() == EMPTY)
            fctrl.control_list.pop_back() ;
          if(!fctrl.control_list.empty())
            fctrl.control_list.back() &= control ;
        } while(!fctrl.control_list.empty() && fctrl.control_list.back() == EMPTY) ;
      }
    }

  }

  executeP recurse_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds ) {
#ifdef PTHREADS
    ++num_total_recursion;
    bool num_threads_counted = false;
#endif
    CPTR<execute_sequence> el = new execute_sequence ;
    if(facts.isDistributed()) {
      list<comm_info> pre_clist =  scheds.get_comm_info_list(recurse_vars, facts, sched_db::RECURSE_PRE_CLIST);
      list<comm_info> pre_plist =  scheds.get_comm_info_list(recurse_vars, facts, sched_db::RECURSE_PRE_PLIST);
     
      execute_comm2::inc_comm_step() ;
      if(!pre_plist.empty()) {
        //executeP exec_commp = new execute_comm(pre_plist, facts);
        executeP exec_commp2 = new execute_comm2(pre_plist, facts);
        el->append_list(exec_commp2) ;
        //el->append_list(exec_commp) ;
      }

      execute_comm2::inc_comm_step() ;
      if(!pre_clist.empty()) {
        //executeP exec_commc = new execute_comm(pre_clist, facts);
        executeP exec_commc2 = new execute_comm2(pre_clist, facts);
        el->append_list(exec_commc2) ;
        //el->append_list(exec_commc) ;
      }
    }

    map<rule, list<entitySet>::const_iterator> rpos ;
    list<vector<pair<variable,entitySet> > >::const_iterator
      sei = recurse_send_entities.begin() ;
    ruleSet::const_iterator ri ;

    size_t idx = 0 ;
    for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri)
      rpos[*ri] = control_set[*ri].control_list.begin() ;

    bool finished = false ;
    do {
      for(variableSet::const_iterator vi=recurse_vars.begin();
          vi!=recurse_vars.end();
          ++vi) {
        vector<list<comm_info> > &commv = send_req_var[*vi] ;
        execute_comm2::inc_comm_step() ;
        if(idx<commv.size() && commv[idx].size() != 0) {
          //executeP exec_commv = new execute_comm(commv[idx],facts);
          executeP exec_commv2 = new execute_comm2(commv[idx],facts);
          el->append_list(exec_commv2) ;
          //el->append_list(exec_commv) ;
        }
      }
      idx++ ;
      for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
        const fcontrol &fctrl = control_set[*ri] ;
        list<entitySet>::const_iterator &li = rpos[*ri] ;

        if(li==fctrl.control_list.end())
          finished = true ;
        else {
          if(li->size() != 0) {
            executeP exec_rule;
#ifdef PTHREADS
            if (!in_internal_query && threading_recursion) {
              int tnum = thread_control->num_threads();
              int minw = thread_control->min_work_per_thread();
              if (!num_threads_counted) {
                ++num_threaded_recursion;
                num_threads_counted = true;
              }
              if (li->size() >= (size_t)tnum*minw)
                exec_rule = new Threaded_execute_rule
                  (*ri, sequence(*li), facts, scheds);
              else
                exec_rule = new execute_rule
                  (*ri, sequence(*li), facts, scheds);
            } else
              exec_rule = new execute_rule
                (*ri,sequence(*li),facts, scheds);
#else
            exec_rule = new execute_rule(*ri,sequence(*li),facts,scheds);
#endif
            el->append_list(exec_rule);
          }
          li++ ;
        }
      }
      if(!finished) {
        if(facts.isDistributed()) {
          list<comm_info> plist = put_precomm_info(*sei, facts) ;
          execute_comm2::inc_comm_step() ;
          if(!plist.empty()) {
            //executeP exec_comm = new execute_comm(plist,facts);
            executeP exec_comm2 = new execute_comm2(plist,facts);
            el->append_list(exec_comm2) ;
            //el->append_list(exec_comm) ;
          }

          // Make sure to request any variables communicated so that
          // the space is allocated.  This is a hack that should be
          // reworked later.
          for(list<comm_info>::const_iterator li=plist.begin();
              li!=plist.end();
              ++li) {
            entitySet all = li->send_set ;
            all += entitySet(li->recv_set) ;

            scheds.variable_request(li->v,all) ;
          }
        }
        sei++ ;
      }
    } while(!finished) ;

    if(facts.isDistributed()) {
       list<comm_info> post_clist =  scheds.get_comm_info_list(recurse_vars, facts, sched_db::RECURSE_POST_CLIST);
      
      execute_comm2::inc_comm_step() ;
      if(!post_clist.empty()) {
        executeP exec_comm2 = new execute_comm2(post_clist, facts);
        el->append_list(exec_comm2) ;
      }
      // Make sure to request any variables communicated so that
      // the space is allocated.  This is a hack that should be
      // reworked later.
      for(list<comm_info>::const_iterator li=post_clist.begin();
          li!=post_clist.end();
          ++li) {
        entitySet all = li->send_set ;
        all += entitySet(li->recv_set) ;

        scheds.variable_request(li->v,all) ;
      }
    }

    if(el->size() == 0)
      return 0 ;
    else
      return executeP(el) ;
  }
}
