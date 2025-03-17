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
#include "comp_tools.h"
#include "dist_tools.h"
#include "loci_globs.h"
#include "thread.h"
#include <sstream>
#include <set>
#include <map>

using std::ostream ;
using std::endl ;
using std::set ;
using std::map ;
using std::ostringstream ;

namespace Loci {
  //#define DYNAMIC_TIMING ;
#ifdef DYNAMIC_TIMING
#include <execute.h>
  // performance measure
  stopWatch sw_expand ;
  stopWatch sw_expand2 ;
  stopWatch sw_expand_start ;
  stopWatch sw_expand_cache ;
  stopWatch sw_expand_collect_img ;
  stopWatch sw_expand_missing ;
  stopWatch sw_context ;
  stopWatch sw_context_nonepd ;
  stopWatch sw_context_nonepd_domt ;
  stopWatch sw_context_nonepd_ints ;
  stopWatch sw_context_epdend ;
  stopWatch sw_context_epdmid ;
  stopWatch sw_context_epdsta ;
  stopWatch sw_output_oh ;
  stopWatch sw_compute ;
  stopWatch sw_erase ;
  stopWatch sw_invalidate ;
  stopWatch sw_keyremoval ;
  stopWatch sw_insertion ;
  stopWatch sw_keyinsert ;
  stopWatch sw_key_dist ;
  stopWatch sw_dist_renumber ;
  stopWatch sw_dist ;
  stopWatch sw_push ;
  stopWatch sw_expand_comm ;
  stopWatch sw_pw_push ;
  stopWatch sw_param_push ;
  stopWatch sw_param_pack ;
  stopWatch sw_param_unpack ;
  stopWatch sw_param_reduce ;
  stopWatch sw_record_erase ;
  stopWatch sw_dctrl ;

  timeAccumulator ta_expand ;
  timeAccumulator ta_expand2 ;
  timeAccumulator ta_expand_start ;
  timeAccumulator ta_expand_cache ;
  timeAccumulator ta_expand_collect_img ;
  timeAccumulator ta_expand_missing ;
  timeAccumulator ta_context ;
  timeAccumulator ta_context_nonepd ;
  timeAccumulator ta_context_nonepd_domt ;
  timeAccumulator ta_context_nonepd_ints ;
  timeAccumulator ta_context_epdend ;
  timeAccumulator ta_context_epdmid ;
  timeAccumulator ta_context_epdsta ;
  timeAccumulator ta_output_oh ;
  timeAccumulator ta_compute ;
  timeAccumulator ta_erase ;
  timeAccumulator ta_invalidate ;
  timeAccumulator ta_keyremoval ;
  timeAccumulator ta_insertion ;
  timeAccumulator ta_keyinsert ;
  timeAccumulator ta_key_dist ;
  timeAccumulator ta_dist_renumber ;
  timeAccumulator ta_dist ;
  int             ta_dist_number = 0 ;
  timeAccumulator ta_push ;
  timeAccumulator ta_expand_comm ;
  timeAccumulator ta_pw_push ;
  timeAccumulator ta_param_push ;
  timeAccumulator ta_param_pack ;
  timeAccumulator ta_param_unpack ;
  timeAccumulator ta_param_reduce ;
  int             ta_param_reduce_num = 0 ;
  timeAccumulator ta_record_erase ;
  timeAccumulator ta_dctrl ;
  int             ta_drule_executes = 0 ;
  int             ta_dctrl_executes = 0 ;
#endif

  extern bool in_internal_query;
  extern bool threading_pointwise;
  extern int num_threaded_pointwise;
  extern int num_total_pointwise;

  int current_rule_id = 0 ;
  int rule_count = 0;

#ifdef DYNAMICSCHEDULING
  // implementation of execute_dynamic_rule
  namespace {
    variableSet
    remove_synonym(const variableSet& vs, fact_db& facts) {
      variableSet ret ;
      for(variableSet::const_iterator vi=vs.begin();vi!=vs.end();++vi)
        ret += facts.remove_synonym(*vi) ;

      return ret ;
    }
  }

  // first define the I/O operator for the expand type
  std::ostream& operator<<(std::ostream& s, const ExpandStartUnit& u) {
    s << u.var ;
    return s ;
  }
  
  std::ostream& operator<<(std::ostream& s,
                           const std::vector<ExpandStartUnit>& vu) {
    if(vu.empty()) {
      s << "()" ;
      return s ;
    }
    std::vector<ExpandStartUnit>::const_iterator vi = vu.begin() ;
    s << "(" << *vi ;
    for(++vi;vi!=vu.end();++vi)
      s << ", " << *vi ;
    s << ")" ;
    return s ;
  }
  
  std::ostream& operator<<(std::ostream& s, const ExpandUnit& u) {
    s << u.var ;
    return s ;
  }
  
  std::ostream& operator<<(std::ostream& s, const ExpandBlock& b) {
    if(b.units.empty()) {
      s << "()" ;
      return s ;
    }
    std::vector<ExpandUnit>::const_iterator vi = b.units.begin() ;
    s << "(" << *vi ;
    for(++vi;vi!=b.units.end();++vi)
      s << ", " << *vi ;
    s << ")" ;
    return s ;
  }
  
  std::ostream& operator<<(std::ostream& s, const ExpandChain& c) {
    s << "(" << c.tag << ") = " ;
    s << c.expand_start << " -> " ;
    for(std::vector<ExpandBlock>::const_iterator vi=c.expand_chain.begin();
        vi!=c.expand_chain.end();++vi)
      s << *vi << " -> " ;
    s << c.expand_end ;
    return s ;
  }  

  std::ostream& operator<<(std::ostream& s, const NonExpandUnit& u) {
    s << u.var ;
    return s ;
  }
  
  std::ostream& operator<<(std::ostream& s,
                           const std::vector<NonExpandUnit>& vu) {
    if(vu.empty()) {
      s << "()" ;
      return s ;
    }
    std::vector<NonExpandUnit>::const_iterator vi = vu.begin() ;
    s << "(" << *vi ;
    for(++vi;vi!=vu.end();++vi)
      s << ", " << *vi ;
    s << ")" ;
    return s ;
  }

  const int init_cache_reserve = 8192 ;
  bool
  collect_expand_chain(const vmap_info& vmi,
                       const std::string& vmi_tag,
                       rule_implP rp, KeySpaceP space,
                       fact_db& facts, sched_db& scheds,
                       std::vector<ExpandChain>& chains,
                       bool* output_cross_space,
                       std::string* other_space_name,
                       const std::vector<entitySet>** output_ptn,
                       int* target_process_rank,
                       MPI_Comm* target_comm) {
    bool set_extra = (output_cross_space && other_space_name &&
                      output_ptn && target_process_rank && target_comm) ;
    
    // default these extra settings as if no mapping exists
    // we'll modify them later should mapping occur.
    if(set_extra) {
      *output_cross_space = false ;
      *other_space_name = space->get_name() ;
      *output_ptn = &(space->get_key_ptn()) ;
      *target_process_rank = space->get_comm_rank() ;
      *target_comm = space->get_mpi_comm() ;
    }

    if(vmi.mapping.empty()) {
      // no mapping, no chain needed, just initialize rule stores
      for(variableSet::const_iterator vi=vmi.var.begin();
          vi!=vmi.var.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        FATAL(srp == 0) ;
        rp->set_store(*vi, srp) ;
      }
      return false ;
    }
    
    variableSet tunnels = space->get_tunnels() ;
    tunnels = remove_synonym(tunnels,facts) ;
    // when there's mapping present, we need to build a chain
    vector<ExpandStartUnit> expand_start ;
    vector<ExpandBlock> expand_chain ;
    ExpandBlock expand_end ;
    
    bool cross_space = false ;
    string cross_space_name ;
    
    variableSet mappings = vmi.mapping[0] ;
    variableSet mappings_unique ;
    // the first block won't go into the chain,
    // we just set up the rule reps and record the
    // first block as the start point
    for(variableSet::const_iterator vi=mappings.begin();
        vi!=mappings.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      FATAL(srp == 0) ;
      rp->set_store(*vi, srp) ;
      expand_start.push_back(ExpandStartUnit(*vi,srp)) ;
    }
    mappings_unique = remove_synonym(mappings,facts) ;
    if(variableSet(mappings_unique & tunnels) != EMPTY) {
      cross_space = true ;
      // NOTE:
      // we assume if cross_space, then the entire block
      // is tunnel, and later blocks are all in the target space
      // this is pending more investigation
      cross_space_name = space->get_tunnel_space(*(mappings.begin())) ;
      FATAL(cross_space_name == "") ;

      if(set_extra) {
        *output_cross_space = true ;
        if(cross_space_name == "main") {
          *other_space_name = "main" ;
          *output_ptn = &(facts.get_init_ptn(0)) ;
          *target_process_rank = Loci::MPI_rank ;
          *target_comm = MPI_COMM_WORLD ;
        } else {
          map<string,KeySpaceP>::const_iterator ki ;
          ki = facts.keyspace.find(cross_space_name) ;
          if(ki == facts.keyspace.end()) {
            cerr << "Error: collect_expand_chain module"
                 << " keyspace tunnel error!" << endl ;
            cerr << "--queried space name: " << cross_space_name << endl ;
            cerr << "--queried from: " << *(mappings.begin()) << endl ;
            Loci::Abort() ;
          }
          KeySpaceP other_space = ki->second ;
          *other_space_name = other_space->get_name() ;
          *output_ptn = &(other_space->get_key_ptn()) ;
          *target_process_rank = other_space->get_comm_rank() ;
          *target_comm = other_space->get_mpi_comm() ;
        }
      } // end if(set_extra)
    }
    
    const Map* self_pack ;
    const dMap* self_unpack ;
    if(space->has_local_number()) {
      self_pack = &(space->get_l2g_map()) ;
      self_unpack = &(space->get_g2l_map()) ;
    } else {
      self_pack = 0 ;
      self_unpack = 0 ;
    }
        
    // chain starts with the second block
    for(size_t i=1;i<vmi.mapping.size();++i) {
      ExpandBlock expand_block ;
      expand_block.dst_pack = self_pack ;
      expand_block.dst_unpack = self_unpack ;
      
      mappings = vmi.mapping[i] ;
 
      if(cross_space) {
        // figure out the src partition
        if(cross_space_name == "main") {
          const Map* src_pack ;
          const dMap* src_unpack ;
          if(facts.is_distributed_start()) {
            fact_db::distribute_infoP df = facts.get_distribute_info() ;
            src_pack = &(df->l2g) ;
            src_unpack = &(df->g2lv[0]) ;
          } else {
            src_pack = 0 ;
            src_unpack = 0 ;
          }
          expand_block.src_pack = src_pack ;
          expand_block.src_unpack = src_unpack ;
          expand_block.src_ptn = &(facts.get_init_ptn(0)) ;
          expand_block.src_comm = MPI_COMM_WORLD ;
        } else {
          // get the corresponding keyspace impl from fact_db
          map<string,KeySpaceP>::const_iterator ki ;
          ki = facts.keyspace.find(cross_space_name) ;
          if(ki == facts.keyspace.end()) {
            cerr << "Error: collect_expand_chain keyspace "
                 << "tunnel error!" << endl ;
            Loci::Abort() ;
          }
          KeySpaceP other_space = ki->second ;
          const Map* src_pack ;
          const dMap* src_unpack ;
          if(other_space->has_local_number()) {
            src_pack = &(other_space->get_l2g_map()) ;
            src_unpack = &(other_space->get_g2l_map()) ;
          } else {
            src_pack = 0 ;
            src_unpack = 0 ;
          }
          expand_block.src_pack = src_pack ;
          expand_block.src_unpack = src_unpack ;
          expand_block.src_ptn = &(other_space->get_key_ptn()) ;
          expand_block.src_comm = other_space->get_mpi_comm() ;
        }

        // then the src rep comes from the fact_db,
        // and the dst rep comes from the keyspace,
        // and the rule reps need to be keyspace reps
        for(variableSet::const_iterator vi=mappings.begin();
            vi!=mappings.end();++vi) {
          storeRepP src_rep = facts.get_variable(*vi) ;
          storeRepP dst_rep = space->get_shadow(*vi) ;
          FATAL(src_rep == 0) ;
          FATAL(dst_rep == 0) ;
          rp->set_store(*vi, dst_rep) ;
          ExpandUnit eu(*vi,src_rep,dst_rep,init_cache_reserve) ;
          expand_block.units.push_back(eu) ;
        }
      } else {
        expand_block.src_pack = self_pack ;
        expand_block.src_unpack = self_unpack ;
        expand_block.src_ptn = &(space->get_key_ptn()) ;
        expand_block.src_comm = space->get_mpi_comm() ;
        // if not crossing space, then all reps come from fact_db
        for(variableSet::const_iterator vi=mappings.begin();
            vi!=mappings.end();++vi) {
          storeRepP srp = facts.get_variable(*vi) ;
          FATAL(srp == 0) ;
          rp->set_store(*vi, srp) ;
          ExpandUnit eu(*vi,srp,srp,init_cache_reserve) ;
          expand_block.units.push_back(eu) ;
        }
      }
      // evaluate if the next block will be crossing keyspace
      if(!cross_space) {
        mappings_unique = remove_synonym(mappings,facts) ;
        if(variableSet(mappings_unique & tunnels) != EMPTY) {
          cross_space = true ;
          cross_space_name = space->get_tunnel_space(*(mappings.begin())) ;
          FATAL(cross_space_name == "") ;

          if(set_extra) {
            *output_cross_space = true ;
            if(cross_space_name == "main") {
              *other_space_name = "main" ;
              *output_ptn = &(facts.get_init_ptn(0)) ;
              *target_process_rank = Loci::MPI_rank ;
              *target_comm = MPI_COMM_WORLD ;
            } else {
              map<string,KeySpaceP>::const_iterator ki ;
              ki = facts.keyspace.find(cross_space_name) ;
              if(ki == facts.keyspace.end()) {
                cerr << "Error: collect_expand_chain module"
                     << " keyspace tunnel error!" << endl ;
                cerr << "--queried space name: "
                     << cross_space_name << endl ;
                cerr << "--queried from: "
                     << *(mappings.begin()) << endl ;
                Loci::Abort() ;
              }
              KeySpaceP other_space = ki->second ;
              *other_space_name = other_space->get_name() ;
              *output_ptn = &(other_space->get_key_ptn()) ;
              *target_process_rank = other_space->get_comm_rank() ;
              *target_comm = other_space->get_mpi_comm() ;
            }
          } // end if(set_extra)
        }
      }

      expand_chain.push_back(expand_block) ;
    }
    // finally we need to process the var
    // to finish up the expand chain
    if(cross_space) {
      expand_end.dst_pack = self_pack ;
      expand_end.dst_unpack = self_unpack ;
      
      // figure out the src partition
      if(cross_space_name == "main") {
        const Map* src_pack ;
        const dMap* src_unpack ;
        if(facts.is_distributed_start()) {
          fact_db::distribute_infoP df = facts.get_distribute_info() ;
          src_pack = &(df->l2g) ;
          src_unpack = &(df->g2lv[0]) ;
        } else {
          src_pack = 0 ;
          src_unpack = 0 ;
        }
        expand_end.src_pack = src_pack ;
        expand_end.src_unpack = src_unpack ;
        expand_end.src_ptn = &(facts.get_init_ptn(0)) ;
        expand_end.src_comm = MPI_COMM_WORLD ;
      } else {
        // get the corresponding keyspace impl from fact_db
        map<string,KeySpaceP>::const_iterator ki ;
        ki = facts.keyspace.find(cross_space_name) ;
        if(ki == facts.keyspace.end()) {
          cerr << "Error: collect_expand_chain keyspace "
               << "tunnel error!" << endl ;
          Loci::Abort() ;
        }
        KeySpaceP other_space = ki->second ;
        const Map* src_pack ;
        const dMap* src_unpack ;
        if(other_space->has_local_number()) {
          src_pack = &(other_space->get_l2g_map()) ;
          src_unpack = &(other_space->get_g2l_map()) ;
        } else {
          src_pack = 0 ;
          src_unpack = 0 ;
        }
        expand_end.src_pack = src_pack ;
        expand_end.src_unpack = src_unpack ;
        expand_end.src_ptn = &(other_space->get_key_ptn()) ;
        expand_end.src_comm = other_space->get_mpi_comm() ;
      }
      
      for(variableSet::const_iterator vi=vmi.var.begin();
          vi!=vmi.var.end();++vi) {
        storeRepP src_rep = facts.get_variable(*vi) ;
        storeRepP dst_rep = space->get_shadow(*vi) ;
        FATAL(src_rep == 0) ;
        FATAL(dst_rep == 0) ;
        rp->set_store(*vi, dst_rep) ;
        ExpandUnit eu(*vi,src_rep,dst_rep,init_cache_reserve) ;
        expand_end.units.push_back(eu) ;
      }
    } else {
      expand_end.src_pack = self_pack ;
      expand_end.src_unpack = self_unpack ;
      expand_end.src_ptn = &(space->get_key_ptn()) ;
      expand_end.src_comm = space->get_mpi_comm() ;
      // if not crossing space, then all reps come from fact_db
      for(variableSet::const_iterator vi=vmi.var.begin();
          vi!=vmi.var.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        FATAL(srp == 0) ;
        rp->set_store(*vi, srp) ;
        ExpandUnit eu(*vi,srp,srp,init_cache_reserve) ;
        expand_end.units.push_back(eu) ;
      }
    }
    
    ExpandChain c ;
    c.tag = vmi_tag ;
    c.expand_start = expand_start ;
    c.expand_chain = expand_chain ;
    c.expand_end = expand_end ;
    
    chains.push_back(c) ;

    return true ;
  }

  bool
  collect_nonexpand_chain(const vmap_info& vmi,
                          rule_implP rp, KeySpaceP space,
                          fact_db& facts, sched_db& scheds,
                          std::vector<NonExpandUnit>& chain) {
    if(!vmi.mapping.empty())
      return false ;

    variableSet space_vars = space->get_control_vars() ;
    for(variableSet::const_iterator vi=vmi.var.begin();
        vi!=vmi.var.end();++vi) {
      if(!space_vars.inSet(*vi)) {
        continue ;              // ignore non local space variables
      }
      
      storeRepP srp = facts.get_variable(*vi) ;
      FATAL(srp == 0) ;

      chain.push_back(NonExpandUnit(*vi,srp)) ;
    }
    return true ;
  }

  // if start set exists we'll use it, otherwise the start
  // is computed from the chain itself by the domain of the
  // start block
  entitySet
  execute_expand_chain(ExpandChain& chain, entitySet* start) {
#ifdef DYNAMIC_TIMING
    sw_expand2.start() ;
#endif
    // first we need to compute the starting request domain
    entitySet exist ;
#ifdef DYNAMIC_TIMING
    sw_expand_start.start() ;
#endif
    for(vector<ExpandStartUnit>::iterator mi=chain.expand_start.begin();
        mi!=chain.expand_start.end();++mi) {
      MapRepP mp = MapRepP(mi->rep->getRep()) ;
      FATAL(mp == 0) ;
      if(start)
        exist += mp->image(*start) ;
      else
        exist += mp->image(mp->domain()) ;
    }
#ifdef DYNAMIC_TIMING
    ta_expand_start.addTime(sw_expand_start.stop(), 1) ;
#endif
    for(vector<ExpandBlock>::iterator
          block=chain.expand_chain.begin();
        block!=chain.expand_chain.end();++block) {
      vector<storeRepP> src, dst ;
      vector<entitySet> dst_dom, dst_needed ;
      entitySet block_missing ;
      entitySet image ;
      
      for(vector<ExpandUnit>::iterator unit=block->units.begin();
          unit!=block->units.end();++unit) {
#ifdef DYNAMIC_TIMING
        sw_expand_missing.start() ;
#endif
        // figure out the missing part for dst store
        entitySet dom = unit->dst->domain() ;
        const entitySet& needed = exist ;
        entitySet missing = needed - dom ;
#ifdef DYNAMIC_TIMING
        ta_expand_missing.addTime(sw_expand_missing.stop(), 1) ;
#endif
        // if no one is missing, then continue with nothing to do.
        if(GLOBAL_AND(missing == EMPTY)) {
#ifdef DYNAMIC_TIMING
          sw_expand_collect_img.start() ;
#endif
          MapRepP mp = MapRepP(unit->dst->getRep()) ;
          FATAL(mp == 0) ;
          image += mp->image(needed) ;
#ifdef DYNAMIC_TIMING
          ta_expand_collect_img.addTime(sw_expand_collect_img.stop(), 1) ;
#endif
          continue ;
        }
#ifdef DYNAMIC_TIMING
        sw_expand_cache.start() ;
#endif
        // otherwise, get the missing part
        // see if the reserve limit is reached
        const entitySet& cached_domain = dom ;
        int available_size = unit->dst_reserve_size -
          cached_domain.size() ;
        if((int)missing.size() > available_size) {
          // need to purge all old data
          entitySet invalid = cached_domain - needed ;
          unit->dst->erase(invalid) ;
          available_size += invalid.size() ;
          // if still not enough space, then we need
          // to increase our reservation limit by a factor of 2
          if((int)missing.size() > available_size) {
            unit->dst_reserve_size = needed.size() * 2 ;
          }
        }
#ifdef DYNAMIC_TIMING
        ta_expand_cache.addTime(sw_expand_cache.stop(), 1) ;
#endif
        block_missing += missing ;
        src.push_back(unit->src) ;
        dst.push_back(unit->dst) ;
        dst_dom.push_back(dom) ;
        dst_needed.push_back(needed) ;
      } // end all Map unit expansion
      // then expand those that need to
#ifdef DYNAMIC_TIMING
      sw_expand_comm.start() ;
#endif
      vector<entitySet> real_missing =
        expand_store2(src,/*src rep*/
                      block->src_pack, block->src_unpack,
                      dst,/*dst rep*/
                      block->dst_unpack, block_missing,
                      *(block->src_ptn), block->src_comm) ;
#ifdef DYNAMIC_TIMING
      ta_expand_comm.addTime(sw_expand_comm.stop(), 1) ;
      
      sw_expand_collect_img.start() ;
#endif
      // gather image (next exist domain)
      for(size_t k=0;k<src.size();++k) {
        MapRepP mp = MapRepP(dst[k]->getRep()) ;
        FATAL(mp == 0) ;
        entitySet not_missing = dst_dom[k] & dst_needed[k] ;
        image += mp->image(not_missing + real_missing[k]) ;
      }
#ifdef DYNAMIC_TIMING
      ta_expand_collect_img.addTime(sw_expand_collect_img.stop(), 1) ;
#endif
      // end of one block expansion
      // then reset exist domain
      exist = image ;
    }   // end all block expansion
    // then we need to expand the end block, if it is not a target chain
    if(chain.tag != "target") {
      vector<storeRepP> src, dst ;
      entitySet block_missing ;
      ExpandBlock& block = chain.expand_end ;
      for(vector<ExpandUnit>::iterator unit=block.units.begin();
          unit!=block.units.end();++unit) {
#ifdef DYNAMIC_TIMING
        sw_expand_missing.start() ;
#endif
        // figure out the missing part for dst store
        const entitySet& needed = exist ;
        entitySet dom = unit->dst->domain() ;
        entitySet missing = needed - dom ;
#ifdef DYNAMIC_TIMING
        ta_expand_missing.addTime(sw_expand_missing.stop(), 1) ;
#endif
        // if no one is missing, then continue with nothing to do.
        //if(GLOBAL_AND(missing == EMPTY))
        //continue ;
#ifdef DYNAMIC_TIMING
        sw_expand_cache.start() ;
#endif
        // otherwise, get the missing part
        // see if the reserve limit is reached
        const entitySet& cached_domain = dom ;
        int available_size = unit->dst_reserve_size -
          cached_domain.size() ;
        if((int)missing.size() > available_size) {
          // need to purge all old data
          entitySet invalid = cached_domain - needed ;
          unit->dst->erase(invalid) ;
          available_size += invalid.size() ;
          // if still not enough space, then we need
          // to increase our reservation limit by a factor of 2
          if((int)missing.size() > available_size) {
            unit->dst_reserve_size = needed.size() * 2 ;
          }
        }
#ifdef DYNAMIC_TIMING
        ta_expand_cache.addTime(sw_expand_cache.stop(), 1) ;
#endif
        block_missing += missing ;
        src.push_back(unit->src) ;
        dst.push_back(unit->dst) ;
      } // end all Map unit expansion
#ifdef DYNAMIC_TIMING
      sw_expand_comm.start() ;
#endif
      expand_store2(src,/*src rep*/ block.src_pack, block.src_unpack,
                    dst,/*dst rep*/ block.dst_unpack, block_missing,
                    *(block.src_ptn), block.src_comm) ;
#ifdef DYNAMIC_TIMING
      ta_expand_comm.addTime(sw_expand_comm.stop(), 1) ;
#endif
    }
    
#ifdef DYNAMIC_TIMING
    ta_expand2.addTime(sw_expand2.stop(),1) ;
#endif
    return exist ;
  }
  
  entitySet
  dynamic_rule_context(const rule& rule_tag, KeySpaceP space,
                       fact_db& facts, sched_db& scheds,
                       const std::vector<ExpandChain>& input_chains,
                       const std::vector<NonExpandUnit>& input_nes,
                       const std::vector<NonExpandUnit>& input_nec) {
    entitySet context ;
    
    entitySet sources = space->get_keys_local() ;
    entitySet constraints = space->get_keys_local() ;
    const entitySet& all_keys = space->get_keys_local() ;

#ifdef DYNAMIC_TIMING
    sw_context_nonepd.start() ;
#endif
    // first process the non expanding variables
    for(vector<NonExpandUnit>::const_iterator vi=input_nes.begin();
        vi!=input_nes.end();++vi) {
#ifdef DYNAMIC_TIMING
      sw_context_nonepd_domt.start() ;
#endif
      entitySet es = vi->var_rep->domain() ;
#ifdef DYNAMIC_TIMING
      ta_context_nonepd_domt.addTime(sw_context_nonepd_domt.stop(), 1) ;

      sw_context_nonepd_ints.start() ;
#endif
      sources &= es ;
#ifdef DYNAMIC_TIMING
      ta_context_nonepd_ints.addTime(sw_context_nonepd_ints.stop(), 1) ;
#endif
    }
    for(vector<NonExpandUnit>::const_iterator vi=input_nec.begin();
        vi!=input_nec.end();++vi) {
#ifdef DYNAMIC_TIMING
      sw_context_nonepd_domt.start() ;
#endif
      entitySet es = vi->var_rep->domain() ;
#ifdef DYNAMIC_TIMING
      ta_context_nonepd_domt.addTime(sw_context_nonepd_domt.stop(), 1) ;

      sw_context_nonepd_ints.start() ;
#endif
      constraints &= es ;
#ifdef DYNAMIC_TIMING
      ta_context_nonepd_ints.addTime(sw_context_nonepd_ints.stop(), 1) ;
#endif
    }
#ifdef DYNAMIC_TIMING
    ta_context_nonepd.addTime(sw_context_nonepd.stop(), 1) ;
#endif
    // then work on the chains
    for(vector<ExpandChain>::const_iterator chain=input_chains.begin();
        chain!=input_chains.end();++chain) {
      // // DEBUG
      // for(vector<ExpandStartUnit>::const_iterator
      //       mi=chain->expand_start.begin();
      //     mi!=chain->expand_start.end();++mi) {
      //   if(chain->tag == "source") {
      //     sources &= mi->rep->domain() ;
      //   } else if(chain->tag == "constraint") {
      //     constraints &= mi->rep->domain() ;
      //   }
      // }
      // continue ;
      // ////////
      // first work on final block
      entitySet request = ~EMPTY ;

#ifdef DYNAMIC_TIMING
      sw_context_epdend.start() ;
#endif
      for(vector<ExpandUnit>::const_iterator
            unit=chain->expand_end.units.begin();
          unit!=chain->expand_end.units.end();++unit) {
        request &= unit->dst->domain() ;
      }
#ifdef DYNAMIC_TIMING
      ta_context_epdend.addTime(sw_context_epdend.stop(), 1) ;
      
      sw_context_epdmid.start() ;
#endif
      // the working on the expanding blocks in reverse order
      for(vector<ExpandBlock>::const_reverse_iterator
            block=chain->expand_chain.rbegin();
          block!=chain->expand_chain.rend();++block) {
        entitySet working = ~EMPTY ;
        for(vector<ExpandUnit>::const_iterator
              unit=block->units.begin();
            unit!=block->units.end();++unit) {
          MapRepP mp = MapRepP(unit->dst->getRep()) ;
          FATAL(mp == 0) ;
          working &= mp->preimage(request).first ;
        }
        request = working ;
      }
#ifdef DYNAMIC_TIMING
      ta_context_epdmid.addTime(sw_context_epdmid.stop(), 1) ;

      sw_context_epdsta.start() ;
#endif
      // then work on the start block
      {
        entitySet working = ~EMPTY ;
        for(vector<ExpandStartUnit>::const_iterator
              mi=chain->expand_start.begin();
            mi!=chain->expand_start.end();++mi) {
          MapRepP mp = MapRepP(mi->rep->getRep()) ;
          FATAL(mp == 0) ;
          working &= mp->preimage(request).first ;
        }
        request = working ;
      }
#ifdef DYNAMIC_TIMING
      ta_context_epdsta.addTime(sw_context_epdsta.stop(), 1) ;
#endif
      if(chain->tag == "source") {
        sources &= request ;
      } else if(chain->tag == "constraint") {
        constraints &= request ;
      }
    }

    // check constraint
    const rule_impl::info& rinfo = rule_tag.get_info().desc ;    
    if(!rinfo.constraints.empty() && constraints != all_keys) {
      if(constraints-sources != EMPTY) {
        if(MPI_processes == 1) {
          cerr << "Error: dynamic rule: " << rule_tag
               << " cannot supply all constraints" << endl ;
          cerr << "--Constraint: " << constraints << endl ;
          cerr << "--Sources   : " << sources << endl ;
          Loci::Abort() ;
        } else {
          debugout << "Error: dynamic rule: " << rule_tag
                   << " cannot supply all constraints" << endl ;
          debugout << "--Constraint: " << constraints << endl ;
          debugout << "--Sources   : " << sources << endl ;
          if(MPI_rank == 0)
            cerr << "Dynamic rule: " << rule_tag
                 << " error, see debug file for more info" << endl ;
          Loci::Abort() ;
        }
      }
      context = sources & constraints ;
    } else
      context = sources ;

    return context ;
  }

  execute_dynamic_rule::
  execute_dynamic_rule(rule r, KeySpaceP kp,
                       fact_db& facts, sched_db& scheds) {
    rule_tag = r ;
    rp = r.get_rule_implP() ;
    space = kp ;

    large_context_size = 2 * facts.get_init_ptn(0)[0].size() ;

    dflag = true ;
    // register itself to the keyspace
    space->register_dcontrol_rule(rule_tag, this) ;

    set<vmap_info>::const_iterator vmsi ;
    // figure out the expand chain
    // and also initialize rule storeReps
    for(vmsi=rule_tag.get_info().desc.sources.begin();
        vmsi!=rule_tag.get_info().desc.sources.end();++vmsi) {
      collect_expand_chain(*vmsi, "source", rp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      // then collect non expanding chains
      collect_nonexpand_chain(*vmsi, rp, space,
                              facts, scheds, input_ne_source) ;
    } // end of source

    for(vmsi=rule_tag.get_info().desc.constraints.begin();
        vmsi!=rule_tag.get_info().desc.constraints.end();++vmsi) {
      collect_expand_chain(*vmsi, "constraint", rp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, rp, space,
                              facts, scheds, input_ne_constraint) ;
    } // end of constraints

    // collect targets
    for(vmsi=rule_tag.get_info().desc.targets.begin();
        vmsi!=rule_tag.get_info().desc.targets.end();++vmsi) {
      variableSet assign ;
      if(!vmsi->assign.empty()) {
        for(size_t i=0;i<vmsi->assign.size();++i) {
          assign += vmsi->assign[i].first ;
        }
      }
      // collect chains of output and output_info
      bool chain_collected =
        collect_expand_chain(*vmsi, "target", rp, space,
                             facts, scheds, output_chains,
                             &output_cross_space,
                             &other_space_name,
                             &output_ptn,
                             &target_process_rank,
                             &target_comm) ;
      // first check to see if output is crossing space
      if(output_cross_space) {
        cerr << "Error: targets of pointwise rule is NOT allowed"
             << " to cross keyspace. Offending rule: " << rule_tag
             << endl ;
        Loci::Abort() ;
      }
      // then construct the output_info
      OutInfo oi ;
      if(chain_collected) {
        oi.expand = true ;
        oi.expand_index = output_chains.size() - 1 ;
      } else {
        oi.expand = false ;
      }
      oi.targets = vmsi->var ;
      variableSet alloc_vars = variableSet(vmsi->var - assign) ;
      // exclude OUTPUT from target alloc
      variableSet output_vars ;
      for(variableSet::const_iterator vi=alloc_vars.begin();
          vi!=alloc_vars.end();++vi) {
        variable sv(*vi, time_ident()) ;
        if(sv == variable("OUTPUT"))
          output_vars += *vi ;
      }
      alloc_vars -= output_vars ;
      oi.targets_alloc = alloc_vars ;
      output_info.push_back(oi) ;
    }
  }

  void execute_dynamic_rule::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;

#ifdef DYNAMIC_TIMING
    ++ta_drule_executes ;
#endif
    // if dflag is not set, then we skip the expansion and context
    if(dflag) {
#ifdef DYNAMIC_TIMING
      ++ta_dctrl_executes ;
      sw_expand.start() ;
#endif
      // first we need to communicate input_chains    
      for(vector<ExpandChain>::iterator chain=input_chains.begin();
          chain!=input_chains.end();++chain) {
        execute_expand_chain(*chain, 0) ;
      }     // end all chain expansion
#ifdef DYNAMIC_TIMING
      ta_expand.addTime(sw_expand.stop(), 1) ;
      
      sw_context.start() ;
#endif
      
      // then we need to compute the execution context
      context =
        dynamic_rule_context(rule_tag, space,
                             facts, scheds, input_chains,
                             input_ne_source, input_ne_constraint) ;
#ifdef DYNAMIC_TIMING
      ta_context.addTime(sw_context.stop(), 1) ;
#endif
    }

    // // give a warning if context is very large
    // size_t context_size = context.size() ;
    // if(context_size > large_context_size) {
    //   if(space->get_comm_rank() == 0)
    //     cerr << "Warning: dynamic rule: " << rule_tag
    //          << " executing over a very large context, potentially"
    //          << " fatal to memory allocation" << endl ;
    // }

#ifdef DYNAMIC_TIMING
    sw_output_oh.start() ;
#endif
    // if dflag is not set, we can skip the process
    if(dflag) {
      // then process the outputs
      for(size_t i=0;i<output_info.size();++i) {
        entitySet target_exist ;
        if(output_info[i].expand) {
          // expand the chain
          ExpandChain& oc = output_chains[output_info[i].expand_index] ;
          target_exist = execute_expand_chain(oc, &context) ;
          
          // figure out the point 2 point comm structure
          vector<P2pCommInfo> send, recv ;
          entitySet push = target_exist -
            (*output_ptn)[target_process_rank] ;
          get_p2p_comm(*output_ptn, push, 0, 0,/*we don't use remap*/
                       target_comm, recv, send) ;
          send.swap(output_info[i].send) ;
          recv.swap(output_info[i].recv) ;
          // add in the recving domain
          for(size_t t=0;t<recv.size();++t)
            target_exist += recv[t].local_dom ;
        } else
          target_exist = context ;
        
        output_info[i].target_exist = target_exist ;
        
        // NOTE, we don't actually allocate the targets, in fact,
        // calling the allocate(target_exist) for every target
        // is WRONG because target_exist is just the domain
        // existed for the target variable for this rule context.
        // allocating it over only this portion will destroy
        // other potential useful domains. And so we don't allocate
        // dynamic targets (since they are dynamic, they can
        // self allocate).
        
        const variableSet& vars = output_info[i].targets_alloc ;
        for(variableSet::const_iterator vi=vars.begin();
            vi!=vars.end();++vi) {
          storeRepP srp = facts.get_variable(*vi) ;
          FATAL(srp == 0) ;
          srp->guarantee_domain(target_exist) ;
        }
        
      }
    } // end if(!dflag)
      
#ifdef DYNAMIC_TIMING
    ta_output_oh.addTime(sw_output_oh.stop(), 1) ;

    sw_compute.start() ;
#endif
    
    // finally call the compute method
    rp->compute(sequence(context));

#ifdef DYNAMIC_TIMING
    ta_compute.addTime(sw_compute.stop(), 1) ;
    
    sw_pw_push.start() ;
#endif
    // then we'll need to see if we will want to
    // communicate the result should there be maps
    // involved in some of the targets
    for(size_t i=0;i<output_info.size();++i) {
      if(output_info[i].expand) {
        const variableSet& vars = output_info[i].targets ;
        const vector<P2pCommInfo>& send = output_info[i].send ;
        const vector<P2pCommInfo>& recv = output_info[i].recv ;
        for(variableSet::const_iterator vi=vars.begin();
            vi!=vars.end();++vi) {
          storeRepP srp = facts.get_variable(*vi) ;
          FATAL(srp == 0) ;
          fill_store(srp, 0/*no remap*/, srp, 0/*no remap*/,
                     send, recv, target_comm) ;
        }
      }
    }
#ifdef DYNAMIC_TIMING
    ta_pw_push.addTime(sw_pw_push.stop(), 1) ;
#endif
    // if dflag is set, then reset it
    if(dflag) dflag = false ;
    
    timer.addTime(s.stop(),0) ;
  }

  void
  execute_dynamic_rule::Print(ostream& s) const {
    printIndent(s) ;
    s << rule_tag << " over [dynamic] sequence " << endl ;
  }

  void
  execute_dynamic_rule::dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "rule: "<<rule_tag ;

    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }
#endif
  
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts, const sched_db &scheds)  {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    exec_size = seq.size() ;
    //
#ifdef PAPI_DEBUG
    papi_events[0] = PAPI_L1_DCM;
    papi_events[1] = PAPI_L2_DCM;
    papi_values[0] = papi_values[1] = 0;
    l1_dcm = l2_dcm = 0;
#endif
  }

  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts,
                             variable v, const storeRepP &p, const sched_db &scheds) {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    rp->set_store(v,p) ;
    exec_seq = seq ;
    exec_size = seq.size() ;
#ifdef PAPI_DEBUG
    papi_events[0] = PAPI_L1_DCM;
    papi_events[1] = PAPI_L2_DCM;
    papi_values[0] = papi_values[1] = 0;
    l1_dcm = l2_dcm = 0;
#endif
  }

  void execute_rule::execute(fact_db &facts, sched_db &scheds) {
#ifdef PAPI_DEBUG
    if( (PAPI_start_counters(papi_events,2)) != PAPI_OK) {
      cerr << "PAPI failed to start counters" << endl;
      Loci::Abort();
    }
#endif

    stopWatch s ;
    s.start() ;
    current_rule_id = rule_tag.ident() ;
#ifdef VERBOSE
    Loci::debugout << "executing " << rule_tag << endl ;
#endif
    //rp->compute(exec_seq);
    execute_prelude(exec_seq);
    execute_kernel(exec_seq);
    execute_postlude(exec_seq);
    current_rule_id = 0 ;
    timer.addTime(s.stop(),exec_size) ;
#ifdef PAPI_DEBUG
    if( (PAPI_stop_counters(papi_values,2)) != PAPI_OK) {
      cerr << "PAPI failed to read counters" << endl;
      Loci::Abort();
    }
    l1_dcm += papi_values[0];
    l2_dcm += papi_values[1];
#endif
  }

  inline void 
  execute_rule::execute_kernel(const sequence& s)
  { rp->compute(s); }
  
  inline void
  execute_rule::execute_prelude(const sequence& s)
  { rp->prelude(s); }

  inline void
  execute_rule::execute_postlude(const sequence& s)
  { rp->postlude(s); }

  void execute_rule::Print(ostream &s) const {
    printIndent(s) ;
    s << rule_tag << " over sequence " ;
    if(verbose || exec_seq.num_intervals() < 4) {
      s << exec_seq << endl ;
    } else {
      s << "[ ... ], l=" << exec_seq.size() << endl ;
    }
  }

  void execute_rule::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "rule: "<<rule_tag ;

    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;

#ifdef PAPI_DEBUG
    data_collector.accumulateSchedMemory(oss.str(),
                                         8*exec_seq.num_intervals());

    data_collector.accumulateDCM(oss.str(),
                                 papi_values[0], papi_values[1]);
#endif
  }

  void impl_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    existential_rule_analysis(impl,facts, scheds) ;
  }

  void impl_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    entitySet exec_seq = process_rule_requests(impl,facts, scheds) ;
    scheds.update_exec_seq(impl, exec_seq);
  }

  executeP impl_compiler::
  create_execution_schedule(fact_db &facts,sched_db &scheds) {
    //    if(GLOBAL_AND(exec_seq.size()==0)) {
    //      return executeP(0) ;
    //    }
    variableSet targets = impl.targets() ;
    WARN(targets.size() == 0) ;

    extern int method ;
    entitySet exec_seq = scheds.get_exec_seq(impl);
    if (impl.get_info().rule_impl->dynamic_schedule_rule() &&
        use_dynamic_scheduling) {
      executeP execute_dynamic =
        new dynamic_schedule_rule(impl,exec_seq,facts, scheds,method) ;

      return execute_dynamic;
    }

#ifdef PTHREADS
    ++num_total_pointwise;
    bool threadable = 
      impl.get_info().rule_impl->thread_rule() &&
      (impl.get_info().rule_impl->get_rule_class() 
                       == rule_impl::POINTWISE ||
       impl.get_info().rule_impl->get_rule_class()
                       == rule_impl::UNIT);
    rule_implP ti = impl.get_rule_implP() ;
    for (variableSet::const_iterator vi=targets.begin();
        vi!=targets.end();++vi) {
      storeRepP tr = ti->get_store(*vi);
      if (isPARAMETER(tr)) {
        threadable = false;
        break;
      }
    }
    if(!in_internal_query && threading_pointwise && threadable) {
      int tnum = thread_control->num_threads();
      int minw = thread_control->min_work_per_thread();
      // if a rule is not for threading, then generate a normal module,
      // also no multithreading if the execution sequence is too small
      if(exec_seq.size() < (size_t)tnum*minw)
        // normal case
        return new execute_rule(impl,sequence(exec_seq),facts, scheds);
      else {
        // generate multithreaded execution module
        ++num_threaded_pointwise;
        return new Threaded_execute_rule(impl, exec_seq, facts, scheds);
      }
    }
#endif
    // normal case
    executeP exec_rule =
      new execute_rule(impl,sequence(exec_seq),facts, scheds);

    return exec_rule;
  }

#ifdef DYNAMICSCHEDULING
  void
  dynamic_impl_compiler::
  set_var_existence(fact_db& facts, sched_db& scheds) {
    // since the dynamic_impl_compiler will
    // evalulate the rules at the runtime, we'll
    // just set everything (existence and request) to
    // be EMPTY. Note: this is also important to interface
    // correctly with the current allocation code and the
    // barrier and reduction compiler. otherwise they'll
    // do some extra job that we don't want.
    variableSet targets = impl.targets() ;
    for(variableSet::const_iterator vi=targets.begin();
        vi!=targets.end();++vi)
      scheds.set_existential_info(*vi, impl, EMPTY) ;
  }

  void
  dynamic_impl_compiler::
  process_var_requests(fact_db& facts, sched_db& scheds) {
    variableSet targets = impl.targets() ;
    for(variableSet::const_iterator vi=targets.begin();
        vi!=targets.end();++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;
    variableSet sources = impl.sources() ;
    for(variableSet::const_iterator vi=sources.begin();
        vi!=sources.end();++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;
  }

  executeP
  dynamic_impl_compiler::create_execution_schedule(fact_db& facts,
                                                   sched_db& scheds) {
    executeP execute = new execute_dynamic_rule(impl,space,facts,scheds) ;
    return execute ;
  }
  
  
  // execute_dynamic_applyrule module implementation
  execute_dynamic_applyrule::
  execute_dynamic_applyrule(rule a, rule u, KeySpaceP kp,
                            fact_db& facts, sched_db& scheds) {    
    apply_tag = a ;
    unit_tag = u ;
    rp = apply_tag.get_rule_implP() ;
    space = kp ;

    large_context_size = 2 * facts.get_init_ptn(0)[0].size() ;

    dflag = true ;
    // register itself to the keyspace
    space->register_dcontrol_rule(apply_tag, this) ;

    FATAL(apply_tag.targets().size() != 1) ;

    set<vmap_info>::const_iterator vmsi ;
    // figure out the expand chain
    // and also initialize rule storeReps

    // build the input chain first
    for(vmsi=apply_tag.get_info().desc.sources.begin();
        vmsi!=apply_tag.get_info().desc.sources.end();++vmsi) {
      collect_expand_chain(*vmsi, "source", rp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, rp, space,
                              facts, scheds, input_ne_source) ;
    } // end of source

    for(vmsi=apply_tag.get_info().desc.constraints.begin();
        vmsi!=apply_tag.get_info().desc.constraints.end();++vmsi) {
      collect_expand_chain(*vmsi, "constraint", rp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, rp, space,
                              facts, scheds, input_ne_constraint) ;
    } // end of constraints

    // for apply rule, we also need to build a output chain
    // we only have one target
    vector<ExpandChain> output_chains ;
    vmsi = apply_tag.get_info().desc.targets.begin() ;
    collect_expand_chain(*vmsi, "target", rp, space,
                         facts, scheds, output_chains,
                         &output_cross_space,
                         &other_space_name,
                         &output_ptn,
                         &target_process_rank,
                         &target_comm) ;
    if(!output_chains.empty()) {
      output_chain = output_chains[0] ;
    }
    if(vmsi->mapping.empty())
      output_map = false ;
    else
      output_map = true ;

    target_var = *(apply_tag.targets().begin()) ;
    target_rep_in_facts = facts.get_variable(target_var) ;
    FATAL(target_rep_in_facts == 0) ;
    if(!output_cross_space)
      rp->set_store(target_var, target_rep_in_facts) ;

    if(output_cross_space) {
      // find out the target remap
      if(other_space_name == "main") {
        if(facts.is_distributed_start()) {
          fact_db::distribute_infoP df = facts.get_distribute_info() ;
          target_pack = &(df->l2g) ;
          target_unpack = &(df->g2lv[0]) ;
        } else {
          target_pack = 0 ;
          target_unpack = 0 ;
        }
      } else {
        map<string,KeySpaceP>::const_iterator ki ;
        ki = facts.keyspace.find(other_space_name) ;
        if(ki == facts.keyspace.end()) {
          cerr << "Error: execute_dynamic_applyrule module keyspace "
               << "tunnel error!" << endl ;
          Loci::Abort() ;
        }
        KeySpaceP other_space = ki->second ;
        if(other_space->has_local_number()) {
          target_pack = &(other_space->get_l2g_map()) ;
          target_unpack = &(other_space->get_g2l_map()) ;
        } else {
          target_pack = 0 ;
          target_unpack = 0 ;
        }
      }
    } else {
      target_pack = 0 ;
      target_unpack = 0 ;
    }
    
  }

  void execute_dynamic_applyrule::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;

    if(dflag) {
#ifdef DYNAMIC_TIMING
      sw_expand.start() ;
#endif
      // first we need to communicate input_chains
      for(vector<ExpandChain>::iterator chain=input_chains.begin();
          chain!=input_chains.end();++chain) {
        execute_expand_chain(*chain, 0) ;
      }     // end all chain expansion
#ifdef DYNAMIC_TIMING
      ta_expand.addTime(sw_expand.stop(), 1) ;
      
      sw_context.start() ;
#endif
      
      // then we need to compute the execution context
      context =
        dynamic_rule_context(apply_tag,space,
                             facts,scheds,input_chains,
                             input_ne_source, input_ne_constraint);
#ifdef DYNAMIC_TIMING
      ta_context.addTime(sw_context.stop(), 1) ;
#endif
    } // end if(dflag)
    // // give a warning if context is very large
    // size_t context_size = context.size() ;
    // if(context_size > large_context_size) {
    //   if(space->get_comm_rank() == 0)
    //     cerr << "Warning: dynamic rule: " << apply_tag
    //          << " executing over a very large context, potentially"
    //          << " fatal to memory allocation" << endl ;
    // }

    // we then go ahead to expand the output chain
    // if there is no mapping in the output
    // then the target exist on context, otherwise
    // we need to expand the output chain to determine
    // what the target exist on
    if(dflag) {
#ifdef DYNAMIC_TIMING
      sw_expand.start() ;
#endif
      if(!output_map)
        target_var_exist = context ;
      else
        target_var_exist = execute_expand_chain(output_chain, &context) ;
#ifdef DYNAMIC_TIMING
      ta_expand.addTime(sw_expand.stop(), 1) ;
#endif
    }

    // we can then figure out the entitySet needs to be sent
    if(dflag) {
#ifdef DYNAMIC_TIMING
      sw_output_oh.start() ;
#endif
      if(output_cross_space) {
        // for cross space apply rule, we need to
        // push everything the target exists since
        // we compute the target locally
        entitySet push = target_var_exist ;
        // get the comm, the direction is reversed since
        // we are doing pushing operation
        get_p2p_comm(*output_ptn, push,
                     target_unpack, 0, target_comm, recv, send) ;
      } else {
        // for local space reduction, we only push those that
        // are outside of local process's keyset
        entitySet push =
          target_var_exist - (*output_ptn)[target_process_rank] ;
        // get the comm without remapping since we are in the
        // local dynamic keyspace
        get_p2p_comm(*output_ptn, push, 0, 0,/*no remap*/
                     target_comm, recv, send) ;
      }
    } // end if(dflag)

    // if the output variable is in another keyspace, then
    // we need to allocate a local copy and then do the reduce
    storeRepP target_rep_local ;
    if(output_cross_space) {
      // initialize the unit value by using the proper thaw method
#ifdef DYNAMICSCHEDULING
      target_rep_local = target_rep_in_facts->thaw(target_var_exist) ;
#else
      cerr << "dynamic schduling feature not enabled" << endl ;
      Abort() ;
#endif
    } else {
      target_rep_local = target_rep_in_facts ;
      // again, we don't want to allocate the targets
      // because here it might be just portion of its total domain
    }
    rp->set_store(target_var, target_rep_local) ;
#ifdef DYNAMIC_TIMING
    ta_output_oh.addTime(sw_output_oh.stop(), 1) ;

    sw_compute.start() ;
#endif
    // ok, perform the compute
    rp->compute(sequence(context)) ;
#ifdef DYNAMIC_TIMING
    ta_compute.addTime(sw_compute.stop(), 1) ;

    sw_push.start() ;
#endif
    // then reduction
    if(output_cross_space) {
      reduce_store(target_rep_local, 0, /*no remap*/
                   target_rep_in_facts, target_unpack,
                   send, recv, rp->get_joiner(), target_comm) ;
    } else {
      reduce_store(target_rep_local, 0, /*no remap*/
                   target_rep_local, 0, /*no remap*/
                   send, recv, rp->get_joiner(), target_comm) ;
    }
#ifdef DYNAMIC_TIMING
    ta_push.addTime(sw_push.stop(), 1) ;
#endif

    if(output_cross_space) {
      target_rep_local->allocate(EMPTY) ;
    }

    if(dflag) dflag = false ;

    timer.addTime(s.stop(),0) ;
  }

  void execute_dynamic_applyrule::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << apply_tag << " over [dynamic] sequence " << endl ;
  }

  void execute_dynamic_applyrule::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "rule: " << apply_tag ;
    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }
  
  // dynamic apply compiler
  void
  dynamic_apply_compiler::
  set_var_existence(fact_db& facts, sched_db& scheds) {
    // since the dynamic_apply_compiler will
    // evalulate the rules at the runtime, we'll
    // just set everything (existence and request) to
    // be EMPTY. Note: this is also important to interface
    // correctly with the current allocation code and the
    // barrier and reduction compiler. otherwise they'll
    // do some extra job that we don't want.
    variableSet targets = apply.targets() ;
    for(variableSet::const_iterator vi=targets.begin();
        vi!=targets.end();++vi)
      scheds.set_existential_info(*vi, apply, EMPTY) ;
  }

  void
  dynamic_apply_compiler::
  process_var_requests(fact_db& facts, sched_db& scheds) {
    variableSet targets = apply.targets() ;
    for(variableSet::const_iterator vi=targets.begin();
        vi!=targets.end();++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;
    variableSet sources = apply.sources() ;
    for(variableSet::const_iterator vi=sources.begin();
        vi!=sources.end();++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;
  }

  executeP
  dynamic_apply_compiler::create_execution_schedule(fact_db& facts,
                                                    sched_db& scheds) {
    // determine if the target is a param or a store
    variable target = *(apply.targets().begin()) ;
    storeRepP target_rep = facts.get_variable(target) ;
    FATAL(target_rep == 0) ;
    if(isPARAMETER(target_rep)) {
      executeP execute =
        new execute_dynamic_applyrule_param(apply,unit_tag,
                                            space,facts,scheds) ;
      return execute ;
    } else {
      executeP execute =
        new execute_dynamic_applyrule(apply,unit_tag,space,facts,scheds) ;

      return execute ;
    }
  }

  // execute_dynamic_applyrule_param module implementation
  CPTR<joiner>
  execute_dynamic_applyrule_param::apply_join_op = 0 ;
  
  execute_dynamic_applyrule_param::
  execute_dynamic_applyrule_param(rule a, rule u, KeySpaceP kp,
                                  fact_db& facts, sched_db& scheds) {    
    apply_tag = a ;
    unit_tag = u ;
    rp = apply_tag.get_rule_implP() ;
    space = kp ;

    large_context_size = 2 * facts.get_init_ptn(0)[0].size() ;

    dflag = true ;
    // register itself to the keyspace
    space->register_dcontrol_rule(apply_tag, this) ;
    
    FATAL(apply_tag.targets().size() != 1) ;

    set<vmap_info>::const_iterator vmsi ;
    // figure out the expand chain
    // and also initialize rule storeReps

    // build the input chain first
    for(vmsi=apply_tag.get_info().desc.sources.begin();
        vmsi!=apply_tag.get_info().desc.sources.end();++vmsi) {
      collect_expand_chain(*vmsi, "source", rp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, rp, space,
                              facts, scheds, input_ne_source) ;
    } // end of source

    for(vmsi=apply_tag.get_info().desc.constraints.begin();
        vmsi!=apply_tag.get_info().desc.constraints.end();++vmsi) {
      collect_expand_chain(*vmsi, "constraint", rp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, rp, space,
                              facts, scheds, input_ne_constraint) ;
    } // end of constraints

    // Make sure there is no mapping in the output
    vmsi=apply_tag.get_info().desc.targets.begin() ;
    if(!vmsi->mapping.empty()) {
      if(space->is_distributed()) {
        debugout << "Dynamic global reduction cannot have maps in target!"
                 << ", Offending rule: " << apply_tag << endl ;
        if(space->get_comm_rank() == 0) {
          cerr << "Dynamic global reduction cannot have maps in target!"
               << ", Offending rule: " << apply_tag << endl ;
        }
        Loci::Abort() ;
      } else {
        cerr << "Dynamic global reduction cannot have maps in target!"
             << ", Offending rule: " << apply_tag << endl ;
        Loci::Abort() ;
      }
    }

    // since this is a param apply rule, we don't need output chain
    target_var = *(apply_tag.targets().begin()) ;
    target_rep = facts.get_variable(target_var) ;
    FATAL(target_rep == 0) ;
    rp->set_store(target_var, target_rep) ;

    // finally we build a MPI Op
    MPI_Op_create(&(execute_dynamic_applyrule_param::create_user_function),
                  0, &mpi_reduction_op) ;
  }

  execute_dynamic_applyrule_param::
  ~execute_dynamic_applyrule_param() {
    MPI_Op_free(&mpi_reduction_op) ;
  }

  void execute_dynamic_applyrule_param::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;

    if(dflag) {
#ifdef DYNAMIC_TIMING
      sw_expand.start() ;
#endif
      // first we need to communicate input_chains
      for(vector<ExpandChain>::iterator chain=input_chains.begin();
          chain!=input_chains.end();++chain) {
        execute_expand_chain(*chain, 0) ;
      }     // end all chain expansion
#ifdef DYNAMIC_TIMING
      ta_expand.addTime(sw_expand.stop(), 1) ;
      
      sw_context.start() ;
#endif
      
      // then we need to compute the execution context
      context =
        dynamic_rule_context(apply_tag,space,
                             facts,scheds,input_chains,
                             input_ne_source, input_ne_constraint);
#ifdef DYNAMIC_TIMING
      ta_context.addTime(sw_context.stop(), 1) ;
#endif
    } // end if(dflag)
    
    // // give a warning if context is very large
    // size_t context_size = context.size() ;
    // if(context_size > large_context_size) {
    //   if(space->get_comm_rank() == 0)
    //     cerr << "Warning: dynamic rule: " << apply_tag
    //          << " executing over a very large context, potentially"
    //          << " fatal to memory allocation" << endl ;
    // }

    // we then go ahead to expand the output chain
    const entitySet& target_var_exist = context ;
    target_rep->guarantee_domain(target_var_exist) ;

#ifdef DYNAMIC_TIMING
    sw_compute.start() ;
#endif
    // ok, perform the compute
    rp->compute(sequence(context)) ;
#ifdef DYNAMIC_TIMING
    ta_compute.addTime(sw_compute.stop(), 1) ;
    sw_param_push.start() ;
#endif
    // then reduction

    // first set the static joiner op
    apply_join_op = rp->get_joiner() ;

    entitySet e ;
    sequence seq ;
    int size = target_rep->pack_size(e) ;

    unsigned char* send_ptr = new unsigned char[size] ;
    int position = 0 ;
#ifdef DYNAMIC_TIMING
    sw_param_pack.start() ;
#endif

    target_rep->pack(send_ptr, position, size, e) ;

#ifdef DYNAMIC_TIMING
    ta_param_pack.addTime(sw_param_pack.stop(), 1) ;
#endif

    unsigned char* result_ptr = new unsigned char[size] ;

#ifdef DYNAMIC_TIMING
    sw_param_reduce.start() ;
#endif
    
    MPI_Allreduce(send_ptr, result_ptr, size,
                  MPI_PACKED, mpi_reduction_op, space->get_mpi_comm()) ;

#ifdef DYNAMIC_TIMING
    ++ta_param_reduce_num ;
    ta_param_reduce.addTime(sw_param_reduce.stop(), 1) ;
    sw_param_unpack.start() ;
#endif

    position = 0 ;
    target_rep->unpack(result_ptr, position, size, seq) ;

    delete[] result_ptr ;
    delete[] send_ptr ;
#ifdef DYNAMIC_TIMING
    ta_param_unpack.addTime(sw_param_unpack.stop(), 1) ;
    ta_param_push.addTime(sw_param_push.stop(), 1) ;
#endif
    
    if(dflag) dflag = false ;
    
    timer.addTime(s.stop(),0) ;
  }

  void execute_dynamic_applyrule_param::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << apply_tag << " over [dynamic] sequence " << endl ;
  }

  void execute_dynamic_applyrule_param::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "rule: " << apply_tag ;
    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }

  void execute_dynamic_applyrule_param::
  create_user_function(void* send_ptr, void* result_ptr,
                       int* size, MPI_Datatype* dptr) {
    entitySet e ;
    sequence seq ;
    int unpack_send_position = 0 ;
    int unpack_result_position = 0 ;
    int pack_result_position = 0 ;
    storeRepP sp, tp ;
    sp = apply_join_op->getTargetRep() ;
    tp = apply_join_op->getTargetRep() ;
    sp->unpack(send_ptr, unpack_send_position, *size, seq) ;
    tp->unpack(result_ptr, unpack_result_position, *size, seq) ;
    apply_join_op->SetArgs(tp, sp) ;
    apply_join_op->Join(seq) ;
    tp->pack(result_ptr, pack_result_position, *size, e) ;
  }
  
  // dclone_invalidate_compiler
  void
  dclone_invalidate_compiler::set_var_existence(fact_db& facts,
                                                sched_db& scheds) {
    // do nothing
  }

  void
  dclone_invalidate_compiler::process_var_requests(fact_db& facts,
                                                   sched_db& scheds) {
    // do nothing
  }

  executeP
  dclone_invalidate_compiler::create_execution_schedule(fact_db& facts,
                                                        sched_db& scheds) {
    executeP execute =
      new execute_dclone_invalidator(var, var_unique, self_clone,
                                     shadow_clone, facts, scheds) ;
    return execute ;
  }

  // execute_dclone_invalidator module code
  execute_dclone_invalidator::
  execute_dclone_invalidator(variable v, variable vu,
                             KeySpaceP sc,
                             const std::vector<KeySpaceP>& sac,
                             fact_db& facts, sched_db& scheds) {
    var = v ;
    var_unique = vu ;
    self_clone = sc ;
    shadow_clone = sac ;
  }

  void execute_dclone_invalidator::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;
#ifdef DYNAMIC_TIMING
    sw_invalidate.start() ;
#endif
    if(self_clone != 0) {
      // get the storeRepP from fact_db
      storeRepP s = facts.get_variable(var) ;
      // get the ptn info from keyspace
      entitySet local_ptn =
        self_clone->get_key_ptn()[self_clone->get_comm_rank()] ;
      entitySet dom = s->domain() ;
      entitySet valid = dom & local_ptn ;
      // invalidate others
      s->invalidate(valid) ;
    }
    for(vector<KeySpaceP>::const_iterator vi=shadow_clone.begin();
        vi!=shadow_clone.end();++vi) {
      // invalidate all shadow clones
      // get the shadow clone
      storeRepP s = (*vi)->get_shadow(var) ;
      s->invalidate(EMPTY) ;
    }
#ifdef DYNAMIC_TIMING
    ta_invalidate.addTime(sw_invalidate.stop(), 1) ;
#endif
    timer.addTime(s.stop(), 0) ;
  }

  void execute_dclone_invalidator::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << "Invalidate dynamic clone of " << var << endl ;
  }

  void execute_dclone_invalidator::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "dynamic clone invalidator: " << var ;
    data_collector.accumulateTime(timer,EXEC_CONTROL,oss.str()) ;
  }
  
  // keyspace_dist_compiler
  void
  keyspace_dist_compiler::set_var_existence(fact_db& facts,
                                            sched_db& scheds) {
    // do nothing
  }

  void
  keyspace_dist_compiler::process_var_requests(fact_db& facts,
                                               sched_db& scheds) {
    // do nothing
  }

  executeP
  keyspace_dist_compiler::create_execution_schedule(fact_db& facts,
                                                    sched_db& scheds) {
    // get the keyspace pointers
    vector<KeySpaceP> spaces ;
    for(vector<string>::const_iterator vi=space_names.begin();
        vi!=space_names.end();++vi) {
      map<string,KeySpaceP>::const_iterator ki ;
      ki = facts.keyspace.find(*vi) ;
      if(ki == facts.keyspace.end()) {
        cerr << "Error: keyspace: " << *vi
             << " does not exist, error from "
             << "keyspace_dist_compiler" << endl ;
        Loci::Abort() ;
      }
      spaces.push_back(ki->second) ;
    }
    
    executeP execute =
      new execute_keyspace_dist(spaces,facts,scheds) ;

    return execute ;
  }

  execute_keyspace_dist::
  execute_keyspace_dist(const std::vector<KeySpaceP>& s,
                        fact_db& facts, sched_db& scheds):spaces(s) {
    FATAL(spaces.empty()) ;
    ostringstream oss ;
    vector<KeySpaceP>::const_iterator vi = spaces.begin() ;
    oss << (*vi)->get_name() ;
    for(++vi;vi!=spaces.end();++vi)
      oss << ", " << (*vi)->get_name() ;
    space_names = oss.str() ;

    for(vector<KeySpaceP>::const_iterator ki=spaces.begin();
        ki!=spaces.end();++ki) {
      variableSet t = (*ki)->get_tunnels() ;
      variableSet allt ;
      for(variableSet::const_iterator vi=t.begin();
          vi!=t.end();++vi) {
        allt += *vi ;
        allt += scheds.try_get_aliases(*vi) ;
        allt += scheds.try_get_antialiases(*vi) ;
        allt += scheds.try_get_synonyms(*vi) ;
        allt += scheds.try_get_rotations(*vi) ;
      }
      t = remove_synonym(allt,facts) ;
      tunnels.push_back(t) ;

      variableSet vars = (*ki)->get_control_vars() ;
      vars = remove_synonym(vars,facts) ;
      space_vars.push_back(vars) ;
    }
  }
  
  void execute_keyspace_dist::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;
#ifdef DYNAMIC_TIMING
    ++ta_dist_number ;
    // Code Tag
    MPI_Barrier(MPI_COMM_WORLD) ;
    ////////
    sw_dist.start() ;
#endif
    for(size_t i=0;i<spaces.size();++i) {
      KeySpaceP sp = spaces[i] ;
      const variableSet& tn = tunnels[i] ;
      const variableSet& vars = space_vars[i] ;
      // first redistribute the keyspace keys
#ifdef DYNAMIC_TIMING
      sw_key_dist.start() ;
#endif
      bool dist = sp->update_key_ptn() ;
#ifdef DYNAMIC_TIMING
      ta_key_dist.addTime(sw_key_dist.stop(), 1) ;
#endif
      if(!dist)
        continue ;
      // get the new key partition
      const vector<entitySet>& ptn = sp->get_key_ptn() ;
#ifdef DYNAMIC_TIMING
      sw_dist_renumber.start() ;
#endif
      // then we renumber the partition
      vector<entitySet> new_ptn ;
      dMap remap ;
      KeyManagerP key_manager = facts.get_key_manager() ;
      key_manager->compact(ptn, remap, new_ptn) ;
#ifdef DYNAMIC_TIMING
      ta_dist_renumber.addTime(sw_dist_renumber.stop(), 1) ;
#endif
      // then redistribute all facts using the renumber
      MPI_Comm comm = sp->get_mpi_comm() ;
      // redistribute vars
      for(variableSet::const_iterator vi=vars.begin();
          vi!=vars.end();++vi) {
        storeRepP s = facts.get_variable(*vi) ;
        FATAL(s == 0) ;
        if(tn.inSet(*vi))
          facts.update_fact(*vi, s->redistribute_omd(ptn,remap,comm)) ;
        else
          facts.update_fact(*vi, s->redistribute(ptn,remap,comm)) ;
      }
      // update the key_ptn in the space
      sp->set_key_ptn(new_ptn) ;
    }
#ifdef DYNAMIC_TIMING
    ta_dist.addTime(sw_dist.stop(), 1) ;
#endif
    timer.addTime(s.stop(), 0) ;
  }

  void execute_keyspace_dist::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << "keyspace distribution: (" << space_names << ")" << endl ;
  }

  void execute_keyspace_dist::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "keyspace distribution: (" << space_names << ")" ;
    data_collector.accumulateTime(timer,EXEC_COMMUNICATION,oss.str()) ;
  }
#endif

#ifdef DYNAMICSCHEDULING
  // execute_keyspace_init module
  execute_init_keyspace::
  execute_init_keyspace(fact_db& facts, sched_db& scheds) {
    map<string,KeySpaceP>::iterator ki ;
    for(ki=facts.keyspace.begin();ki!=facts.keyspace.end();++ki) {
      KeySpaceP space = ki->second ;
      const variableSet& critical_vars = space->get_critical_vars() ;
      for(variableSet::const_iterator vi=critical_vars.begin();
          vi!=critical_vars.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        if(srp == 0) {
          cerr << "Error: keyspace: " << ki->first << " critical"
               << " structure: " << *vi << " has no setup in fact_db"
               << endl ;
          Loci::Abort() ;
        } else
          space->setup_csn_rep(*vi, srp) ;
      }
      // then we'll expand all the controlling variables
      // to include all possible aliases, synonyms, rotations etc
      // since at this point, we know all these information
      variableSet control_vars = space->get_control_vars() ;
      variableSet expand_control_vars ;
      for(variableSet::const_iterator vi=control_vars.begin();
          vi!=control_vars.end();++vi) {
        expand_control_vars += *vi ;
        expand_control_vars += scheds.get_aliases(*vi) ;
        expand_control_vars += scheds.get_antialiases(*vi) ;
        expand_control_vars += scheds.get_synonyms(*vi) ;
        expand_control_vars += scheds.get_rotations(*vi) ;
      }
      space->add_control_vars(expand_control_vars) ;

      // integraty check for this keyspace
      string err ;
      if(!space->integrity_check(err)) {
        cerr << "Error: " << err << endl ;
        Loci::Abort() ;
      }
    } // finish all keyspace init
    // finally init the key manager in the facts
    facts.init_key_manager() ;
  }

  void execute_init_keyspace::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;
    map<string,KeySpaceP>::iterator ki ;
    for(ki=facts.keyspace.begin();ki!=facts.keyspace.end();++ki) {
      KeySpaceP space = ki->second ;
      space->set_dcontrol() ;
    }
    timer.addTime(s.stop(), 0) ;
  }
  
  void execute_init_keyspace::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << "Initialize all keyspaces" << endl ;
  }

  void execute_init_keyspace::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "Initialize all keyspaces" ;
    data_collector.accumulateTime(timer,EXEC_CONTROL,oss.str()) ;
  }
#endif
  
  // blackbox_compiler code
  void
  blackbox_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    //    existential_blackboxrule_analysis(impl, facts, scheds) ;
    // set UNIVERSE existence for all targets
     variableSet targets = impl.targets() ;
     for(variableSet::const_iterator vi=targets.begin();vi!=targets.end();++vi)
       scheds.set_existential_info(*vi, impl, ~EMPTY) ;
  }

  void blackbox_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    //    exec_seq = process_blackboxrule_requests(impl, facts, scheds) ;

    //everyone will need to request for their existence
    variableSet targets = impl.targets() ;
    for(variableSet::const_iterator vi = targets.begin(); vi != targets.end(); ++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;

    variableSet sources = impl.sources() ;
    for(variableSet::const_iterator vi=sources.begin(); vi != sources.end(); ++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;
    entitySet exec_seq = ~EMPTY ;
    scheds.update_exec_seq(impl, exec_seq);
  }

  executeP blackbox_compiler::
  create_execution_schedule(fact_db& facts, sched_db& scheds) {
    entitySet exec_seq = scheds.get_exec_seq(impl);
    executeP execute = new execute_rule(impl,
                                        sequence(exec_seq), facts, scheds);
    return execute ;
  }

  // superRule_compiler code
  void
  superRule_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    CPTR<super_rule> rp(impl.get_rule_implP()) ;
    rp->process_existential(impl,facts,scheds) ;
  }

  void superRule_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    CPTR<super_rule> rp(impl.get_rule_implP()) ;
    rp->process_requests(impl,facts,scheds) ;
  }

  executeP superRule_compiler::create_execution_schedule(fact_db& facts, sched_db& scheds) {
    executeP execute = new execute_rule(impl, ~EMPTY, facts, scheds);
    return execute;
  }

#ifdef DYNAMICSCHEDULING

  // here comes the insertion and deletion rule impls
  executeP insertion_rule_compiler::
  create_execution_schedule(fact_db& facts, sched_db& scheds) {
    executeP execute = new execute_insertion(rule_tag, facts, scheds) ;

    return execute ;
  }

  executeP deletion_rule_compiler::
  create_execution_schedule(fact_db& facts, sched_db& scheds) {
    executeP execute = new execute_key_destruction(rule_tag, facts, scheds) ;

    return execute ;
  }

  execute_insertion::
  execute_insertion(const rule& r, fact_db& facts, sched_db& scheds) {
    rule_tag = r ;
    rule_implP tagp = r.get_rule_implP() ;
    // first check to make sure that no output targets
    // have maps involved. because for an insertion rule,
    // that does not make sense.
    set<vmap_info>::const_iterator vmsi ;
    for(vmsi=rule_tag.get_info().desc.targets.begin();
        vmsi!=rule_tag.get_info().desc.targets.end();++vmsi) {
      if(!vmsi->mapping.empty()) {
        cerr << "Error: targets of insertion rule cannot "
             << "have mappings involved! Offending rule: "
             << rule_tag << endl ;
        Loci::Abort() ;
      }
    }
    // we will need to check to see if the rule is in
    // a dynamic keyspace
    string space_name = tagp->get_keyspace_tag() ;
    map<string,KeySpaceP>::const_iterator
      mi = facts.keyspace.find(space_name) ;
    if(mi == facts.keyspace.end()) {
      cerr << "Error: insertion rule: " << r << " exists in"
           << " unknown keyspace: " << space_name << endl ;
      Loci::Abort() ;
    }
    space = mi->second ;
    if(space->get_dynamism() != DYNAMIC) {
      cerr << "Error: only dynamic keyspace rule is allowed"
           << " to have insertion rule. Offending rule: "
           << r << " in keyspace: " << space_name << endl ;
      Loci::Abort() ;
    }
    rp = dynamic_cast<insertion_rule_interface*>(&(*tagp)) ;
    if(rp == 0) {
      cerr << "Error: execute_insertion module does not have an"
           << " insertion rule associated!" << endl ;
      Loci::Abort() ;
    }
    key_manager = facts.get_key_manager() ;
    rp->initialize(facts) ;
    // then initialize the key manager in the rule
    rp->set_key_manager(key_manager) ;
    // there is no need to set up expand chains for insertion rule
  }

  void execute_insertion::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;

#ifdef DYNAMIC_TIMING
    sw_insertion.start() ;
#endif
    // first call the compute to compelete the insertion
    rp->compute(EMPTY);
    // then obtain the keys and modify the keyspace accordingly
    entitySet new_keys = rp->get_keys_inserted() ;
#ifdef DYNAMIC_TIMING
    sw_keyinsert.start() ;
#endif
    space->add_keys(new_keys) ;
#ifdef DYNAMIC_TIMING
    ta_keyinsert.addTime(sw_keyinsert.stop(), 1) ;
    ta_insertion.addTime(sw_insertion.stop(), 1) ;
#endif

    timer.addTime(s.stop(),0) ;
  }

  void execute_insertion::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << "Insert data into: " << rule_tag.targets() << endl ;
  }

  void execute_insertion::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "Insert data into: " << rule_tag.targets() ;
    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }

  execute_key_destruction::
  execute_key_destruction(const rule& r,
                          fact_db& facts, sched_db& scheds) {
    rule_tag = r ;
    rule_implP tagp = r.get_rule_implP() ;

    dflag = true ;
    // first check to make sure that no output targets
    // have maps involved. because for an deletion rule,
    // that does not make sense.
    set<vmap_info>::const_iterator vmsi ;

    if(rule_tag.get_info().desc.targets.size() != 1) {
      cerr << "Error: key destruction rule can only have one"
           << " target of parameter<bool>! Offending rule: "
           << rule_tag << endl ;
      Loci::Abort() ;
    }
    vmsi = rule_tag.get_info().desc.targets.begin() ;
    if(vmsi->var.size() != 1) {
      cerr << "Error: key destruction rule can only have one"
           << " target of parameter<bool>! Offending rule: "
           << rule_tag << endl ;
      Loci::Abort() ;
    } else {
      target = *(vmsi->var.begin()) ;
      target_rep = tagp->get_store(target) ;
      if(!isPARAMETER(target_rep)) {
        cerr << "Error: key destruction rule can only have one"
             << " target of parameter<bool>! Offending rule: "
             << rule_tag << endl ;
        Loci::Abort() ;
      } else {
        // initialize the rule's store
        target_rep = facts.get_variable(target) ;
        FATAL(target_rep == 0) ;
        tagp->set_store(target, target_rep) ;
      }
    }
    // check for mapping
    if(!vmsi->mapping.empty()) {
      cerr << "Error: targets of key destruction rule cannot "
           << "have mappings involved! Offending rule: "
           << rule_tag << endl ;
      Loci::Abort() ;
    }

    // we will need to check to see if the rule is in
    // a dynamic keyspace
    string space_name = tagp->get_keyspace_tag() ;
    map<string,KeySpaceP>::const_iterator
      mi = facts.keyspace.find(space_name) ;
    if(mi == facts.keyspace.end()) {
      cerr << "Error: key destruction rule: " << r << " exists in"
           << " unknown keyspace: " << space_name << endl ;
      Loci::Abort() ;
    }
    space = mi->second ;
    if(space->get_dynamism() != DYNAMIC) {
      cerr << "Error: only dynamic keyspace rule is allowed"
           << " to have key destruction rule. Offending rule: "
           << r << " in keyspace: " << space_name << endl ;
      Loci::Abort() ;
    }
    // register the module to the keyspace
    space->register_dcontrol_rule(rule_tag,this) ;
    rp = dynamic_cast<deletion_rule*>(&(*tagp)) ;
    if(rp == 0) {
      cerr << "Error: execute_key_destruction module does not have a"
           << " deletion rule associated!" << endl ;
      Loci::Abort() ;
    }
    key_manager = facts.get_key_manager() ;
    // then initialize the key manager in the rule
    rp->set_key_manager(key_manager) ;

    large_context_size = 2 * facts.get_init_ptn(0)[0].size() ;

    // then build input chains, if any

    // figure out the expand chain
    // and also initialize rule storeReps
    for(vmsi=rule_tag.get_info().desc.sources.begin();
        vmsi!=rule_tag.get_info().desc.sources.end();++vmsi) {
      collect_expand_chain(*vmsi, "source", tagp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, tagp, space,
                              facts, scheds, input_ne_source) ;
    } // end of source

    for(vmsi=rule_tag.get_info().desc.constraints.begin();
        vmsi!=rule_tag.get_info().desc.constraints.end();++vmsi) {
      collect_expand_chain(*vmsi, "constraint", tagp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, tagp, space,
                              facts, scheds, input_ne_constraint) ;
    } // end of constraints
    // since we only have one parameter target and have
    // already processed it, we don't need to do anything here
    // execpt for initializing the rep in the rule, which we
    // already did in the previous step    
    
    // finally we need to build a vector of stores that this
    // keyspace is controlling, we need to erase the keys in all
    // of them
    variableSet control_vars = space->get_control_vars() ;
    control_vars = remove_synonym(control_vars, facts) ;
    for(variableSet::const_iterator vi=control_vars.begin();
        vi!=control_vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      controlled_stores.push_back(srp) ;
    }
  }

  void execute_key_destruction::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;

    if(dflag) {
#ifdef DYNAMIC_TIMING
      sw_expand.start() ;
#endif
      // first we need to communicate input_chains
      for(vector<ExpandChain>::iterator chain=input_chains.begin();
          chain!=input_chains.end();++chain) {
        execute_expand_chain(*chain, 0) ;
      }     // end all chain expansion
      
#ifdef DYNAMIC_TIMING
      ta_expand.addTime(sw_expand.stop(), 1) ;
      
      sw_context.start() ;
#endif
      
      // then we need to compute the execution context
      context =
        dynamic_rule_context(rule_tag, space,
                             facts, scheds, input_chains,
                             input_ne_source, input_ne_constraint) ;
#ifdef DYNAMIC_TIMING
      ta_context.addTime(sw_context.stop(), 1) ;
      
      sw_erase.start() ;
#endif
    }

    // // give a warning if context is very large
    // if(context.size() > large_context_size) {
    //   if(space->get_comm_rank() == 0)
    //     cerr << "Warning: dynamic rule: " << rule_tag
    //          << " executing over a very large context, potentially"
    //          << " fatal to memory allocation" << endl ;
    // }

    // allocate the target on all the local keys
    target_rep->allocate(space->get_keys_local()) ;

    // then call the compute method
    rp->compute(sequence(context)) ;

    // then we need to destroy the keys
    entitySet keys_to_destroy = rp->get_keys_deleted() ;
    // first all remove keys from all the controlling stores
    for(vector<storeRepP>::iterator vi=controlled_stores.begin();
        vi!=controlled_stores.end();++vi) {
      (*vi)->erase(keys_to_destroy) ;
    }
    // finally remove keys from keyspace
#ifdef DYNAMIC_TIMING
    sw_keyremoval.start() ;
#endif
    space->remove_keys(keys_to_destroy) ;
#ifdef DYNAMIC_TIMING
    ta_keyremoval.addTime(sw_keyremoval.stop(), 1) ;
#endif

    param<bool> target_param(target_rep) ;
    *target_param = true ;

#ifdef DYNAMIC_TIMING
    ta_erase.addTime(sw_erase.stop(), 1) ;
#endif
    if(dflag) dflag = false ;
    
    timer.addTime(s.stop(),0) ;
  }

  void execute_key_destruction::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << "Destroy keys according to: " << rule_tag.sources() << endl ;
  }

  void execute_key_destruction::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "Destroy keys according to: " << rule_tag.sources() ;
    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }

  executeP erase_rule_compiler::
  create_execution_schedule(fact_db& facts, sched_db& scheds) {
    executeP execute = new execute_erase(rule_tag, facts, scheds) ;

    return execute ;
  }

  execute_erase::
  execute_erase(const rule& r,
                fact_db& facts, sched_db& scheds) {
    rule_tag = r ;
    rule_implP tagp = r.get_rule_implP() ;

    dflag = true ;
    // we will need to check to see if the rule is in
    // a dynamic keyspace
    string space_name = tagp->get_keyspace_tag() ;
    map<string,KeySpaceP>::const_iterator
      mi = facts.keyspace.find(space_name) ;
    if(mi == facts.keyspace.end()) {
      cerr << "Error: erase rule: " << r << " exists in"
           << " unknown keyspace: " << space_name << endl ;
      Loci::Abort() ;
    }
    space = mi->second ;
    if(space->get_dynamism() != DYNAMIC) {
      cerr << "Error: only dynamic keyspace rule is allowed"
           << " to have erase rule. Offending rule: "
           << r << " in keyspace: " << space_name << endl ;
      Loci::Abort() ;
    }
    // register the module to the keyspace
    space->register_dcontrol_rule(rule_tag,this) ;
    rp = dynamic_cast<erase_rule*>(&(*tagp)) ;
    if(rp == 0) {
      cerr << "Error: execute_erase module does not have an"
           << " erase rule associated!" << endl ;
      Loci::Abort() ;
    }
    large_context_size = 2 * facts.get_init_ptn(0)[0].size() ;

    set<vmap_info>::const_iterator vmsi ;
    // figure out the expand chain
    // and also initialize rule storeReps
    for(vmsi=rule_tag.get_info().desc.sources.begin();
        vmsi!=rule_tag.get_info().desc.sources.end();++vmsi) {
      collect_expand_chain(*vmsi, "source", tagp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, tagp, space,
                              facts, scheds, input_ne_source) ;
    } // end of source

    for(vmsi=rule_tag.get_info().desc.constraints.begin();
        vmsi!=rule_tag.get_info().desc.constraints.end();++vmsi) {
      collect_expand_chain(*vmsi, "constraint", tagp, space,
                           facts, scheds, input_chains, 0, 0, 0, 0, 0) ;
      collect_nonexpand_chain(*vmsi, tagp, space,
                              facts, scheds, input_ne_constraint) ;
    } // end of constraints
    
    // check to make sure no mapping is involved in the targets
    std::vector<ExpandChain> tmp_chain ;
    for(vmsi=rule_tag.get_info().desc.targets.begin();
        vmsi!=rule_tag.get_info().desc.targets.end();++vmsi) {
      collect_expand_chain(*vmsi, "target", tagp, space,
                           facts, scheds, tmp_chain, 0, 0, 0, 0, 0) ;
      if(!tmp_chain.empty()) {
        cerr << "Error: erase rule cannot have mappings in the"
             << " output. Offending rule: " << rule_tag << endl ;
        Loci::Abort() ;
      }
      collect_nonexpand_chain(*vmsi, tagp, space,
                              facts, scheds, outputs) ;
    }
    
    targets = rule_tag.targets() ;
  }

  void execute_erase::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;

    if(dflag) {
#ifdef DYNAMIC_TIMING
      sw_expand.start() ;
#endif
      // first we need to communicate input_chains
      for(vector<ExpandChain>::iterator chain=input_chains.begin();
          chain!=input_chains.end();++chain) {
        execute_expand_chain(*chain, 0) ;
      }     // end all chain expansion
      
#ifdef DYNAMIC_TIMING
      ta_expand.addTime(sw_expand.stop(), 1) ;
      
      sw_context.start() ;
#endif
      
      // then we need to compute the execution context
      context =
        dynamic_rule_context(rule_tag, space,
                             facts, scheds, input_chains,
                             input_ne_source, input_ne_constraint) ;
#ifdef DYNAMIC_TIMING
      ta_context.addTime(sw_context.stop(), 1) ;
      
      sw_erase.start() ;
#endif
    }

    // targets don't need to be allocated

    // then call the compute method
    rp->compute(sequence(context)) ;

#ifdef DYNAMIC_TIMING
    sw_record_erase.start() ;
#endif
    // then we need to erase the records from the targets
    entitySet record_to_erase = rp->get_erased_record() ;
    // first all remove keys from all the controlling stores
    for(vector<NonExpandUnit>::iterator unit=outputs.begin();
        unit!=outputs.end();++unit) {
      unit->var_rep->erase(record_to_erase) ;
    }
#ifdef DYNAMIC_TIMING
    ta_record_erase.addTime(sw_record_erase.stop(), 1) ;
#endif

    if(dflag) dflag = false ;

    timer.addTime(s.stop(),0) ;
  }

  void execute_erase::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << "Erase records from: " << targets << endl ;
  }

  void execute_erase::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "Erase records from: " << targets ;
    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
  }

  // dcontrol_reset_compiler
  void
  dcontrol_reset_compiler::set_var_existence(fact_db& facts,
                                             sched_db& scheds) {
    // do nothing
  }

  void
  dcontrol_reset_compiler::process_var_requests(fact_db& facts,
                                                sched_db& scheds) {
    // do nothing
  }

  executeP
  dcontrol_reset_compiler::create_execution_schedule(fact_db& facts,
                                                     sched_db& scheds) {
    vector<KeySpaceP> register_kp ;
    for(set<string>::const_iterator si=register_keyspaces.begin();
        si!=register_keyspaces.end();++si) {
      KeySpaceP kp = facts.get_keyspace(*si) ;
      if(kp == 0) {
        cerr << "Error: keyspace: " << *si
             << " does not exist, error from dcontrol_reset_compiler"
             << endl ;
        Loci::Abort() ;
      }
      register_kp.push_back(kp) ;
    }
    
    executeP execute =
      new execute_dcontrol_reset(var, var_unique,
                                 register_kp, facts, scheds) ;
    return execute ;
  }

  // execute_dcontrol_reset
  execute_dcontrol_reset::
  execute_dcontrol_reset(variable v, variable vu,
                         const std::vector<KeySpaceP>& ks,
                         fact_db& facts, sched_db& scheds)
    :var(v), var_unique(vu), register_keyspaces(ks) {
    // register itself to the keyspaces
    for(vector<KeySpaceP>::iterator ki=register_keyspaces.begin();
        ki!=register_keyspaces.end();++ki) {
      (*ki)->register_dcontrol_var(v, this) ;
    }
  }

  void execute_dcontrol_reset::
  execute(fact_db& facts, sched_db& scheds) {
    stopWatch s ;
    s.start() ;
#ifdef DYNAMIC_TIMING
    sw_dctrl.start() ;
#endif
    for(vector<execute_druleP>::iterator rp=drules.begin();
        rp!=drules.end();++rp) {
      (*rp)->reset_ctrl() ;
    }
#ifdef DYNAMIC_TIMING
    ta_dctrl.addTime(sw_dctrl.stop(), 1) ;
#endif
    timer.addTime(s.stop(), 0) ;
  }

  void execute_dcontrol_reset::
  Print(std::ostream& s) const {
    printIndent(s) ;
    s << "Reset Drule control flag by: " << var << endl ;
  }

  void execute_dcontrol_reset::
  dataCollate(collectData& data_collector) const {
    ostringstream oss ;
    oss << "Reset Drule control flag by: " << var ;
    data_collector.accumulateTime(timer,EXEC_CONTROL,oss.str()) ;
  }
#endif  

  ///////////////////////////////////////////////
  // Lets set up some common super rule functions
  namespace {
    const char *systemVarDocs[] = {
      "SYSTEM\000NOT(X)\000param<bool>\000Exists over complement of X",
      "SYSTEM\000OR(X,Y)\000param<bool>\000Exists over union of X and Y",
      "SYSTEM\000OR(X,Y,Z)\000param<bool>\000Exists over union of X, Y and Z",
      "SYSTEM\000OR(W,X,Y,Z)\000param<bool>\000Exists over union of W, X, Y and Z",
      "SYSTEM\000AND(X,Y)\000param<bool>\000Exists over intersection of X and Y",
      "SYSTEM\000AND(X,Y,Z)\000param<bool>\000Exists over intersection of X, Y and Z",
      "SYSTEM\000AND(W,X,Y,Z)\000param<bool>\000Exists over intersection of W, X, Y and Z",
      "\000\000\000"
    } ;
  }
  class NOT_rule : public super_rule {
    param<bool> NOT ;
  public:
    NOT_rule() {
      comments("Internal System Rule.  Use with care to prevent infinite regression!") ;
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("NOT(X)",0) ;
      name_store("NOT(X)",NOT) ;
      constraint("X") ;
      output("NOT(X)") ;
    }
    void compute(const sequence &seq) {
      *NOT = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = ~EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints &= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      //      debugout << "constraints = " << constraints << endl ;
      // Now complement (entities we own)
      constraints = (~constraints) & my_entities ;
      //      debugout << "constraints & my_entities =" << constraints << endl ;

      variableSet output = r.targets() ;
      //      debugout << "setting " << output << endl ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<NOT_rule> register_NOT_rule ;

  class OR_rule : public super_rule {
    param<bool> OR ;
  public:
    OR_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("OR(X,Y)",1) ;
      name_store("OR(X,Y)",OR) ;
      constraint("X,Y") ;
      output("OR(X,Y)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR_rule> register_OR_rule ;

  class OR3_rule : public super_rule {
    param<bool> OR ;
  public:
    OR3_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("OR(X,Y,Z)",2) ;
      name_store("OR(X,Y,Z)",OR) ;
      constraint("X,Y,Z") ;
      output("OR(X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR3_rule> register_OR3_rule ;

  class OR4_rule : public super_rule {
    param<bool> OR ;
  public:
    OR4_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("OR(W,X,Y,Z)",3) ;
      name_store("OR(W,X,Y,Z)",OR) ;
      constraint("W,X,Y,Z") ;
      output("OR(W,X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR4_rule> register_OR4_rule ;

  class AND2_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND2_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("AND(X,Y)",4) ;
      name_store("AND(X,Y)",AND) ;
      constraint("X,Y") ;
      output("AND(X,Y)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND2_rule> register_AND2_rule ;

  class AND3_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND3_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("AND(X,Y,Z)",5) ;
      name_store("AND(X,Y,Z)",AND) ;
      constraint("X,Y,Z") ;
      output("AND(X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND3_rule> register_AND3_rule ;

  class AND4_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND4_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("AND(W,X,Y,Z)",6) ;
      name_store("AND(W,X,Y,Z)",AND) ;
      constraint("W,X,Y,Z") ;
      output("AND(W,X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND4_rule> register_AND4_rule ;
    
  
}

// ... the end ...
