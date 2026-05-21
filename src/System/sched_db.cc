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
#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>

#include <sched_db.h>
#include <fact_db.h>
using std::string ; 
using std::map ;
using std::make_pair ;
using std::vector ;
using std::list ;
using std::sort ;
using std::pair ;
using std::make_pair ;
#include <Tools/parse.h>

using std::cerr ;
using std::endl ;
using std::ofstream ;

//#define VERBOSE

namespace Loci {

  extern ofstream debugout ;
  extern bool rule_has_mapping_in_output(rule r);
  sched_db::sched_db(fact_db &facts) {
    detected_errors = false ;
    init(facts) ;
  }
  
  sched_db::sched_db() {detected_errors = false ;}
  sched_db::~sched_db() {}

  void
  sched_db::init(fact_db &facts) {
    variableSet tmp_all_vars = facts.get_typed_variables() ;
    for(size_t i = 0; i < tmp_all_vars.size(); i++)
      sched_infov.push_back(sched_data()) ;
    int sched_alloc = 0 ;
    for(variableSet::const_iterator vi = tmp_all_vars.begin(); vi != tmp_all_vars.end(); ++vi) {
      sched_info si ; 
      si.sched_info_ref = sched_alloc++ ;
      storeRepP rp = facts.get_variable(*vi) ;
      si.existence = rp->domain() ;
      // si.fact_installed = si.existence ;/*!fact_installed assigned, never used*/
      si.synonyms  += *vi ;
      vmap[*vi] = si ;
      variable v = variable(*vi) ;
      storeRepP st = rp->getRep() ;
      sched_data sd(v, st) ;
      sd.aliases += *vi ;
      sched_infov[si.sched_info_ref] = sd ;
    }

    //    synonyms = facts.synonyms ;
    //    free_set = facts.free_set ;

    all_vars = tmp_all_vars ;
  }
  
  void sched_db::install_sched_info(variable v, sched_info info) {
    vmap_type::iterator mi = vmap.find(v) ;
    if(mi != vmap.end()) {
      cerr << "error:  reinstalling fact in database for variable "
	   << v << endl ;
      abort() ;
    }
    info.synonyms += v ;
    vmap[v] = info ;
    sched_infov[info.sched_info_ref].aliases += v ;
    all_vars += v ;
  }
  /*! since free_set is never assigned, else block will never be executed, commented out here */
  void sched_db::install_sched_data(variable v, sched_data data) {
    sched_infov.push_back(data) ;
    install_sched_info(v,sched_info(sched_infov.size()-1)) ;
  }
  sched_db::sched_info &sched_db::get_sched_info(variable v) {
    vmap_type::iterator mi = vmap.find(remove_synonym(v)) ;
    if(mi == vmap.end()) 
      return get_sched_info(variable("EMPTY")) ; /*! assume variable("EMPTY") will be definitely in vmap, other infinite loop*/  
    return mi->second ;
  }

  void sched_db::set_variable_rotations(variableSet vars) {
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      get_sched_data(*vi).rotations += vars ;
    }
  }
  
  void sched_db::alias_variable(variable v, variable alias) {

    if(all_vars.inSet(v)) {
      if(all_vars.inSet(alias)) {
        if(get_exec_size() == 1) {
          cerr << "alias already in fact_db!" << endl ;
          cerr << "error found in alias_variable("<<v<<","<< alias<<")" << endl ;
        } else {
          debugout << "alias already in fact_db!" << endl ;
          debugout << "error found in alias_variable("<<v<<","<< alias<<")" << endl ;
        }
	detected_errors = true ;
        return ;
      }
      
      int ref = get_sched_info(v).sched_info_ref ;
      sched_infov[ref].aliases += alias ;
      install_sched_info(alias,sched_info(ref)) ;
      ref = get_sched_info(alias).sched_info_ref ;
      sched_infov[ref].antialiases += v ;
      
    } else if(all_vars.inSet(alias)) {
      alias_variable(alias,v) ;
    } else {
      if(get_exec_size() == 1)
        cerr << "neither variable " << v << ", nor " << alias << " exist in db, cannot create alias" << endl ;
      else
        debugout << "neither variable " << v << ", nor " << alias << " exist in db, cannot create alias" << endl ;
        
      detected_errors = true ;
    }
    
  }

  void sched_db::alias_variable(variable v, variable alias, fact_db &facts) {

    facts.synonym_variable(v,alias) ;

    if(all_vars.inSet(v)) {
      if(all_vars.inSet(alias)) {
        if(facts.get_comm_size() == 1) {
          cerr << "alias already in fact_db!" << endl ;
          cerr << "error found in alias_variable("<<v<<","<< alias<<")" << endl ;
        } else {
          debugout << "alias already in fact_db!" << endl ;
          debugout << "error found in alias_variable("<<v<<","<< alias<<")" << endl ;
        }
        detected_errors = true ;
        return ;
      }
      
      int ref = get_sched_info(v).sched_info_ref ;
      sched_infov[ref].aliases += alias ;
      install_sched_info(alias,sched_info(ref)) ;
      ref = get_sched_info(alias).sched_info_ref ;
      sched_infov[ref].antialiases += v ;
      
    } else if(all_vars.inSet(alias)) {
      alias_variable(alias,v, facts) ;
    } else {
      if(facts.get_comm_size() == 1)
        cerr << "neither variable " << v << ", nor " << alias << " exist in db, cannot create alias" << endl ;
      else
        debugout << "neither variable " << v << ", nor " << alias << " exist in db, cannot create alias" << endl ;
        
      detected_errors = true ;
    }
    
  }
  
  void sched_db::synonym_variable(variable v, variable synonym) {
    v = remove_synonym(v) ;
    vmap_type::iterator vmi ;
    if((vmi = vmap.find(synonym)) != vmap.end()) {
      if(get_exec_size() == 1) {
        cerr << "synonym already in fact_db!" << endl ;
        cerr << "error found in synonym_variable("<<v<<","<<synonym<<")"<<endl;
      } else {
        debugout << "synonym already in fact_db!" << endl ;
        debugout << "error found in synonym_variable("<<v<<","<<synonym<<")"<<endl;
      }
      //      detected_errors = true ;
      //      return ;
    }

    sched_info &vfinfo = get_sched_info(v) ;
    
    vfinfo.synonyms += synonym ;
    synonyms[synonym] = v ;
    all_vars += synonym ;
    sched_info &sfinfo = get_sched_info(synonym) ;
    sfinfo.synonyms += v ;
  }
  
  void sched_db::synonym_variable(variable v, variable synonym, fact_db &facts) {
    facts.synonym_variable(v, synonym) ;
    v = remove_synonym(v) ;
    vmap_type::iterator vmi ;
    if((vmi = vmap.find(synonym)) != vmap.end()) {
      if(facts.get_comm_size() == 1) {
        cerr << "synonym already in fact_db!" << endl ;
        cerr << "error found in synonym_variable("<<v<<","<<synonym<<")"<<endl;
      } else {
        debugout << "synonym already in fact_db!" << endl ;
        debugout << "error found in synonym_variable("<<v<<","<<synonym<<")"<<endl;
      }
      detected_errors = true ;
      return ;
    }

    sched_info &vfinfo = get_sched_info(v) ;
    
    vfinfo.synonyms += synonym ;
    synonyms[synonym] = v ;
    all_vars += synonym ;
    sched_info &sfinfo = get_sched_info(synonym) ;
    sfinfo.synonyms += v ;
  }
  
 
  
  void sched_db::set_existential_info(variable v, rule f, entitySet x) {
#ifdef VERBOSE
    debugout << "set_existential_info(v="<< v << ",rule="<< f << ",set=" << x << endl ;
#endif
    sched_info &finfo = get_sched_info(v) ;
#ifdef VERBOSE
    debugout << "finfo.existence = " << finfo.existence << endl ;
#endif
    if((finfo.existence & x) != EMPTY) {
      ruleSet conflicts ;
      map<rule,existential_info>::iterator mi ;
      for(mi=finfo.exist_map.begin();mi!=finfo.exist_map.end();++mi) {
        if(mi->first == f) {
          continue ;
        }
        if((mi->second.exists & x) != EMPTY) {
          const vector<string> &p1 = v.get_info().priority ;
          const vector<string> &p2 = mi->second.v.get_info().priority ;
          for(int i=0,j=p1.size()-1,k=p2.size()-1;
              i<min(int(p1.size()),int(p2.size()));
              ++i,--j,--k)
            
            if(p1[j] != p2[k]) {
              conflicts += mi->first ;
              if(get_exec_size() == 1)
                cerr << "adding to conflicts because " << p1[i]
                     << "!=" << p2[i]<<endl ;
              else
                debugout << "adding to conflicts because " << p1[i]
                         << "!=" << p2[i]<<endl ;
                
              detected_errors = true ;
            }
          if(p1.size() == p2.size()) {
            conflicts += mi->first ;
          } else if(p1.size() > p2.size()) {
            mi->second.exists -= x ;
            //cerr << f << " has priority over " << mi->first << endl ;
          } else {
            x -= mi->second.exists ;
            //cerr << mi->first << " has priority over " << f << endl ;
          }
        }
      }
      if(conflicts != EMPTY && v.get_info().name != string("OUTPUT")) {
        if(get_exec_size() == 1) {
          cerr << "rule " << f << " conflicts with " << conflicts << endl ;
          cerr << "conflicting entities are " << (finfo.existence & x) << endl ;
        } else {
          debugout << "rule " << f << " conflicts with " << conflicts << endl ;
          debugout << "conflicting entities are " << (finfo.existence & x) << endl ;
        }          
        detected_errors = true ;
        //        debugger_() ;
        //        exit(-1) ;
      }
    }
    existential_info &einfo = finfo.exist_map[f] ;

#ifdef VERBOSE
    debugout << "einfo.v=" << einfo.v << ",v=" << v << ",x=" << x << ",rule=" << f << endl  ;
#endif
    einfo.v = v ;
    einfo.exists += x ;
    
    finfo.existence += x ;
  }
  
  entitySet sched_db::variable_existence(variable v) {
    if(v.get_info().tvar)  // if time variable, variable exists everywhere
      return ~EMPTY ;
    return get_sched_info(v).existence ;
  }
  
  ruleSet sched_db::get_existential_rules(variable v) {
    std::map<rule,existential_info>::const_iterator mi ;
     sched_info &finfo = get_sched_info(v) ;
    ruleSet rules ;
    for(mi=finfo.exist_map.begin();mi!=finfo.exist_map.end();++mi)
      rules += mi->first ;
    return rules ;
  }
  
  void sched_db::variable_request(variable v, entitySet e) {
    if(v.get_info().tvar)
      return ;
    get_sched_info(v).requested += e ;
  }
  
  entitySet sched_db::get_variable_request(rule f, variable v) {
    sched_info &finfo = get_sched_info(v) ;
    map<rule,existential_info>::iterator mi = finfo.exist_map.find(f) ;
    if(mi == finfo.exist_map.end())
      return EMPTY ;
    return mi->second.exists & finfo.requested ;
  }

  void sched_db::set_variable_type(variable v, storeRepP st, fact_db &facts) {
    // creates an intensional fact since this is Loci deduced fact
    facts.create_intensional_fact(v, st) ;
    if(!all_vars.inSet(v)) {
      install_sched_data(v, sched_data(v, st)) ;
    }
  }

  entitySet sched_db::image(variable v, entitySet e) {
    sched_data &fdata = get_sched_data(v) ;
    if(!fdata.ismap)
      return EMPTY ;
    map<entitySet,entitySet>::const_iterator ii = fdata.imageMap.find(e) ;
    if(ii != fdata.imageMap.end())
      return ii->second ;
    else
      return fdata.imageMap[e] = get_sched_data(v).minfo->image(e) ;
  }
   
  pair<entitySet,entitySet> sched_db::preimage(variable v, entitySet e) {
    sched_data &fdata = get_sched_data(v) ;
    if(!fdata.ismap)
      return make_pair(EMPTY,EMPTY) ;
    map<entitySet,pair<entitySet,entitySet> >::const_iterator ii ;
    ii = fdata.preimageMap.find(e) ;
    if(ii != fdata.preimageMap.end())
      return ii->second ;
    else
      return fdata.preimageMap[e] = get_sched_data(v).minfo->preimage(e) ;
  }
  
  void sched_db::add_policy(variable v, duplicate_policy p) {
    unsigned int policy = get_policy(v);
    unsigned int temp = 1;
    switch(p) {
    case NEVER:
      break;
    case ALWAYS:
      temp = temp << 1;
      break;
    case MODEL_BASED:
      temp = temp << 2;
      break;
    }
    policy |= temp;
    variableSet syns = get_synonyms(v);
    for(variableSet::const_iterator vi = syns.begin(); vi != syns.end(); vi++)
      set_policy(*vi, policy);
  }

  bool sched_db::is_policy(variable v, duplicate_policy p) {
    unsigned int policy = get_policy(v);
    unsigned int temp = 1;
    switch(p) {
    case NEVER:
      break;
    case ALWAYS:
      temp = temp << 1;
      break;
    case MODEL_BASED:
      temp = temp << 2;
      break;
    }
    return (temp & policy);
  }

  void sched_db::add_model_info(double comm_ts, double comm_tw,
				const map<rule, pair<double, double> > &comp_info) {
    comm_model.ts = comm_ts;
    comm_model.tw = comm_tw;
    
    for(map<rule, pair<double, double> >::const_iterator mi = comp_info.begin();
	mi != comp_info.end(); mi++) {
      model tmpModel(mi->second.first, mi->second.second);
      comp_model[mi->first] = tmpModel;
    }
  }

  std::ostream &sched_db::print_summary(fact_db &facts, std::ostream &s) {
    s << "Summary of Existential deduction:" << endl ;
    std::map<variable,sched_info>::const_iterator mi ;
    for(mi=vmap.begin();mi!=vmap.end();++mi) {
      storeRepP sp = facts.get_variable(mi->first) ;
      if(isMAP(sp))
	s << " Container = MAP " << endl ;
      else if(isSTORE(sp))
	s << " Container = STORE " << endl ;
      else if(isPARAMETER(sp))
	s << " Container = PARAMETER " << endl ;
      else if(isBLACKBOX(sp))
	s << " Container = BLACKBOX " << endl;
      else if(isCONSTRAINT(sp))
	s << " Container = CONSTRAINT " << endl ;
      
      //      double size = 0 ;
      //      size = sp->pack_size(mi->second.requested) / 1000000 ;
      s << mi->first << " " <<mi->second.synonyms << " "<< mi->second.existence
        << " request= " << mi->second.requested <<endl ; //<< "  size requested = " << size << " MB " << endl ;
    }
    return s ;
  }

  std::vector<std::pair<variable,entitySet> > sched_db::get_send_entities(variableSet eset, send_entities_type e) {

    std::vector<std::pair<variable,entitySet> > re;
    variableSet::const_iterator vi;
    std::map<variable, entitySet>::const_iterator mi;
    
    variableSet send_vars;
    for(vi = eset.begin(); vi !=eset.end(); vi++){
      variable v = variable(*vi);
      ruleSet rs = get_existential_rules(v) ;
      for(ruleSet::const_iterator rsi = rs.begin(); rsi != rs.end(); ++rsi) {
        if(rule_has_mapping_in_output(*rsi)) {
          send_vars += v;
          break;
        }
      }
    }

    if(e == BARRIER){
      for(vi = send_vars.begin(); vi != send_vars.end(); vi++){
        variable v = variable(*vi);
        mi = barrier_send_entities_map.find(v);
        if(mi != barrier_send_entities_map.end()){
          re.push_back(make_pair(v, mi->second));
        }else{
          debugout << "WARNING: for variable " << v << " barrier_send_entities_map is read before it's written" << endl;
        }
      }
    }else if( e == RECURSE_PRE){
      for(vi = send_vars.begin(); vi !=send_vars.end(); vi++){
        variable v = variable(*vi);
        mi = recurse_pre_send_entities_map.find(v);
        if(mi != recurse_pre_send_entities_map.end()){
          re.push_back(make_pair(v, mi->second));
        }else{
          debugout << "WARNING: for variable " << v << " recurse_pre_send_entities_map is read before it's written" << endl;
        }
      }
    }else{
      debugout<<"ERROR: unrecognized send_entities_type in get_send_entities() " << endl;
    }
    return re;
  }
  

    
  void sched_db::update_send_entities( const std::vector<std::pair<variable,entitySet> >& evec, send_entities_type e){
    
    if(e == BARRIER){
      for(unsigned int vi = 0; vi < evec.size(); vi++){
        variable v = evec[vi].first;
        map<variable, entitySet>::const_iterator mi = barrier_send_entities_map.find(v);
        if(mi != barrier_send_entities_map.end()){
          debugout << "WARNING: for variable " << v << " barrier_send_entities_map is written more then once" << endl;
        }
        barrier_send_entities_map[v] = (evec[vi]).second;
      }
    }else if(e == RECURSE_PRE){
      for(unsigned int vi = 0; vi < evec.size(); vi++){
        variable v = evec[vi].first;
        map<variable, entitySet>::const_iterator mi = recurse_pre_send_entities_map.find(v);
        if(mi != recurse_pre_send_entities_map.end()){
          debugout << "WARNING: for variable " << v << " barrier_send_entities_map is written more then once" << endl;
        }
        recurse_pre_send_entities_map[v] = (evec[vi]).second;
      }
    }else{
      debugout<<"ERROR: unrecognized send_entities_type in update_send_entities" << endl;
    }
  }
  
  /*!sort_comm_map() composes a list from the map so that in the list all send info is in front of  recv info
    and the order of comm info is
    for(pi = procs.begin(); pi != procs.end(); pi++){
    for(vi = vset.begin(); vi!= vset.end(); vi++){
    comm_info(*pi, *vi);
    }}
  */
  std::list<comm_info> sort_comm_map(variableSet vset,
                                     const  std::map<variable, std::list<comm_info> >& clist_m,
                                     fact_db& facts,
                                     sched_db::list_type e){
    
    
    vector<std::pair<int,std::vector<send_var_info> > > send_info ;
    vector<std::pair<int,std::vector<recv_var_info> > > recv_info ;


    // First collect information from slist

    HASH_MAP(int,vector<send_var_info>) send_data ;
    HASH_MAP(int,vector<recv_var_info>) recv_data ;
    intervalSet send_procs, recv_procs ;
    for(variableSet::const_iterator vi = vset.begin(); vi != vset.end(); vi++){
      variable va = variable(*vi);
      list<comm_info> slist;
      std::map<variable, std::list<comm_info> >::const_iterator mi = clist_m.find(va);
      if(mi != clist_m.end()){
        //        debugout <<"get variable " << va <<" from list " << e << endl;
        slist = mi->second;
        list<comm_info>::const_iterator cli ;
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
    const int nrecv = recv_info.size() ;
    const int nsend = send_info.size() ;
   
    // Pack the buffer for sending
    for(int i=0;i<nsend;++i) {
      for(size_t j=0;j<send_info[i].second.size();++j) {
        comm_info ci ;
        ci.v = send_info[i].second[j].v ;
        ci.processor = send_info[i].first ;
        ci.send_set = send_info[i].second[j].set ;
        clist.push_back(ci) ;
      }
    }
    for(int i=0;i<nrecv;++i) {
      for(size_t j=0;j<recv_info[i].second.size();++j) {
       	comm_info ci ;
        ci.v = recv_info[i].second[j].v ;
        ci.processor = recv_info[i].first ;
        ci.recv_set = recv_info[i].second[j].seq ;
	clist.push_back(ci) ;
      }
    }
    return clist ;
  }

  std::list<comm_info> sched_db::get_comm_info_list(variableSet eset, fact_db& facts, list_type e) const{
    list<comm_info> re;
    if(e==BARRIER_CLIST){
      re =  sort_comm_map(eset, barrier_clist_map, facts, e);
    }else if(e==BARRIER_PLIST){
      re = sort_comm_map(eset, barrier_plist_map, facts, e);
    }else if(e==REDUCE_RLIST){
      re = sort_comm_map(eset, reduce_rlist_map, facts, e);
    }else if(e==REDUCE_CLIST){
      re = sort_comm_map(eset, reduce_clist_map, facts, e);
    }else if(e==LOOP_ADVANCE_LIST){
      re = sort_comm_map(eset, loop_advance_list_map, facts, e);
    }else if(e==RECURSE_CLIST){
      re = sort_comm_map(eset, recurse_clist_map, facts, e);
    }else if(e==RECURSE_PRE_CLIST){
      re = sort_comm_map(eset, recurse_pre_clist_map, facts, e);
    }else if(e==RECURSE_POST_CLIST){
      re = sort_comm_map(eset, recurse_post_clist_map, facts, e);
    }else if(e==RECURSE_PRE_PLIST){
      re = sort_comm_map(eset, recurse_pre_plist_map, facts, e);
    }else{
      debugout << "ERROR: unrecognized list_type " << e << " in get_comm_info_list" << endl;
    }
    return re;
  }
  void sched_db::update_comm_info_list(const std::list<comm_info>& elist, list_type e){
    list<comm_info>::const_iterator cli ;
    // map<variable, list<comm_info> >::iterator mi; 
    for(cli=elist.begin();cli!=elist.end();++cli) {
      variable v = cli->v ;
      if(e==BARRIER_CLIST){
        barrier_clist_map[v].push_back( *cli); 
      }else if(e==BARRIER_PLIST){
        barrier_plist_map[v].push_back( *cli);
      }else if(e==REDUCE_RLIST){
        reduce_rlist_map[v].push_back( *cli); 
      }else if(e==REDUCE_CLIST){
        reduce_clist_map[v].push_back( *cli); 
      }else if(e==LOOP_ADVANCE_LIST){
        loop_advance_list_map[v].push_back( *cli); 
      }else if(e==RECURSE_CLIST){
        recurse_clist_map[v].push_back( *cli); 
      }else if(e==RECURSE_PRE_CLIST){
        recurse_pre_clist_map[v].push_back( *cli); 
      }else if(e==RECURSE_POST_CLIST){
        recurse_post_clist_map[v].push_back( *cli); 
      }else if(e==RECURSE_PRE_PLIST){
        recurse_pre_plist_map[v].push_back( *cli);  
      } else{
        debugout << "ERROR: unrecognized list_type " << e << " in update_comm_info_list" << endl;
      }
        
    }     
  }
}
  
