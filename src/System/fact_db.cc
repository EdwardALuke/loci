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
#include <typeinfo>
#include <Tools/except.h>

#include <fact_db.h>
#include <rule.h>
#include <constraint.h>
#include <Tools/debugger.h>
#include <DStore.h>
#include "dist_tools.h"

#include <Loci_Datatypes.h>

using std::string ; 
using std::map ;
using std::make_pair ;
using std::vector ;
using std::list ;
using std::sort ;
using std::pair ;

using std::istream ;
using std::ostream ;
using std::ostringstream ;
using std::endl ;
using std::ios ;

#include <Tools/parse.h>

namespace Loci {
  extern int MPI_processes ;
  extern int MPI_rank ;
  extern fact_db *exec_current_fact_db;
  int factdb_allocated_base = 0 ;

  fact_db::fact_db() {
    distributed_info = 0 ;
    maximum_allocated = factdb_allocated_base+0 ;
    vector<entitySet> baseline_ptn ;
    for(int i = 0; i < MPI_processes; ++i) {
      baseline_ptn.push_back(EMPTY) ;
    }
    init_ptn.push_back(baseline_ptn) ;
    //    exec_current_fact_db = this;
    dist_from_start = 0 ;
    constraint EMPTY_constraint ;
    EMPTY_constraint = EMPTY ;
    create_fact("EMPTY",EMPTY_constraint) ;
    constraint UNIVERSE_constraint ;
    UNIVERSE_constraint = ~EMPTY ;
    create_fact("UNIVERSE",UNIVERSE_constraint) ;
  }

  fact_db::~fact_db() {}

  void
  fact_db::copy_all_from(const fact_db& f) {
    init_ptn = f.init_ptn ;
    gmax_alloc = f.gmax_alloc ;
    keyDomainName = f.keyDomainName ;
    global_comp_entities = f.global_comp_entities ;
    synonyms = f.synonyms ;
    maximum_allocated = f.maximum_allocated ;
    dist_from_start = f.dist_from_start ;
    fmap = f.fmap ;
    tmap = f.tmap ;
    nspace_vec = f.nspace_vec ;
    extensional_facts = f.extensional_facts ;

#ifdef DYNAMICSCHDEDULING
    // create new keyspaces
    keyspace.clear() ;
    const std::map<std::string,KeySpaceP>& fks = f.keyspace ;
    for(std::map<std::string,KeySpaceP>::const_iterator
          mi=fks.begin();mi!=fks.end();++mi) {
      keyspace[mi->first] = mi->second->new_keyspace() ;
    }      
    // here is something very important, we'll need to
    // assign the current fact_db's synonym record to
    // all the keyspaces here
    for(std::map<std::string,KeySpaceP>::iterator
          mi=keyspace.begin();mi!=keyspace.end();++mi) {
      (mi->second)->set_synonyms(&synonyms) ;
    }
    // deep copy the key manager
    if(f.key_manager != 0) {
      key_manager = f.key_manager->clone() ;
    }
#endif
    /* we cannot use the following direct assignment
       to copy the distributed_info from f since
       distributed_info is a NPTR pointer and is
       reference counted this would not be a true copy
    */
    // distributed_info = f.distributed_info ;
    distribute_infoP& df = distributed_info ;
    const distribute_infoP& fdf = f.distributed_info ;
    if(fdf == 0) {
      df = 0 ;
      return ;
    }
    df = new distribute_info ;
    df->myid = fdf->myid ;
    df->isDistributed = fdf->isDistributed ;
    // we make a deep copy of the maps
    entitySet l2g_alloc = fdf->l2g.domain() ;
    df->l2g.allocate(l2g_alloc) ;
    df->l2f.allocate(l2g_alloc) ;
    df->key_domain.allocate(l2g_alloc) ;
    for(entitySet::const_iterator ei=l2g_alloc.begin();
        ei!=l2g_alloc.end();++ei) {
      df->l2g[*ei] = fdf->l2g[*ei] ;
      df->l2f[*ei] = fdf->l2f[*ei] ;
      df->key_domain[*ei] = fdf->key_domain[*ei] ;
    }

    // num_keydomains
    int num_keydomains = fdf->g2lv.size() ; ;
    df->g2lv = vector<dMap>(num_keydomains) ;
    for(int i=0;i<num_keydomains;++i) {
      dMap tmp ;
      df->g2lv[i].setRep(tmp.Rep()) ;
    }

    for(int i=0;i<num_keydomains;++i) {
      entitySet g2l_alloc = fdf->g2lv[i].domain() ;
      df->g2lv[i].allocate(g2l_alloc) ;
      for(entitySet::const_iterator ei=g2l_alloc.begin();
	  ei!=g2l_alloc.end();++ei)
	df->g2lv[i][*ei] = fdf->g2lv[i][*ei] ;
    }
    df->my_entities = fdf->my_entities ;
    df->comp_entities = fdf->comp_entities ;
    df->copy = fdf->copy ;
    df->xmit = fdf->xmit ;
    df->copy_total_size = fdf->copy_total_size ;
    df->xmit_total_size = fdf->xmit_total_size ;
    //      df->remap = fdf->remap ;
#ifdef LOCI_COMPAT_MODE1
    df->g2fv = fdf->g2fv ;
    df->g2l.setRep(df->g2lv[0].Rep()) ;
    df->g2f.setRep(df->g2fv[0].Rep()) ;
#endif
  }
  
  void fact_db::set_maximum_allocated(int i) {
    maximum_allocated = i ;
  }
  
  int fact_db::getKeyDomain(std::string name) {
    for(size_t i=0;i<keyDomainName.size();++i)
      if(name == keyDomainName[i])
	return i ;
    gmax_alloc.push_back(maximum_allocated) ;
    keyDomainName.push_back(name) ;
    if(dist_from_start) {
      while(distributed_info->g2fv.size() < keyDomainName.size())
	distributed_info->g2fv.push_back(dMap()) ;
      while(distributed_info->g2lv.size() < keyDomainName.size())
	distributed_info->g2lv.push_back(dMap()) ;
    }
    vector<entitySet> baseline_ptn ;
    for(int i = 0; i < MPI_processes; ++i) {
      baseline_ptn.push_back(EMPTY) ;
    }
    while(init_ptn.size() < keyDomainName.size())
      init_ptn.push_back(baseline_ptn) ;
    return keyDomainName.size()-1 ;
  }

  void fact_db::synonym_variable(variable v, variable synonym) {
    // Find all variables that should be synonymous with v
    variableSet synonym_set ;
    std::map<variable,variable>::const_iterator mi ;
    while((mi=synonyms.find(v)) != synonyms.end()) {
      synonym_set += v ;
      v = mi->second ;
    }
    variable s = synonym;
    while((mi=synonyms.find(s)) != synonyms.end()) {
      synonym_set += s ;
      s = mi->second ;
    }
    synonym_set += s ;
 
    // If the two are already synonymous, we are done
    if(s == v)
      return ;

    // Sanity check, make sure v exists
    std::map<variable,fact_info>::iterator vmi, vmj ;
    if((vmi = fmap.find(v)) == fmap.end()) {
      cerr << "WARNING: synonym_variable("<<v<<","<<synonym<<")"<<endl ;
      cerr << "WARNING: type not known for target of synonym, ignoring operation" << endl ;
      //      abort() ;
      return ;
    }
    // If the synonym already points to a different variable instance,
    // remove it
    if((vmj = fmap.find(s)) != fmap.end()) {
      fact_info &finfo = vmj->second ;
      if((finfo.data_rep->domain() != EMPTY &&
          !isPARAMETER(finfo.data_rep->getRep()))) {
        cerr << "unable to define synonym variable " << synonym
             << " when variable already created in db. "  << endl ;
        cerr << "variable v = " << v << endl ;
        abort() ;
      }
      remove_variable(synonym) ;
    }
    
    // Add new synonyms so that they point to v
    for(variableSet::const_iterator vi = synonym_set.begin();
        vi!=synonym_set.end();
        ++vi) {
      synonyms[*vi] = v ;
    }
  }
  
  
  void fact_db::update_fact(variable v, storeRepP st) {
    //if st is STORE or MAP, update maximum_allocated
    if(isMAP(st) || isSTORE(st)) {
      int max_val = st->domain().Max() ;
      if(gmax_alloc.size()==0)
	getKeyDomain("Main") ;
      int kd = st->getDomainKeySpace() ;
      gmax_alloc[kd] = max(gmax_alloc[kd],max_val+1) ;
    }
    //add namespace
    variable tmp_v = add_namespace(v) ;

    //tmp_v should not have synonyms
    warn(synonyms.find(tmp_v) != synonyms.end()) ;
    std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
    //tmp_v should be in fmap
    if(mi != fmap.end()) {
      mi->second.data_rep->setRep(st->getRep()) ;
    } else
      cerr << "warning: update_fact: fact does not exist for variable " << tmp_v
	   << endl ;
  }

  variable fact_db::add_namespace(variable v) const {
    variable tmp_v ;
    if((nspace_vec.size() > 0) && (v.get_info().namespac.size() == 0)) {
      tmp_v = v ;
      for(size_t i = 0; i < nspace_vec.size(); ++i)
        tmp_v = tmp_v.add_namespace(nspace_vec[i]) ;
    } else
      tmp_v = v ;
    return tmp_v ;
  }
  
  void fact_db::create_pure_fact(const variable& v, storeRepP st) {
    if(isMAP(st) || isSTORE(st)) {
      int max_val = st->domain().Max() ;
      if(gmax_alloc.size()==0)
	getKeyDomain("Main") ;
      int kd = st->getDomainKeySpace() ;
      gmax_alloc[kd] = max(gmax_alloc[kd],max_val+1) ;
    }
    variable tmp_v ;
    tmp_v = v ;
    
    if(synonyms.find(tmp_v) != synonyms.end()) {
      tmp_v = remove_synonym(tmp_v) ;
      std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
      if(mi==fmap.end()) {
        fmap[tmp_v].data_rep = new store_ref ;
        fmap[tmp_v].data_rep->setRep(st->getRep()) ;
      } else {
        if(typeid(st->getRep()) != typeid(mi->second.data_rep->getRep())) {
          cerr << "set_variable_type() method of fact_db changing type for variable " << tmp_v << endl ;
        }
        mi->second.data_rep->setRep(st->getRep()) ;
      }
      return ;
    }
    
    std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
    if(mi != fmap.end()) {
      cerr << "WARNING: fact_db::set_variable_type retyping variable "
	   << tmp_v << endl ;
      mi->second.data_rep->setRep(st->getRep()) ;
    } else {
      fmap[tmp_v].data_rep = new store_ref ;
      fmap[tmp_v].data_rep->setRep(st->getRep()) ;
    }
  } 

  /*! remove from synonym and fmap */
  void fact_db::remove_variable(variable v) {
    std::map<variable, variable>::iterator si ;
    std::map<variable, fact_info>::iterator mi ;
    if((si=synonyms.find(v)) != synonyms.end()) {
      variable real_var = remove_synonym(v) ;
      synonyms.erase(si) ;
      remove_variable(real_var) ;
    } else if((mi=fmap.find(v)) != fmap.end()) {
      // First remove any synonyms to this variable.
      variableSet syn_vars ;
      vector<map<variable,variable>::iterator > lrm ;
      for(si=synonyms.begin();si!=synonyms.end();++si)
        if(si->second == v)
          syn_vars += si->first ;
      for(variableSet::const_iterator vi=syn_vars.begin();
          vi!=syn_vars.end();++vi) {
        if((si=synonyms.find(*vi)) != synonyms.end()) {
          synonyms.erase(si) ;
        }
      }

      // Now erase the variable
      fmap.erase(mi) ;
    }
  }
  
  variableSet fact_db::get_typed_variables() const {
    std::map<variable, fact_info>::const_iterator mi ;
    std::map<variable, variable>::const_iterator si ;
    variableSet all_vars ;
    for(mi=fmap.begin();mi!=fmap.end();++mi)
      all_vars += mi->first ;
    
    for(si=synonyms.begin();si!=synonyms.end();++si)
      all_vars += si->first ;
    return all_vars ;
  }

  void fact_db::make_extensional_fact(const variable& v) {
    // if it is already extensional_fact, then we do nothing
    if(extensional_facts.inSet(v))
      return ;
    std::map<variable, fact_info>::iterator mi = fmap.find(v) ;
    // the fact must exist in the fact_db, otherwise we do nothing
    if(mi != fmap.end())
      extensional_facts += v ;
  }
  void fact_db::make_intensional_fact(const variable& v) {
    // We perform action only if it is already an extensional_fact
    // which implies that the fact exists in the fact_db
    if(extensional_facts.inSet(v))
      extensional_facts -= v ;
  }

  int fact_db::get_max_alloc(int kd) {
    if(gmax_alloc.size() == 0) {
      getKeyDomain("Main") ;
    }
    return gmax_alloc[kd] ;
  }
  std::pair<entitySet, entitySet> fact_db::get_dist_alloc(int size, size_t kd) {

    if(gmax_alloc.size() == 0) {
      getKeyDomain("Main") ;
    }
    if(MPI_processes > 1) {
      if(!dist_from_start) {
	dist_from_start = 1 ;
	distributed_info = new distribute_info;
	int num_keyspace = gmax_alloc.size() ;
	for(int i=0;i<num_keyspace;++i) {
	  distributed_info->g2fv.push_back(dMap()) ;
	  distributed_info->g2lv.push_back(dMap()) ;
	}
#ifdef LOCI_COMPAT_MODE1
	distributed_info->g2f.setRep(distributed_info->g2fv[0].Rep()) ;
	distributed_info->g2l.setRep(distributed_info->g2lv[0].Rep()) ;
#endif
      }


      int* send_buf = new int[MPI_processes] ;
      int* size_send = new int[MPI_processes] ;
      int* size_recv = new int[MPI_processes] ;
      int* recv_buf = new int[MPI_processes] ;
      for(int i = 0; i < MPI_processes; ++i) {
	send_buf[i] = gmax_alloc[kd] ;
	size_send[i] = size ;
      } 
      MPI_Alltoall(send_buf, 1, MPI_INT, recv_buf, 1, MPI_INT, MPI_COMM_WORLD) ;
      MPI_Alltoall(size_send, 1, MPI_INT, size_recv, 1, MPI_INT, MPI_COMM_WORLD) ;
      std::sort(recv_buf, recv_buf+MPI_processes) ;
      gmax_alloc[kd] = recv_buf[MPI_processes-1] ;
      int local_max = gmax_alloc[kd] ;
      int global_max = 0 ;
      for(int i = 0; i < MPI_rank; ++i)
	local_max += size_recv[i] ;
      for(int i = 0; i < MPI_processes; ++i) 
	global_max += size_recv[i] ;
      
      for(int i = 0 ; i < MPI_processes; ++i) {
	int local = gmax_alloc[kd] ;
	for(int j = 0; j < i; ++j)
	  local += size_recv[j] ;
	if(size_recv[i] > 0 )
	  init_ptn[kd][i] += interval(local, local+size_recv[i]-1) ;
      }
      entitySet local_ivl, global_ivl ;
      if(size > 0 )
	local_ivl = entitySet(interval(local_max, local_max + size - 1)) ;
      else {
	local_ivl = EMPTY ;
      }
      
      global_ivl = entitySet(interval(gmax_alloc[kd], gmax_alloc[kd]+global_max-1)) ;

      global_comp_entities += local_ivl;
      
      delete [] send_buf ;
      delete [] recv_buf ;
      delete [] size_send ;
      delete [] size_recv ;
      gmax_alloc[kd] = max(gmax_alloc[kd],global_ivl.Max()+1) ;
      return(make_pair(local_ivl, global_ivl)) ;
    }
    entitySet alloc = entitySet(interval(gmax_alloc[kd],gmax_alloc[kd]+size-1)) ;
    gmax_alloc[kd] += size ;
    
    init_ptn[kd][0] += alloc ;
    global_comp_entities += alloc;
    return (make_pair(alloc, alloc)) ;
  }
    
  void fact_db::update_remap(const std::vector<std::pair<int, int> > &remap_update, size_t kd) {
    if(Loci::MPI_processes > 1) {
      warn(!dist_from_start);
      fatal(distributed_info == NULL);
      
      for(std::vector<std::pair<int, int> >::const_iterator vi = remap_update.begin(); vi != remap_update.end(); vi++) {
        //	distributed_info->remap[vi->first] = vi->second;
        distributed_info->g2fv[kd][vi->second] = vi->first ;
      }
#ifdef LOCI_COMPAT_MODE1
      distributed_info->g2f.setRep(distributed_info->g2fv[0].Rep()) ;
      //      distributed_info->g2l.setRep(distributed_info->g2lv[0].Rep()) ;
#endif
    }
  }

  std::pair<entitySet, entitySet> fact_db::get_distributed_alloc(int size, size_t kd) {
    pair<entitySet, entitySet> allocation = get_dist_alloc(size, kd);
    vector<pair<int, int> > remap_update;
   
    FORALL(allocation.first, ai) {
      remap_update.push_back(make_pair(ai, ai));
    }ENDFORALL;

    update_remap(remap_update,kd);    
    return allocation;
  }

  std::pair<entitySet, entitySet> fact_db::get_distributed_alloc(const std::vector<int> &remap_entities, size_t kd) {
    pair<entitySet, entitySet> allocation = get_dist_alloc(remap_entities.size(),kd);
    vector<pair<int, int> > remap_update;
    
    int i = 0;   
    FORALL(allocation.first, ai) {
      remap_update.push_back(make_pair(remap_entities[i], ai));
      i++;
    }ENDFORALL;

    update_remap(remap_update,kd);
    return allocation;
  }

  std::pair<entitySet, entitySet> fact_db::get_distributed_alloc(int size, int offset,size_t kd) {  
    pair<entitySet, entitySet> allocation = get_dist_alloc(size,kd);
    vector<pair<int, int> > remap_update;
    
    FORALL(allocation.first, ai) {
      remap_update.push_back(make_pair(ai+offset, ai));
    }ENDFORALL;
    
    update_remap(remap_update,kd);
    return allocation;
  }
  
  storeRepP fact_db::get_variable(variable v) {
    variable tmp_v = add_namespace(v) ;
    tmp_v = remove_synonym(tmp_v) ;
    std::map<variable, fact_info>::iterator mi =
      fmap.find(remove_synonym(tmp_v)) ;
    if(mi == fmap.end()) {
      return storeRepP(0) ;
    }
    else
      return storeRepP(mi->second.data_rep) ;
  }
  
  fact_db::distribute_infoP fact_db::get_distribute_info() {
    return(distributed_info);
  }
  
  void fact_db::put_distribute_info(distribute_infoP dp) {
    distributed_info = dp ;
  }
 
  bool fact_db::isDistributed() {
    if(distributed_info == 0)
      return 0 ;
    else 
      return 1 ;
  }
  
  void fact_db::rotate_vars(const std::list<variable> &lvars) {
    list<variable>::const_iterator jj ;
    jj = lvars.begin() ;
    storeRepP cp = fmap[remove_synonym(*jj)].data_rep->getRep() ;
    ++jj ;
    if(jj != lvars.end()) {
      for(;jj!=lvars.end();++jj) {
        fact_info &fd = fmap[remove_synonym(*jj)] ;
        storeRepP tmp = fd.data_rep->getRep() ;
        fd.data_rep->setRep(cp) ;
        cp = tmp ;
      }
    }
    fmap[remove_synonym(lvars.front())].data_rep->setRep(cp) ;
  }

#ifdef DYNAMICSCHEDULING
  void fact_db::adjust_rotation_vars(const std::list<variable>& lvars) {
    list<variable>::const_iterator jj ;
    jj = lvars.begin() ;
    storeRepP cur = fmap[remove_synonym(*jj)].data_rep->getRep() ;
    entitySet cur_dom = cur->domain() ;

    ++jj ;
    if(jj != lvars.end()) {
      for(;jj!=lvars.end();++jj) {
        fact_info& fd = fmap[remove_synonym(*jj)] ;
        storeRepP his = fd.data_rep->getRep() ;
        // first we will erase the domains outside of cur_dom
        entitySet his_dom = his->domain() ;
        entitySet out = his_dom - cur_dom ;
        if(out != EMPTY) {
          his->erase(out) ;
        }
        // then copy those missing
        entitySet missing = cur_dom - his_dom ;
        if(missing != EMPTY) {
          his->copy(cur, missing) ;
        }
      }
    }
  }
#endif
  ostream &fact_db::write(ostream &s) const {
    std::map<variable, fact_info>::const_iterator vmi ;
    for(vmi=fmap.begin();vmi!=fmap.end();++vmi) {
      variable v=vmi->first;
      storeRepP storeRep = storeRepP(vmi->second.data_rep) ;
      entitySet en=storeRep->domain();
      std::string groupname = (v.get_info()).name;
      s << groupname << ":" ;
      storeRep->Print(s);
    }
    return s ;
  }
  
  istream &fact_db::read(istream &s) {
    bool syntax_error = false ;
    try {
      string vname ;
      parse::kill_white_space(s) ;

      if(s.peek()!='{') {
        throw StringError("format error in fact_db::read, missing '{'") ;
        return s ;
      }
      s.get() ;
      
      variableSet read_vars ;
      for(;;) {
        parse::kill_white_space(s) ;
        if(s.peek() == '}') {
          s.get() ;
          break ;
        }
        if(s.peek() == std::char_traits<char>::eof()) {
          throw StringError("unexpected EOF in fact_db::read") ;
        }
        parse::kill_white_space(s) ;
        if(parse::is_name(s)) 
          vname = parse::get_name(s) ;
        else {
          throw StringError("syntax error in fact_db::read") ;
        }
        parse::kill_white_space(s) ;
	if(s.peek() == '@') {
	  vname += '@' ;
	  s.get() ;
	  parse::kill_white_space(s) ;
	  if(parse::is_name(s)) 
	    vname += parse::get_name(s) ;
	  else {
	    ostringstream oss ;
	    oss << "syntax error in fact_db::read, no variable name after namspace '@' separator for namespace '" << vname << "'" ;
	    throw StringError(oss.str()) ;
	  }
	}
        variable var(vname) ;
        if(!parse::get_token(s,":")) {
          ostringstream oss ;
          oss << "syntax error in fact_db::read, no ':' separator for variable '" << var << "'" ;
          throw StringError(oss.str()) ;
        }

        if(read_vars.inSet(var)) {
          cerr << "WARNING: Variable '" << var << "' is redefined in fact_db input!!!" << endl ;
        }
        read_vars += var ;
        try {
          storeRepP vp = get_variable(var) ;
          if(vp == 0) {
            vp = get_variable_type(var) ;
            if(vp != 0) {
              create_fact(var,vp) ;
            }
            vp = get_variable(var) ;
          }
          if(vp == 0) {
            ostringstream oss ;
            oss << "variable named '" << vname
                << "' not found in database in fact_db::read." << endl ;
            throw StringError(oss.str()) ;
          }
          vp->Input(s) ;
          parse::kill_white_space(s) ;
          if(s.peek() != '}' && !parse::is_name(s)) {
            ostringstream oss ;
            oss << "syntax error in fact_db::read while reading variable '" << vname << "'" ;
            throw StringError(oss.str()) ;
          }
        } catch(const BasicException &err) {
          err.Print(cerr) ;
          cerr << "input failed for variable " << vname << endl ;
          cerr << "near the following text: \"" ;
          while(s.peek()!=EOF&&s.peek()!='\n'&&s.peek()!='\r')
            cerr << char(s.get()) ;
          cerr << "\"" << endl ;
         
 
          syntax_error = true ;
        }
      }
    } catch(const BasicException &err) {
      err.Print(cerr) ;
      throw StringError("read fact_db failed") ;
    }
    if(syntax_error) {
      cerr << "syntax error reading fact db" << endl ;
      throw StringError("read fact_db failed") ;
    }
    return s ;
  }

  void fact_db::setupDefaults(const rule_db &rdb) {
    ruleSet special_rules = rdb.get_default_rules() ;
    // first we process the default rules
    variableSet working_vars ;
    std::map<variable,std::pair<variable,rule> > var_def_rules ;
    std::map<variable,std::pair<variable,rule> >::const_iterator mi ;

    for(ruleSet::const_iterator ri=special_rules.begin();
        ri!=special_rules.end();++ri) {
      // first we need to create the facts in the fact_db
      variableSet targets = ri->targets() ;
      //      rule_implP rp = ri->get_rule_implP() ;
      if(targets.size() > 1) // ignore all default rules with more than 1 arg
        continue ;

      variable v = *targets.begin() ;
      // only process rules in current namspace
      if(nspace_vec.size()>0 &&v.get_info().namespac != nspace_vec) 
	continue ;
      variable dpv = v ;
      while(dpv.get_info().priority.size() != 0)
	dpv = dpv.drop_priority() ;
      mi = var_def_rules.find(dpv) ;
      std::pair<variable,rule> defrule(v,*ri) ;
      if(mi == var_def_rules.end()) {
	var_def_rules[dpv] = defrule ;
      } else {
	// check priority
	const variable vcheck = mi->second.first ;
	bool replace = false ;
	variable vsearch = v ;
	while(vsearch.get_info().priority.size() != 0) {
	  vsearch = vsearch.drop_priority() ;
	  if(vcheck == vsearch)
	    replace = true ;
	}
	if(replace) {
	  var_def_rules[dpv] = defrule ;
	}
      }

    }

    for(mi=var_def_rules.begin();mi!=var_def_rules.end();++mi) {
      // we need to get the storeRep for this variable
      rule_implP rp = mi->second.second.get_rule_implP() ;
      variable v = mi->first ;
      storeRepP srp = rp->get_store(mi->second.first) ;
      if(srp == 0) {
        ostringstream oss ;
        oss << "default rule " << mi->second.second 
	    << " unable to provide type for " << v
            << endl ;
	cerr << oss.str() << endl ;
        throw StringError(oss.str()) ;
      }
      std::map<variable, fact_info>::iterator fmi = fmap.find(v) ;
      if(fmi == fmap.end()) { // does not exist, so install
	//	srp->allocate(~EMPTY) ;
	create_fact(v,srp) ;
	rp->set_store(mi->second.first,srp) ;
        rp->compute(sequence(EMPTY)) ;
      }
    }
  }

  // this read in routine reads the .vars file and set
  // the facts according to the default and optional rules
  std::istream& fact_db::read_vars(std::istream& s, const rule_db& rdb) {
    bool syntax_error = false ;
    try {
      // first of all, we need to process the default and optional rules

      setupDefaults(rdb) ;
      // now process the optional rules
      ruleSet special_rules = rdb.get_optional_rules() ;
      for(ruleSet::const_iterator ri=special_rules.begin();
          ri!=special_rules.end();++ri) {
        // first we need to create the facts in the fact_db
        variableSet targets = ri->targets() ;
	//	Loci::debugout << "targets=" << targets << endl ;
        bool UseRule = true ;
        for(variableSet::const_iterator vi=targets.begin();
            vi!=targets.end();++vi) 
          if(nspace_vec.size()>0 && vi->get_info().namespac != nspace_vec) {
            UseRule = false ;
          }
        if(!UseRule)

          continue ;


        rule_implP rp = ri->get_rule_implP() ;
        for(variableSet::const_iterator vi=targets.begin();
            vi!=targets.end();++vi) {
          // we need to get the storeRep for this variable
          storeRepP srp = rp->get_store(*vi) ;
          if(srp == 0) {
            ostringstream oss ;
            oss << "rule " << *ri << " unable to provide type for " << *vi
                << endl ;
            throw StringError(oss.str()) ;
          }
          // here we only need to set up the variable type in
          // the fact_db
          set_variable_type(*vi,srp) ;
        }
      }
    
      string vname ;
      parse::kill_white_space(s) ;
      if(s.peek()!='{') {
        throw StringError("format error in fact_db::read") ;
        return s ;
      }
      s.get() ;
      
      variableSet read_var_set ;
      for(;;) {
        parse::kill_white_space(s) ;
        if(s.peek() == '}') {
          s.get() ;
	  parse::kill_white_space(s) ;
	  if(s.peek() == '@') { // ok now we are reading in vars for a namespace
	    s.get() ;
	    parse::kill_white_space(s) ;
	    vector<std::string> nspace ;
	    while(parse::is_name(s)) {
	      string name=parse::get_name(s) ;
	      debugout << "name=" << name << endl ;
	      nspace.push_back(name) ;
	      parse::kill_white_space(s) ;
	      if(s.peek() =='@') {
		s.get() ;
		parse::kill_white_space(s) ;
	      }
	    }
	    // Now expecting an open brace for next var space
	    if(s.peek() == '{') {
	      vector<std::string> save = nspace_vec ;
	      nspace_vec = nspace ;
	      read_vars(s,rdb) ;
	      nspace_vec=save ;
	    } else {
	      throw StringError("expecting '{' to follow namespace") ;
	    }
	  }
          break ;
        }
        if(s.peek() == std::char_traits<char>::eof()) {
          throw StringError("unexpected EOF in fact_db::read") ;
        }
        parse::kill_white_space(s) ;
        if(parse::is_name(s)) 
          vname = parse::get_name(s) ;
        else {
          throw StringError("syntax error in fact_db::read") ;
        }
        parse::kill_white_space(s) ;
	if(s.peek() == '@') {
	  vname += '@' ;
	  s.get() ;
	  parse::kill_white_space(s) ;
	  if(parse::is_name(s)) 
	    vname += parse::get_name(s) ;
	  else {
	    ostringstream oss ;
	    oss << "syntax error in fact_db::read, no variable name after namspace '@' separator for namespace '" << vname << "'" ;
	    throw StringError(oss.str()) ;
	  }
	}
        variable var(vname) ;
        try {
          parse::kill_white_space(s) ;
          if(!parse::get_token(s,":")) {
            throw StringError("syntax error in fact_db::read, no ':' separator") ;
          }
          
          variable var(vname) ;
          if(read_var_set.inSet(var)) {
            cerr << "WARNING: Redefining variable '" << var << "' while reading in fact_db!!!!!" << endl ;
          }
          read_var_set += var ;
          storeRepP vp = get_variable(var) ;
          if(vp == 0) {
            vp = get_variable_type(var) ;
            if(vp != 0) {
              create_fact(var,vp) ;
            }
            vp = get_variable(var) ;
          }
          if(vp == 0) {
            ostringstream oss ;
            oss << "variable named '" << vname
                << "' not found in database in fact_db::read." << endl ;
            throw StringError(oss.str()) ;
          }
          vp->Input(s) ;
          parse::kill_white_space(s) ;
          if(s.peek() != '}' && !parse::is_name(s)) {
            ostringstream oss ;
            oss << "syntax error in fact_db::read while reading variable '" << vname << "'" ;
            throw StringError(oss.str()) ;
          }
        } catch(const BasicException &err) {
          err.Print(cerr) ;
          cerr << "input failed for variable " << vname << endl ;
          cerr << "near the following text: \"" ;
          while(s.peek()!=EOF&&s.peek()!='\n'&&s.peek()!='\r')
            cerr << char(s.get()) ;
          cerr << "\"" << endl ;
         
 
          syntax_error = true ;
        }
      }
    } catch(const BasicException &err) {
      err.Print(cerr) ;
      throw StringError("read fact_db failed") ;
    }
    if(syntax_error) {
      cerr << "syntax error reading fact db"  << endl ;
      throw StringError("read fact_db failed") ;
    }
    return s ;
  }

  void fact_db::Print_diagnostics() {
    std::map<variable, fact_info>::iterator mi ;
    ostringstream oss ;
    oss << "memory." ;
    oss << MPI_rank ;
    string file_name = oss.str() ;
    std::ofstream ofile(file_name.c_str(), ios::out) ;
    double total_size = 0 ;
    double total_wasted = 0 ;
    entitySet dom, total, unused ;
    for(mi = fmap.begin(); mi != fmap.end(); ++mi) {
      fact_info &finfo = mi->second ;
      if(isSTORE(finfo.data_rep->getRep())) {
	dom = finfo.data_rep->domain() ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	total_size += finfo.data_rep->pack_size(dom) ; 
	total_wasted += finfo.data_rep->pack_size(unused) ;
      }
    }
    for(mi = fmap.begin(); mi != fmap.end(); ++mi) {
      fact_info &finfo = mi->second ;
      if(isSTORE(finfo.data_rep->getRep())) {
	dom = finfo.data_rep->domain() ;
	double size = finfo.data_rep->pack_size(dom) ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	double wasted_space = finfo.data_rep->pack_size(unused) ;
	ofile << " ****************************************************" << endl ;
	ofile << " Total_size = " << total_size << endl ;
	ofile << "Variable = "  << mi->first << endl ;
	ofile << "Domain = " << dom << endl ;
	ofile << "Size allocated = " << size << endl ; 
	if( isDistributed() )  {
	  Loci::fact_db::distribute_infoP d ;
	  d   = Loci::exec_current_fact_db->get_distribute_info() ;
	  entitySet my_entities = d->my_entities ; 
	  entitySet clone = dom - my_entities ;
	  double clone_size = finfo.data_rep->pack_size(clone) ;
	  ofile << "----------------------------------------------------" << endl;
	  ofile << " My_entities = " << my_entities << endl ;
	  ofile << " Clone entities = " << clone << endl ;
	  ofile << "Memory required for the  clone region  = " << clone_size << endl ;
	  ofile << "Percentage of clone memory required (of size allocated)  = " << double(double(100*clone_size) / size)<< endl ;
	  ofile << "Percentage of clone memory required (of total size allocated)  = " << double(double(100*clone_size) / total_size)<< endl ;
	  ofile << "----------------------------------------------------" << endl;
	}
	ofile << "Percentage of total memory allocated  = " << double(double(100*size) / (total_size+total_wasted)) << endl ;
	ofile << "----------------------------------------------------" << endl;
	ofile << "Total wasted size = " << total_wasted << endl ;
	ofile << "Unused entities = " << unused << endl ;
	ofile << "Wasted space = " << wasted_space << endl ;
	ofile << "Percentage of total memory wasted  = " << double(double(100*wasted_space) / (total_size + total_wasted)) << endl ;
	ofile << " ***************************************************" << endl << endl << endl ;
      }
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////
  
  void reorder_facts(fact_db &facts, dMap &remap) {
    variableSet vars = facts.get_typed_variables() ;
    if(facts.is_distributed_start()) {
      fact_db::distribute_infoP df = facts.get_distribute_info()  ;
      for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
	storeRepP p = facts.get_variable(*vi) ;
	if(p->domain() == ~EMPTY) {
	  // For universal set, keep set universal
	  facts.replace_fact(*vi,p->freeze()) ;
	  p = facts.get_variable(*vi) ;
	} else if(!isMAP(p)) {
	  int kd = p->getDomainKeySpace() ;
	  entitySet rdom = p->domain() & df->g2lv[kd].domain() ;
	  storeRepP fp = (p->remap(df->g2lv[kd]))->freeze() ;
	  facts.replace_fact(*vi,fp) ;
	  p = facts.get_variable(*vi) ;
	} else {
	  MapRepP mp = MapRepP(p->getRep()) ;
	  int kd = p->getDomainKeySpace() ;
	  int rd = mp->getRangeKeySpace() ;
	  storeRepP fp = (mp->MapRemap(df->g2lv[kd],
				       df->g2lv[rd]))->freeze() ;
	  facts.replace_fact(*vi,fp) ;
	  p = facts.get_variable(*vi) ;
	  mp = MapRepP(p->getRep()) ;
	}
      }
#ifdef LOCI_COMPAT_MODE1
      df->g2f.setRep(df->g2fv[0].Rep()) ;
      df->g2l.setRep(df->g2lv[0].Rep()) ;
#endif

    } else {
      for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
	storeRepP p = facts.get_variable(*vi) ;
	if(p->domain() == ~EMPTY) {
	  // For universal set, keep set universal
	  facts.update_fact(*vi,p->freeze()) ;
	} else if(!isMAP(p)) {
	  facts.update_fact(*vi,(p->remap(remap))->freeze()) ;
	} else {
	  MapRepP mp = MapRepP(p->getRep()) ;
	  facts.update_fact(*vi,(mp->MapRemap(remap,remap))->freeze()) ;
	}
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  
  void serial_freeze(fact_db &facts) {
    variableSet vars = facts.get_typed_variables() ;

    fact_db::distribute_infoP df = new fact_db::distribute_info;
    df->myid = 0 ;
    df->isDistributed = 0 ;
    entitySet dom ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP p = facts.get_variable(*vi) ;
      p = p->freeze() ;
      facts.replace_fact(*vi, p) ;
      if(p->domain()!=~EMPTY)
	dom += p->domain() ;
    }
    Map l2g ;
    store<unsigned char> key_domain ;
    l2g.allocate(dom) ;
    key_domain.allocate(dom) ;
    dMap g2l ;
    FORALL(dom,ii) {
      l2g[ii] = ii ;
      key_domain[ii] = 0 ;
      g2l[ii] = ii ;
    } ENDFORALL ;
    df->l2g = l2g.Rep() ;
    df->key_domain = key_domain.Rep() ;
    df->l2f = l2g.Rep() ;
#ifdef LOCI_COMPAT_MODE1
    df->g2l = g2l.Rep() ;
    df->g2f = g2l.Rep() ;
#endif
    df->my_entities = ~EMPTY ;
    df->g2lv.push_back(g2l) ;
    df->g2fv.push_back(g2l) ;
    facts.put_distribute_info(df) ;
  }

  void fact_db::set_variable_type(variable v, storeRepP st) {
    v = add_namespace(v) ;
    tmap[v] = storeRepP(st->new_store(EMPTY)) ;
  }

  storeRepP fact_db::get_variable_type(variable v) const {
    v = add_namespace(v) ;
    map<variable,storeRepP>::const_iterator mi ;
    if((mi=tmap.find(v)) != tmap.end())
      return storeRepP(mi->second->new_store(EMPTY)) ;
    else
      return storeRepP(0) ;
  }
  void fact_db::write_all_hdf5(const char *filename) {
    variableSet vars = get_typed_variables() ;
    write_hdf5(filename, vars) ;
  }
  void fact_db::read_all_hdf5(const char *filename) {
    variableSet vars = get_typed_variables() ;
    read_hdf5(filename, vars) ; 
  }

  void fact_db::write_hdf5(const char *filename, variableSet &vars) {
    hid_t  file_id=0 ;
    file_id = Loci::hdf5CreateFile(filename, H5F_ACC_TRUNC,
                                   H5P_DEFAULT, H5P_DEFAULT) ;
       
    for(variableSet::const_iterator vi = vars.begin(); vi != vars.end(); ++vi) {
      storeRepP  p = get_variable(*vi) ;
      if(isSTORE(p)) {
        writeContainer(file_id,variable(*vi).get_info().name,p,*this) ;
      }
    }
    hdf5CloseFile(file_id);
  }
  
  
 

  void fact_db::read_hdf5(const char *filename, variableSet &vars) {
    hid_t  file_id=0 ;
    file_id =  hdf5OpenFile(filename,  H5F_ACC_RDONLY, H5P_DEFAULT);
    for(variableSet::const_iterator vi = vars.begin(); vi != vars.end(); ++vi) {
      storeRepP  p = get_variable(*vi) ;
      if(isSTORE(p)) {
        readContainer(file_id,variable(*vi).get_info().name,p,EMPTY,*this) ;
      }
    }
    hdf5CloseFile(file_id);
  }


#ifdef DYNAMICSCHEDULING
  // experimental code to create keyspace from the
  // global registered keyspace list
  // returns "true" to indicate the methods succeeded,
  // "false" to indicate an error.
  bool
  fact_db::create_keyspace(KeySpaceList& global_list) {
    for(KeySpaceList::Iterator ki=global_list.begin();
        ki!=global_list.end();++ki) {
      KeySpaceP kp = ki.get_p()->rr->get_space() ;
      // first get the space name
      if(!kp->named_space()) {
        if(Loci::MPI_rank == 0) {
	  auto &kpdr = *kp ;
          cerr << "fact_db Error: Initializing Unnamed Keyspace!"
               << " typeid = " << typeid(kpdr).name() << endl ;
	}
        return false ;
      }
      string name = kp->get_name() ;
      map<string,KeySpaceP>::const_iterator mi = keyspace.find(name) ;
      if(mi!=keyspace.end()) {
        if(Loci::MPI_rank == 0)
          cerr << "fact_db Error: Duplicated Keyspace: "
               << name << endl ;
        return false ;
      }
      kp->set_synonyms(&synonyms) ;
      keyspace[name] = kp ;
    }
    return true ;
  }

  KeySpaceP
  fact_db::get_keyspace(const string& kname) {
    map<string,KeySpaceP>::const_iterator mi = keyspace.find(kname) ;
    if(mi == keyspace.end())
      return KeySpaceP(0) ;
    else
      return mi->second ;
  }
  
  void
  fact_db::init_key_manager() {
    int max_alloc = get_max_alloc(0) ;
    int global_max = 0 ;
    MPI_Allreduce(&max_alloc, &global_max, 1,
                  MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
    key_manager = new KeyManager(global_max+1) ;
  }
#endif  
}

  
