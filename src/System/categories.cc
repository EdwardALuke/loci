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
#include "dist_tools.h"
#include <Tools/debug.h>
#include <entitySet.h>
#include "Tools/debugger.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <map>
using std::map ;
#include <list>
using std::list ;
#include <string>
using std::string ;

#include <iostream>
using std::cout ;
using std::cerr ;
using std::endl ;
using std::ios ;
using std::ostream ;

#include <algorithm>
using std::sort ;
using std::pair ;
using std::make_pair ;

using std::ostringstream ;
//#define VERBOSE


namespace Loci {

  inline bool compare_var_sizes(const pair<variable,int> &p1, const pair<variable,int> &p2) {
    return p1.second > p2.second ;
  }

  inline bool compare_groups(const pair<variableSet,string> &p1,
                             const pair<variableSet,string> &p2) {
    return p1.second < p2.second ;
  }


  namespace {
    // Compute associations between variables and their domains
    // Note, we also create image variables to represent the images of maps
    //   over their domains.
    void getVariableAssociations(map<variable,entitySet> &vm, fact_db &facts, int kd) {
      vector<entitySet> ptn = facts.get_init_ptn(kd) ;
      entitySet total_entities ;
      entitySet my_entities = ptn[MPI_rank] ;
      vm.clear() ;
      variableSet vars = facts.get_typed_variables() ;
      for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
        storeRepP p = facts.get_variable(*vi) ;

        if(p->getDomainKeySpace() == kd) {
          if(isMAP(p)) {
            entitySet tmp = my_entities & p->domain() ;
            vm[*vi] = tmp ;
            total_entities += tmp ;
          } else if(isSTORE(p)) {
            entitySet tmp = my_entities & p->domain() ;
            vm[*vi] = tmp ;
            total_entities += tmp ;
          } else {
            if(p->domain() != ~EMPTY) {
              entitySet tmp = my_entities & p->domain() ;
              vm[*vi] = tmp ;
              total_entities += tmp ;
            }
          }
        }
      }
      for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
        storeRepP p = facts.get_variable(*vi) ;
        if(isMAP(p)) {
          // Add any map image that refers to entities not already identified.
          MapRepP mp = MapRepP(p->getRep()) ;

          if(mp->getRangeKeySpace() == kd) {
            entitySet dom = my_entities & p->domain() ;
            entitySet image_dom = mp->image(dom) ;
            image_dom += distribute_entitySet(image_dom, ptn) ;
            entitySet testSet = image_dom - total_entities ;
            int size = testSet.size() ;
            int tsize = 0 ;
            MPI_Allreduce(&size, &tsize, 1, MPI_INT,
                          MPI_SUM, facts.get_comm()) ;
            if(tsize != 0) {
              std::string name = "image_" ;
              name.append(vi->get_info().name) ;
              variable v = variable(name) ;
              vm[v] = image_dom ;
              total_entities += image_dom ;
            }
          }
        }
      }
      variable v_leftout("SENTINEL_LEFTOUT") ;
      vm[v_leftout] = ptn[MPI_rank]-total_entities ;
    }

    // Compute the distict intervals that are required to describe all of
    // the variable domain sets.
    void getLocalIntervals(vector<interval> &pvec,
                           const map<variable,entitySet> &vm) {
      pvec.clear() ;
      vector<int> vals ;
      set<entitySet>::iterator si ;
      map<variable,entitySet>::const_iterator mi ;
      entitySet totSet ;
      for(mi=vm.begin();mi!=vm.end();++mi) {
        entitySet s = mi->second ;
        totSet += s ;
        for(size_t i = 0;i < s.num_intervals(); ++i) {
          vals.push_back(s[i].first-1) ;
          vals.push_back(s[i].first) ;
          vals.push_back(s[i].second) ;
          vals.push_back(s[i].second+1) ;
        }
      }

      if(vals.size() > 0) {
        std::sort(vals.begin(),vals.end()) ;

        vector<int>::iterator uend = std::unique(vals.begin(),vals.end()) ;
        
        vector<int>::iterator ii = vals.begin() ;
        
        ++ii ;
        while(ii+1 != uend) {
          int i1 = *ii ;
          int i2 = *(ii+1)  ;
          if(i1+1 == i2)
            i2 = i1 ;
          else
            ++ii ;
          ++ii ;
          interval iv(i1,i2) ;
          if((entitySet(iv)&totSet) != EMPTY)
            pvec.push_back(iv) ;
        }
      }
    }

    // Identify categories by the association of attributes (variables)
    // with entities.  These categories are only for entities owned by
    // the current processors (i.e. local)
    void getLocalCategories(map<variableSet, entitySet> &cmap,
                            const map<variable, entitySet> &vm) {
      cmap.clear() ;
      vector<interval> pvec ;
      getLocalIntervals(pvec,vm) ;

      vector<variableSet> cat_names(pvec.size()) ;
#ifdef VERBOSE
      debugout << "pvec.size() = " << pvec.size()

               << "vm.size() = " << vm.size() << endl ;
#endif

      map<entitySet, variableSet> vt ;
      map<variable,entitySet>::const_iterator mi ;
      for(mi=vm.begin();mi!=vm.end();++mi) {
        vt[mi->second] += mi->first ;
      } ;
#ifdef VERBOSE
      debugout << "pvec.size() = " << pvec.size()
               << " vm.size() = " << vm.size()
               << " vt.size() = " << vt.size() << endl ;
#endif
      map<entitySet,variableSet>::const_iterator mt ;

      for(mt=vt.begin();mt!=vt.end();++mt) {
        variableSet vs = mt->second ;
        entitySet e = mt->first ;
        for(size_t i=0;i<pvec.size();++i) {
          if(e.inSet(pvec[i].first)) {
#ifdef DEBUG
            entitySet isect = entitySet(pvec[i]) & e ;
            WARN(isect != entitySet(pvec[i]) ) ;
#endif
            cat_names[i] += vs ;
          }
        }
      }

      for(size_t i=0;i<pvec.size();++i)
        cmap[cat_names[i]] += pvec[i] ;

    }



    // Unify local categories across all processors, return in a universal order
    // variables and categories.
    void getGlobalCategories(vector<variableSet> &clist,
                             vector<variable> &vlist,
                             const map<variableSet,entitySet> &cmap,
                             const map<variable,entitySet> &vm) {
      vector<string> vnlist ;
      map<variable,entitySet>::const_iterator ii ;
      for(ii=vm.begin();ii!=vm.end();++ii) {
        ostringstream oss ;
        oss << ii->first ;
        vnlist.push_back(oss.str()) ;
      }
      sort(vnlist.begin(),vnlist.end()) ;
      vlist.clear() ;
      map<variable,int> vident ;
      for(size_t i=0;i<vnlist.size();++i) {
        variable v(vnlist[i]) ;
        vlist.push_back(v) ;
        vident[v] = i ;
      }

      set<entitySet> cat_set ;
      map<variableSet,entitySet>::const_iterator ci ;
      for(ci=cmap.begin();ci!=cmap.end();++ci) {
        entitySet tmp ;
        variableSet s = ci->first ;
        for(variableSet::const_iterator vi=s.begin();vi!=s.end();++vi)
          tmp += vident[*vi] ;
        cat_set.insert(tmp) ;
      }

      set<entitySet> cat_set_unified ;
      int comm_root = 0 ;
      while(comm_root != -1) {
        vector<int> buffer ;
        if(comm_root == MPI_rank) {
          set<entitySet>::const_iterator si ;
          for(si=cat_set.begin();si!=cat_set.end();++si) {
            entitySet s = *si ;
            buffer.push_back(s.size()) ;
            for(entitySet::const_iterator ei=s.begin();ei!=s.end();++ei)
              buffer.push_back(*ei) ;
          }
        }
        int buf_size = buffer.size() ;
        MPI_Comm comm = get_exec_comm() ;
        MPI_Bcast(&buf_size, 1, MPI_INT, comm_root, comm) ;
        if(comm_root != MPI_rank)
          buffer = vector<int>(buf_size) ;
        MPI_Bcast(&buffer[0], buf_size, MPI_INT, comm_root,
                  comm) ;

        for(size_t ex=0;ex<buffer.size();) {
          int num_ents = buffer[ex] ;
          ex++ ;
          entitySet s ;
          for(int i=0;i<num_ents;++i) {
            s += buffer[ex] ;
            ex++ ;
          }
          cat_set_unified.insert(s) ;
          cat_set.erase(s) ;
        }
        int r = -1 ;
        if(cat_set.size() != 0)
          r = MPI_rank ;
        MPI_Allreduce(&r, &comm_root, 1, MPI_INT, MPI_MAX,
                      comm) ;
      }
      clist.clear() ;
      set<entitySet>::const_iterator si ;
      for(si=cat_set_unified.begin();si!=cat_set_unified.end();++si) {
        entitySet s = *si ;
        variableSet vs ;
        for(entitySet::const_iterator ei=s.begin();ei!=s.end();++ei)
          vs += vlist[*ei] ;
        clist.push_back(vs) ;
      }
    }

    // Collect the sizes of the distributed variable sets
    void getVarSizes(vector<pair<variable,int> > &var_sizes,
                     const vector<variable> &var_list,
                     const map<variable,entitySet> &vm) {
      var_sizes.clear() ;
      for(size_t i=0;i<var_list.size();++i) {
        variable v = var_list[i] ;
        map<variable,entitySet>::const_iterator vinfo = vm.find(v) ;
        int local_size = vinfo->second.size() ;
        int vsize = 0 ;
        MPI_Allreduce(&local_size, &vsize, 1, MPI_INT, MPI_SUM,
                      get_exec_comm()) ;
        var_sizes.push_back(pair<variable,int>(v,vsize)) ;
      }
      std::stable_sort(var_sizes.begin(),var_sizes.end(),compare_var_sizes) ;
    }

    void getCategoryKeys(vector<pair<variableSet,string> > &cat_keys,
                         const vector<pair<variable,int> > &var_sizes,
                         const vector<variableSet> &cat_list) {
      cat_keys.clear() ;
      for(size_t i=0;i<cat_list.size();++i) {
        variableSet vset = cat_list[i] ;
        string key ;
        for(size_t j=0;j<var_sizes.size();++j) {
          key += vset.inSet(var_sizes[j].first)?'1':'0' ;
        }
        cat_keys.push_back(pair<variableSet,string>(vset,key)) ;
      }
      std::stable_sort(cat_keys.begin(),cat_keys.end(),compare_groups) ;
    }
  }
  void categories(fact_db &facts,vector<entitySet> &pvec, int kd) {
    map<variable, entitySet> vm ;
    getVariableAssociations(vm,facts,kd) ;

    if(vm.size() == 1) // Empty kd
      return ;

    map<variableSet, entitySet> cmap ;
    getLocalCategories(cmap,vm) ;

    vector<variableSet> cat_list ;
    vector<variable> var_list ;
    getGlobalCategories(cat_list,var_list,cmap,vm) ;

    vector<pair<variable,int> > var_sizes ;
    getVarSizes(var_sizes,var_list,vm) ;

    vector<pair<variableSet,string> > cat_keys ;
    getCategoryKeys(cat_keys,var_sizes,cat_list) ;

    pvec.clear() ;
#ifdef VERBOSE
    debugout << "num_categories = " << cat_list.size() <<endl ;
#endif
    for(size_t i=0;i<cat_keys.size();++i) {
      variableSet vset = cat_keys[i].first ;
      entitySet eset = cmap[vset] ;
      pvec.push_back(eset) ;
#ifdef VERBOSE
      debugout << "category " << cat_keys[i].first << " = " << cat_keys[i].second << ", pvec_n = " << pvec[i].num_intervals() <<  endl ;
#endif

    }

  }
}
