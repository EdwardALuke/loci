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
#include <keyspace.h>
#include <key_manager.h>
#include <iostream>
#include <limits>

using std::cout ;
using std::cerr ;
using std::endl ;
using std::vector ;
using std::string ;

#include <sstream>
using std::ostringstream ;
#include <map>
using std::map ;
#include <set>
using std::set ;

#include <parSampleSort.h>      // for balanceDistribution()
#include <distribute.h>
#include <Tools/except.h>

#include <execute.h>
#include <rule.h>
#include "sched_tools.h"

#ifdef DYNAMICSCHEDULING
namespace Loci {
  extern void ORBPartition(const std::vector<vector3d<float> >&,
                           std::vector<int>&, MPI_Comm) ;

  // utility functions
  namespace {
    // this function partitions keys according to their
    // spatial coordinates, returns a new KeySet for each process
    KeySet
    orb_partition_keys(vector<vector3d<float> >& coords,
                       vector<Key>& keys, MPI_Comm comm) {
      FATAL(coords.size() != keys.size()) ;

      KeySet new_keys = Loci::create_intervalSet(keys.begin(),keys.end()) ;
      
      unsigned int len = coords.size() ;
      unsigned int max_len = 0 ;
      MPI_Allreduce(&len, &max_len, 1, MPI_UNSIGNED, MPI_MAX, comm) ;

      if(max_len == 0)
        return new_keys ;       // no elements at all

      int np = 0 ;
      MPI_Comm_size(comm, &np) ;

      unsigned int min_len = 0 ;
      MPI_Allreduce(&len, &min_len, 1, MPI_UNSIGNED, MPI_MIN, comm) ;

      if(min_len < static_cast<unsigned int>(np)) {
        // if some processes have less than "np" items,
        // we would need to perform a balance first
        balanceDistribution(coords, comm) ;
        balanceDistribution(keys, comm) ;
        FATAL(coords.size() != keys.size()) ;

        new_keys = Loci::create_intervalSet(keys.begin(),keys.end()) ;

        len = coords.size() ;
        // check to see if the minimum sequence length on
        // any process is still less than "np" or not.
        // in case of yes, just return without doing
        // anything, since we cannot perform the orb
        // partition because it requires that the items
        // on any process must be >= "np"
        MPI_Allreduce(&len, &min_len, 1, MPI_UNSIGNED, MPI_MIN, comm) ;

        if(min_len < static_cast<unsigned int>(np))
          return new_keys ;
      }
      // we then proceed to use orb partition to determine the distribution
      vector<int> i2p(len) ;
      Loci::ORBPartition(coords, i2p, comm) ;

      // we now have the item to process mapping, we need
      // to package and communicate them (but we only need to
      // communicate the "keys" vector since we are only computing
      // key distributions)

      vector<int> send_counts(np, 0) ;
      for(size_t i=0;i<len;++i)
        send_counts[i2p[i]] += 1 ;

      vector<int> recv_counts(np, 0) ;
      MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                   &recv_counts[0], 1, MPI_INT, comm) ;

      vector<int> send_displs(np, 0) ;
      vector<int> recv_displs(np, 0) ;
      for(int i=1;i<np;++i) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
      }

      int tot_send_counts = 0 ;
      int tot_recv_counts = 0 ;
      for(int i=0;i<np;++i) {
        tot_send_counts += send_counts[i] ;
        tot_recv_counts += recv_counts[i] ;
      }
      // create and pack the send buffer
      vector<Key> send_buffer(tot_send_counts) ;
      vector<int> buffer_idx = send_displs ;
      for(size_t i=0;i<len;++i) {
        int pid = i2p[i] ;
        int& offset = buffer_idx[pid] ;
        send_buffer[offset] = keys[i] ;
        ++offset ;
      }

      vector<Key> recv_buffer(tot_recv_counts) ;

      for(int i=0;i<np;++i) {
        send_counts[i] *= sizeof(Key) ;
        send_displs[i] *= sizeof(Key) ;
        recv_counts[i] *= sizeof(Key) ;
        recv_displs[i] *= sizeof(Key) ;
      }
      
      MPI_Alltoallv(&send_buffer[0], &send_counts[0],
                    &send_displs[0], MPI_BYTE,
                    &recv_buffer[0], &recv_counts[0],
                    &recv_displs[0], MPI_BYTE, comm) ;

      return Loci::create_intervalSet(recv_buffer.begin(),
                                      recv_buffer.end()) ;
    }

    // this function gathers the keyset on each process and
    // form a distribution vector
    vector<KeySet>
    gather_keyset(const KeySet& keys, MPI_Comm comm) {
      int np = 0 ;
      MPI_Comm_size(comm, &np) ;

      bool pack_interval = false ;
      int elem_size = keys.size() ;
      int interval_size = keys.num_intervals() ;
      pack_interval = (2*interval_size) < elem_size ;

      int local_size = 0 ;
      if(keys != EMPTY) {
        if(pack_interval)
          local_size = 2*interval_size + 1 ;
        else
          local_size = elem_size + 1 ;
      }

      vector<int> recv_counts(np,0) ;
      MPI_Allgather(&local_size, 1, MPI_INT,
                    &recv_counts[0], 1, MPI_INT, comm) ;
      vector<int> recv_displs(np,0) ;
      for(int i=1;i<np;++i)
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

      int tot_size = recv_displs[np-1] + recv_counts[np-1] ;

      // pack the local keys into an array
      vector<int> local_keys(local_size) ;
      if(local_size > 0) {
        if(pack_interval) {
          local_keys[0] = 1 ;   // signals we're packing intervals
          int count = 1 ;
          for(size_t i=0;i<keys.num_intervals();++i) {
            local_keys[count++] = keys[i].first ;
            local_keys[count++] = keys[i].second ;
          }
        } else {
          local_keys[0] = 0 ;   // signals we're packing elements
          int count = 1 ;
          for(KeySet::const_iterator ki=keys.begin();
              ki!=keys.end();++ki,++count)
            local_keys[count] = *ki ;
        }
      }

      // allocate the entire array for all data from all processors
      vector<int> global_keys(tot_size) ;
      MPI_Allgatherv(&local_keys[0], local_size, MPI_INT,
                     &global_keys[0], &recv_counts[0],
                     &recv_displs[0], MPI_INT, comm) ;
      // release useless buffer
      vector<int>().swap(local_keys) ;
      // unpack the raw buffer into a vector<KeySet>
      vector<KeySet> out(np) ;

      for(int i=0;i<np;++i) {
        int b = recv_displs[i] ;
        int e = b + recv_counts[i] ;
        if(b==e) {
          out[i] = EMPTY ;
          continue ;
        }
        KeySet& k = out[i] ;
        int flag = global_keys[b] ;
        ++b ;
        if(flag == 1) {
          // unpack intervals
          for(;b<e;b+=2) {
            int f = global_keys[b] ;
            int s = global_keys[b+1] ;
            k += interval(f,s) ;
          }
        } else if(flag == 0) {
          // unpack elements
          for(;b<e;++b)
            k += global_keys[b] ;
        } else {
          // impossible unless something is wrong
          int p ;
          MPI_Comm_rank(comm, &p) ;
          cerr << "Function gather_keyset failed on process "
               << p << ", with unpack flag = " << flag << endl ;
          Loci::Abort() ;
        }
      }

      return out ;
      
    } // end gather keys
    
  } // end of un-named namespace

  vector<KeySet>
  OrbKeySpace::orb_partition() {
    KeySet new_keys = orb_partition_keys(coords,coord_keys,comm) ;
    return gather_keyset(new_keys,comm) ;
  }

  void
  KeySpace::create_shadow(const variable& v, storeRepP sp) {
    variable true_v = remove_synonym(v) ;
    std::map<variable,store_refP>::iterator
      mi = shadow.find(true_v) ;
    if(mi != shadow.end()) {
      if(typeid(sp->getRep()) != typeid(mi->second->getRep()))
        std::cerr << "Retyping variable " << true_v
                  << " in the shadow of keyspace " << space_name << endl ;
      mi->second->setRep(sp->getRep()) ;
    } else {
      store_refP s = new store_ref ;
      s->setRep(sp->getRep()) ;
      shadow[true_v] = s ;
    }
  }

  // get the rep information of a shadow variable
  storeRepP
  KeySpace::get_shadow(const variable& v) const {
    variable true_v = remove_synonym(v) ;

    std::map<variable,store_refP>::const_iterator
      mi = shadow.find(true_v) ;
    if(mi == shadow.end()) {
      return storeRepP(0) ;
    } else {
      return storeRepP(mi->second) ;
    }
  }
  // NOTE: we would need to come back and revise these
  // methods. If a keyspace uses local numbering, then
  // we also need to convert the added new keys into
  // a local number system and also expand the g2l and
  // l2g mapping. BUT! again I think only dynamic keyspace
  // is ever going to be allowed to have insertion and
  // deletion rules and we only use the global numbering
  // scheme for dynamic keyspace. So it might be okay to
  // ignore the local numbering completely here.
  void KeySpace::
  add_keys(const KeySet& ks) {
    keys += ks ;
    my_keys += ks ;
    // we need to recompute the key partition by
    // re-gathering all keys on each process
    //key_ptn = gather_keyset(keys, comm) ;
  }

  void KeySpace::
  remove_keys(const KeySet& ks) {
    keys -= ks ;
    my_keys -= ks ;
    //key_ptn = gather_keyset(keys, comm) ;
  }
  
  // the following codes provide facilities to set up user defined
  // KeySpaces. the codes mirror the currently design and implementation
  // of rule list and the register rule mechanism.

  KeySpaceP
  KeySpace::new_keyspace() const {
    cerr << "KeySpace::new_keyspace() was called:" << endl
         << "this method should never be called, use CopyKeySpace<> to"
         << endl << " generate a copyable key space" << endl ;
    abort() ;
    return 0 ;
  }
  // KeySpaceP
  // KeySpace::clone_keyspace() const {
  //   cerr << "KeySpace::clone_keyspace() was called:" << endl
  //        << "this method should never be called, use CopyKeySpace<> to"
  //        << endl << " generate a copyable key space" << endl ;
  //   abort() ;
  //   return 0 ;
  // }

  void KeySpace::
  set_dcontrol_map(const variable& v, const rule& r) {
    ostringstream oss ;
    oss << r ;
    dcontrol_map[v].insert(oss.str()) ;
  }

  void KeySpace::
  register_dcontrol_var(const variable& v, execute_modules* ep) {
    dcontrol_vars[v] = ep ;
  }

  void KeySpace::
  register_dcontrol_rule(const rule& r, execute_modules* ep) {
    ostringstream oss ;
    oss << r ;
    dcontrol_rules[oss.str()] = ep ;
  }

  void KeySpace::
  set_dcontrol() const {
#ifdef DYNAMICSCHEDULING
    for(map<variable,set<string> >::const_iterator
          mi=dcontrol_map.begin();mi!=dcontrol_map.end();++mi) {
      const variable& var = mi->first ;
      // get the execute_dcontrol_reset* first
      map<variable, execute_modules*>::const_iterator
        var_epi = dcontrol_vars.find(var) ;
      if(var_epi == dcontrol_vars.end()) {
        continue ;
        // if not found, then we skip to next var
      }
      execute_dcontrol_reset* dcrep =
        dynamic_cast<execute_dcontrol_reset*>(var_epi->second) ;
      if(dcrep == 0) {
        cerr << "Error: dynamic cast execute_dcontrol_reset* failed!"
             << " in keyspace: " << get_name() << " for variable: "
             << var << endl ;
        debugout << "Error: dynamic cast execute_dcontrol_reset* failed!"
                 << " in keyspace: " << get_name() << " for variable: "
                 << var << endl ;
        debugout.flush() ;
        Loci::Abort() ;
      }
      const set<string>& drules = mi->second ;
      for(set<string>::const_iterator si=drules.begin();
          si!=drules.end();++si) {
        // search through the dcontrol_rule records
        map<string, execute_modules*>::const_iterator
          rule_epi = dcontrol_rules.find(*si) ;
        if(rule_epi == dcontrol_rules.end()) {
          cerr << "Error: cannot find execute_drule*"
               << " in keyspace: " << get_name() << " for rule: "
               << *si << endl ;
          debugout << "Error: cannot find execute_drule*"
                   << " in keyspace: " << get_name() << " for rule: "
                   << *si << endl ;
          debugout.flush() ;
          Loci::Abort() ;
        }
        execute_drule_module* drm =
          dynamic_cast<execute_drule_module*>(rule_epi->second) ;
        if(drm == 0) {
          cerr << "Error: dynamic cast execute_drule* failed!"
               << " in keyspace: " << get_name() << " for rule: "
               << *si << endl ;
          debugout << "Error: dynamic cast execute_drule* failed!"
                   << " in keyspace: " << get_name() << " for rule: "
                   << *si << endl ;
          debugout.flush() ;
          Loci::Abort() ;
        }
        // everything okay, we just add the drule_module
        dcrep->add_drule(drm) ;
      } // end for(drules)
    }
#else
    Loci::Abort() ;
#endif
  }

  // definition of global key space list
  RegisterKeySpaceList register_key_space_list ;
  KeySpaceList global_key_space_list ;

  KeySpaceList::~KeySpaceList() {
    KeySpaceListEnt *p, *v ;
    for(p=list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
  }

  void
  KeySpaceList::clear() {
    KeySpaceListEnt *p, *v ;
    for(p=list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
    list = 0 ;
  }

  void
  KeySpaceList::push_space(RegisterKeySpaceType* p) {
    KeySpaceListEnt* flp = new KeySpaceListEnt(p,list) ;
    list = flp ;
  }

  void
  KeySpaceList::copy_space_list(const KeySpaceList& sl) {
    KeySpaceListEnt *p, *v ;
    for(p=sl.list;p!=0;p=v) {
      push_space(p->rr) ;
      v = p->next ;
    }
  }

  void
  KeySpaceList::copy_space_list(const RegisterKeySpaceList& sl) {
    KeySpaceListEnt *p, *v ;
    for(p=sl.global_list;p!=0;p=v) {
      push_space(p->rr) ;
      v = p->next ;
    }
  }

  // declaration of static variable global_list
  KeySpaceList::KeySpaceListEnt* RegisterKeySpaceList::global_list = 0 ;

  RegisterKeySpaceList::~RegisterKeySpaceList() {
    KeySpaceListEnt *p, *v ;
    for(p=global_list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
  }

  void
  RegisterKeySpaceList::clear() {
    KeySpaceListEnt *p, *v ;
    for(p=global_list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
    global_list = 0 ;
  }

  bool
  RegisterKeySpaceList::empty() {
    return (global_list == 0) ;
  }

  void
  RegisterKeySpaceList::push_space(RegisterKeySpaceType* p) {
    KeySpaceListEnt* flp = new KeySpaceListEnt(p,global_list) ;
    global_list = flp ;
  }

  // KeyManager implementation
  KeyManager::KeyManager(Key global_start) {
    // first we obtain the limits on the Key type
    // we would like it to be slightly smaller than
    // the true max value
    Key global_max_key = std::numeric_limits<Key>::max() - 5 ;
    Key global_range = global_max_key - global_start ;
    // get the number of processes
    int comm_size = 0 ;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size) ;
    Key local_range = global_range / comm_size ;

    rank = 0 ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    min_key = rank * local_range + global_start ;
    max_key = min_key + local_range - 1 ;

    int np = 0 ;
    MPI_Comm_size(MPI_COMM_WORLD, &np) ;
    min_dist.resize(np) ;
    max_dist.resize(np) ;
    MPI_Allgather(&min_key, 1, MPI_INT,
                  &min_dist[0], 1, MPI_INT, MPI_COMM_WORLD) ;
    MPI_Allgather(&max_key, 1, MPI_INT,
                  &max_dist[0], 1, MPI_INT, MPI_COMM_WORLD) ;
    range_size = max_key - min_key + 1 ;

    // initially all are free keys
    freeset = interval(min_key, max_key) ;
  }

  Key
  KeyManager::generate_key() {
    if(freeset == EMPTY) {
      std::stringstream err ;
      err << "KeyManager Error: process " << rank
          << " is running out of keys!" ;
      throw StringError(err.str()) ;
    }
    Key k = freeset.Min() ;
    freeset -= k ;

    return k ;
  }

  KeySet
  KeyManager::generate_key(size_t size) {
    if( (size_t)freeset.size() < size) {
      std::stringstream err ;
      err << "KeyManager Error: process " << rank
          << " is running out of keys!" ;
      throw StringError(err.str()) ;
    }
    KeySet k ;
    for(KeySet::const_iterator i=freeset.begin();
        size!=0;++i,--size)
      k += *i ;

    freeset -= k ;

    return k ;
  }

  void
  KeyManager::recycle_key(Key k) {
    freeset += k ;
  }

  void
  KeyManager::recycle_key(const KeySet& keys) {
    freeset += keys ;
  }

  KeyManagerP
  KeyManager::clone() const {
    return new KeyManager(*this) ;
  }

  void
  KeyManager::compact(const std::vector<KeySet>& key_ptn,
                      dMap& remap, std::vector<KeySet>& new_key_ptn) {
    new_key_ptn.resize(key_ptn.size()) ;
    // since we know the distribution of min and max on each
    // process, we can renumber all of the key partitions
    // without communication. we just assume every process all
    // starts from the min in the renumber
    for(size_t i=0;i<key_ptn.size();++i) {
      // check range on local process
      KeyRangeType cur_range_size = key_ptn[i].size() ;
      if(i == (size_t)rank) {
        if(cur_range_size > range_size) {
          std::stringstream err ;
          err << "KeyManager Error: process " << rank
              << " insufficient range to hold all current keys!" ;
          throw StringError(err.str()) ;
        } else {
          freeset = interval(min_key+cur_range_size, max_key) ;
        }
      }
      const KeySet& ks = key_ptn[i] ;
      Key f = min_dist[i] ;
      for(KeySet::const_iterator j=ks.begin();j!=ks.end();++j,++f) {
        remap[*j] = f ;
      }
      // get the new_key_ptn[i] (min_dist[i], f-1)
      if(cur_range_size > 0)
        new_key_ptn[i] =
          interval(min_dist[i], min_dist[i]+cur_range_size-1) ;
      else
        new_key_ptn[i] = EMPTY ;
    }
  }

} // end of namespace Loci

#endif
// ... end of file ...
