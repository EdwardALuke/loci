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
#ifndef KEY_MANAGER_H
#define KEY_MANAGER_H

#include <vector>
#include <Tools/intervalSet.h>
#include <Tools/cptr.h>
#include <DMap.h>

namespace Loci {

  typedef Entity Key ;
  typedef entitySet KeySet ;

  // here we define a key manager that will be used to generate and
  // (possibly) recycle keys used during the program run.
  // In the current design, we try to keep a single global number
  // scheme, and also in order to store the key partitions as compact
  // as possible, we pre-allocate the segment of keys that can be
  // generated on each process. for example, when using 32-bit integers,
  // we can have total of 2^32 keys, if we have p processes, we will
  // pre-allocate 2^32/p keys on each process. this way, each process
  // can have a compact representation of the keys in use. When a keyspace
  // is redistributed, we also remap the keys so that they again fit
  // in the segment pre-allocated on each process. also if on some
  // processes, keys are generated and recycled frequently, "holes"
  // may be created in the key range segment, we may also want to
  // periodically compact the key range on each process. Note, in the
  // current design, if a process runs out of keys, we simply halt
  // the program entirely. a more sophisticated way is to "borrow"
  // keys from other processes (assuming there are others that have
  // not run out of keys). this could be a costly operation and in
  // this prototype version, we are not doing it, but it may be
  // considered later. Also if there are multiple keyspaces defined,
  // we will also need to try to keep keys as compact as possible
  // within a keyspace. therefore we may also want to partition the
  // local key range into different segments for different keyspaces.
  // this is to be evaluated later.
  typedef int KeyRangeType ;
  class KeyManager: public CPTR_type {
    Key min_key ;              // the base key on this process
    Key max_key ;              // the maximum key id allowed on this process
    std::vector<Key> min_dist ; // distribution of min on each process
    std::vector<Key> max_dist ; // distribution of max on each process
    KeyRangeType range_size ;   // total keys can be used on this process
    // thus the interval [min_key, max_key] defines the range
    // of keys that are possible on this local process.
    KeySet freeset ;            // the keys that are not used
    int rank ;
  public:
    typedef CPTR<KeyManager> KeyManagerP ;
    //KeyManager(Key min, Key max):min_key(min), max_key(max) {}
    // this constructor takes a starting key and then builds
    // the min_key and max_key on each process. i.e., the "start"
    // key is the global minimum allowed key.
    KeyManager(Key global_start) ;
    // generate a single key
    Key
    generate_key() ;
    // generate a set of keys, whose size is passed in as an argument
    KeySet
    generate_key(size_t size) ;
    // recycle methods
    void
    recycle_key(Key k) ;
    void
    recycle_key(const KeySet& keys) ;

    KeyManagerP
    clone() const ;

    // this method takes a partition of KeySet and then tries to
    // compact the entire keys covered it returns a remap that
    // renumbers the passed in KeySet. This is used mainly to
    // compact the keyspace keys after many key generation,
    // recycles, and redistribution, it also generates a new
    // key_ptn
    void
    compact(const std::vector<KeySet>& key_ptn,
            dMap& remap, std::vector<KeySet>& new_key_ptn) ;
  } ;
  typedef KeyManager::KeyManagerP KeyManagerP ;

} // end of namespace Loci

#endif
