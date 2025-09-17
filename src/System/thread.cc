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

// disable this file by default for now
#ifdef PTHREADS

#include "thread.h"
#include <sstream>
using std::stringstream;

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::set;
using std::pair;
using std::istream;
using std::ostream;
using std::deque;

namespace Loci {

  ThreadControl* thread_control = 0;
  pthread_spinlock_t** ReductionLock::locks = 0;

  // this implements the leading execution predicate
  bool is_leading_execution()
  {
    bool lp = MPI_processes>1 ? MPI_rank==0 : true;
    bool lt = thread_control ? thread_control->is_leading_thread() : true;
    return lp && lt;
  }

  // these two functions are intended to provide a global atomic
  // region of execution (the ideal use of them is to use the 
  // Loci preprocessor to hide these calls from the users)
  void global_atomic_region_begin()
  {
    if (thread_control)
      thread_control->atomic_begin();
  }

  void global_atomic_region_end()
  {
    if (thread_control)
      thread_control->atomic_end();
  }
  
  // this function generates a unique name on each process/thread
  std::string gen_local_name(const string& prefix,
                             const string& suffix)
  {
    string mpi_id = "";
    string thread_id = "";
    if (MPI_processes > 1) {
      stringstream ss; ss << MPI_rank;
      mpi_id += ss.str();
    }
    if (thread_control) {
      thread_id += thread_control->get_local_id();
    }
    string name = prefix;
    if (!mpi_id.empty()) {
      name += ".";
      name += mpi_id;
    }
    if (!thread_id.empty()) {
      name += ".";
      name += thread_id;
    }
    if (!suffix.empty()) {
      name += ".";
      name += suffix;
    }
    return name;
  }

  // this function partitions a sequence into n sequential parts
  std::vector<sequence> sequential_equal_cut(const sequence& s, int n)
  {
    int psize = s.size() ;
    int div = psize/n ;
    int rem = psize%n ;
    vector<sequence> partition(n) ;

    sequence::const_iterator si = s.begin();
    for(int p=0;p<n;++p) {
      sequence& sp = partition[p];
      for(int k=0;k<div;++k,++si)
        sp += *si;
      if(rem > 0) {
        sp += *si;
        ++si;
        --rem;
      }
    }

    return partition;    
  }
  // this is the entitySet version.
  std::vector<entitySet> sequential_equal_cut(const entitySet& s, int n)
  {
    int psize = s.size() ;
    int div = psize/n ;
    int rem = psize%n ;
    vector<entitySet> partition(n) ;

    entitySet::const_iterator si = s.begin();
    for(int p=0;p<n;++p) {
      entitySet& sp = partition[p];
      for(int k=0;k<div;++k,++si)
        sp += *si;
      if(rem > 0) {
        sp += *si;
        ++si;
        --rem;
      }
    }

    return partition;    
  }

  sequence randomize_seq(const sequence& s)
  {
    vector<int> v;
    for(sequence::const_iterator si=s.begin();si!=s.end();++si)
      v.push_back(*si);
    random_shuffle(v.begin(),v.end());
    sequence o;
    for(vector<int>::const_iterator vi=v.begin();vi!=v.end();++vi)
      o += *vi;
    return o;
  }

  ////////////////////////////////////////////////////////////////////
  // implementation of the thread partitioners
  ////////////////////////////////////////////////////////////////////
  //
  // this function generates a reduction lock to be used with threads
  // blocks within the reduction and some of the chomping rules.  this
  // function has several steps involved:
  // 1) for all the blocks, we need to run them through the maps to
  //    calculate the target set they will compute, we will maintain
  //    such a map from blocks to their targets
  // 2) we will calculate all the pairs of block targets that are 
  //    actually in conflict (the target sets intersect or overlap),
  //    refering back to the blocks, if two block targets overlap,
  //    then these two blocks are in conflict and cannot be scheduled
  //    at the same time together.
  // 3) once we know the conflict block pairs, we can then generate a
  //    lock pattern to be used inside the ReductionLock type.
  //    the ReductionLock works by assigning one spin lock to each of
  //    the thread blocks.  each thread also maintains a list of conflict
  //    blocks they have, before a block is scheduled to a thread,
  //    all of the spin locks corresponding to the blocks on the list
  //    must be acquired, once successful, the block is then scheduled,
  //    upon completion, all the spin locks will be released.
  //    since each block will have the potential to acquire multi-locks,
  //    in order to avoid problems like livelock, we will sort all
  //    conflict blocks in each list according to their IDs, this way
  //    every one will acquire the spin locks in the same order, thus
  //    avoiding the potential thread synch problems.
  //
  // a comparison utility for block IDs
  inline bool block_id_cmp(const TBlockID& bid1, const TBlockID& bid2)
  { return bid1.first==bid2.first ? bid1.second<bid2.second
                                  : bid1.first<bid2.first; }

  // this creates an empty conflict lock, just for the purpose of debugging
  ReductionLock ThreadPartition_simple::generate_empty_conflict_lock() const
  {
    ReductionLock rl;
    rl.lock_pattern.resize(thread_num);
    for(int t=0;t<thread_num;++t) {
      rl.lock_pattern[t].resize(real_block_num[t]);
    }
    return rl;
  }

  // NOTE: this is just a temporary implementation.  Once the proper
  // infrastructure is in place (the Loci grid reader with partitioner
  // and maps on all levels, we will need to rewrite this part, hopefully
  // in a simpler and more efficient way)
  //
  // Note: we take in a vector of map_info for the reason that a thread
  // might execute a series of reduction rules together as a group.
  // this situation would happen for example, in the chomping graphs
  ReductionLock ThreadPartition_simple::generate_conflict_lock
    (const std::vector<std::vector<int> >& blocks_id,
     const std::vector<vmap_info>& map_info,
     fact_db& facts, sched_db& scheds) const
  {
    // first create a vector for all the block mapping results.
    map<TBlockID,entitySet> bmapr;
    // NOTE: the first dimension "t" ITSELF represents the thread id
    // the second dimension "b" however is just an index, the actual
    // value: blocks_id[t][b] represents the block id.
    for(size_t t=0;t<blocks_id.size();++t)
      for(size_t b=0;b<blocks_id[t].size();++b) {
        int tb = blocks_id[t][b];
        entitySet block_reach;

        for(size_t i=0;i<map_info.size();++i)
          block_reach += 
            vmap_target_exist(map_info[i],facts,blocks[t][tb],scheds);

        bmapr[make_pair(t,tb)] = block_reach;
      }
    // this is a naive n^2 algorithm that is slow, we will replace it
    // with a faster algorithm once we have all the parts in place.
    ReductionLock rl;
    rl.lock_pattern.resize(thread_num);
    for(int t=0;t<thread_num;++t) {
      rl.lock_pattern[t].resize(real_block_num[t]);
    }
#ifdef SHOW_LOCK_STAT
    int max_conflicts = 0;
    int max_non_local_conflicts = 0;
#endif
    for(map<TBlockID,entitySet>::const_iterator
        mi=bmapr.begin();mi!=bmapr.end();++mi) {
      vector<TBlockID> conflicts;
      TBlockID self = mi->first;
      const entitySet& self_targets = mi->second;
      for(map<TBlockID,entitySet>::const_iterator
          mi2=bmapr.begin();mi2!=bmapr.end();++mi2) {
        TBlockID other = mi2->first;
        const entitySet& other_targets = mi2->second;
        if( (self_targets & other_targets) != EMPTY)
          conflicts.push_back(other); 
      }
      sort(conflicts.begin(), conflicts.end());
#ifdef SHOW_LOCK_STAT
      if(max_conflicts < conflicts.size())
        max_conflicts = conflicts.size();
      int non_local_conflicts = 0;
#endif
      for(size_t c=0;c<conflicts.size();++c) {
#ifdef SHOW_LOCK_STAT
        if(conflicts[c].first != self.first)
          ++non_local_conflicts;
#endif
        rl.lock_pattern[self.first][self.second].
          push_back(&(ReductionLock::locks[conflicts[c].first]
                                          [conflicts[c].second]));
      }
#ifdef SHOW_LOCK_STAT
      if(max_non_local_conflicts < non_local_conflicts)
        max_non_local_conflicts = non_local_conflicts;
#endif
    }

#ifdef SHOW_LOCK_STAT
    std::cerr << "======== Thread Unit Lock Stat:" << endl;
    std::cerr << "== max conflicts: " << max_conflicts << endl;
    std::cerr << "== max non-local conflicts: " << max_non_local_conflicts
      << endl;
#endif

    return rl;
  }

  vector<sequence>
  ThreadPartition_simple::partition_by_domain(const sequence& s) const
  {
    vector<sequence> p(thread_num);
    // these intermediate conversions are needed because we cannot
    // perform set operations on sequence directly
    //
    // @PROBLEM:
    // we will probably need to take a look at the various operations
    // performed on entitySet and sequence in the code sometime.  many
    // of these would have to convert an entitySet to a sequence, or
    // convert a sequence to an entitySet, and then finally convert
    // the result back into one form.  this seems to be inefficient
    // and at least many of these operations can be made without such
    // conversions if perhaps there are some internal methods available
    // for a class of such operations.
    vector<entitySet> e(thread_num);
    entitySet ss = entitySet(s);
    for(int i=0;i<thread_num;++i) {
      e[i] = ss & ptn[i];
    }
    // It is possible in the current Loci setup for a passed in
    // sequence to have entities OUTSIDE of the thread domain
    // partition.  For example, the UNIT rule will have an execution
    // sequence including those clone region entities.  In this case,
    // we will need to check to make sure that all entities are
    // partitioned and returned, otherwise a rule may not execute
    // on all intended entities.  This may be a problem to fix at
    // the system level in the future.
    //
    // Right now, we will just cut the extra entities equally and
    // append them to the end of the original domain partition.
    // 
    entitySet extra = ss - whole_entities;
    vector<sequence> extra_seq_ptn =
      sequential_equal_cut(sequence(extra), thread_num);
    // finally convert to sequence again
    for(int i=0;i<thread_num;++i)
      p[i] = sequence(e[i]) + extra_seq_ptn[i];

    return p;
  }

  vector<int>
  ThreadPartition_simple::partition_by_blocks(int t, const sequence& s)
  const 
  {
    const std::vector<entitySet>& bs = blocks[t];
    vector<int> pb;
    entitySet se = entitySet(s);
    for(size_t i=0;i<bs.size();++i) {
      if( (se & bs[i]) != EMPTY)
        pb.push_back(i);
    }

    return pb;
  }

  sequence
    ThreadPartition_simple::refine_seq_by_block
    (int t, int b, const sequence& s) const
    {
      entitySet e = entitySet(s);
      return sequence(blocks[t][b] & e);
    }
  entitySet
    ThreadPartition_simple::refine_seq_by_block
    (int t, int b, const entitySet& s) const
    {
      return blocks[t][b] & s;
    }

  namespace {
    // small comparison utility used in thread entities partition
    struct SortNodePtnSize {
      bool operator()(const pair<int,int>& p1,
                      const pair<int,int>& p2) const
      { return p1.second < p2.second; }
    };
  }

  // this is a helper function to create thread partition and blocks.
  // it is however the central function in the process.  the idea is
  // that giving the necessary inputs (several stores and maps from
  // which we can get related information on cells, faces, nodes, etc.),
  // this function will create "p" partitions that tries to maximize
  // locality according to the mesh topology.
  //
  // first of all, the cells are cut sequentially into "p" parts, this
  // assumes that cells have already been organized in a spatially 
  // related way (such as by using a space-filling curve typed arrangement).
  //
  // then each cell partition will also own the faces and nodes that
  // belong to the cells in such a partition.  the only source of problem
  // comes from the fact that some of the faces and nodes will be shared
  // among multiple cell partitions.  we will resolve this issue by
  // assigning these "partition boundary" faces and nodes to cell partitions
  // that tries to balance the number of total entities in all of the
  // partitions.
  //
  // Note: this function is used to initially create "p" partitions of
  // the original MPI entity set, "p" equals the number of threads
  // multiply the number of blocks per thread.  the reason for doing so
  // is because we want to partition the thread domain as well as the
  // thread blocks both in such a way to try to minimize inter-domain
  // and inter-block communications.  after we created the "p" partitions,
  // we will then group "b" such partitions together to form a thread
  // partition, where "b" is the number of blocks per thread.  notice 
  // that we will sequentially group each "b" partitions for reasons 
  // that such grouping will also preserve the cells partition locality
  // (and hence the faces and nodes as well) so that the grouped thread
  // domain locality is also preserved.
  //
  // Also note: this function only partitions the cells, faces, and nodes,
  // it does not try to partition any other entities, such as the 
  // boundary faces, etc.  all other entities are partitioned afterwards
  // in a separate step (mainly because the topological relationship
  // between those entities is less clear for a locality type partition).
  // we will mainly just use a simple round robin type partition for those
  // entities.
  //
  void ThreadPartition_simple::partition_topology
    (constraint& geom_cells, const Map& cl, const Map& cr,
     const multiMap& face2node,
     // limit is the total entities that this MPI process owns, i.e.,
     // we cannot include entities within any of the partitions 
     // that are beyond this limit.
     // p is the number of partitions, ptn is the output of partition
     const entitySet& limit, int p, vector<entitySet>& ptn)
  {
    entitySet cells = *geom_cells;
    cells &= limit;
    vector<entitySet> cell_ptn(p);

    // compute per thread size 
    int psize = cells.size();
    int div = psize/p;
    int rem = psize%p;

    // during the cell partition, we will build a cell -> partition map,
    // i.e., a structure that records which cell belongs to which partition
    map<Entity,int> cell_2_ptn;
    entitySet::const_iterator ei=cells.begin();
    for(int i=0;i<p;++i) {
      for(int k=0;k<div;++k,++ei) {
        cell_ptn[i] += *ei;
        cell_2_ptn[*ei] = i;
      }
      if(i<rem) {
        cell_ptn[i] += *ei;
        cell_2_ptn[*ei] = i;
        ++ei;
      }
    }

    // now partition the faces, first get all the local faces
    entitySet faces = cl.domain() & cr.domain();
    faces &= limit;
    vector<entitySet> face_ptn(p);
    map<Entity,int> face_2_ptn;

    // we will first loop over all the faces and check to see if
    // the cells on both sides of a face belong to the same thread?
    // if it is, then clearly that face belongs to the same thread,
    // if not, then the face sits on the partition boundary, we will
    // record them at this step.
    entitySet partition_boundary_faces;
    for(ei=faces.begin();ei!=faces.end();++ei) {
      map<Entity,int>::const_iterator mi;
      Entity lc = cl[*ei];
      Entity rc = cr[*ei];
      int plc=-1, prc=-1;
      mi = cell_2_ptn.find(lc);
      if(mi != cell_2_ptn.end())
        plc = mi->second;
      mi = cell_2_ptn.find(rc);
      if(mi != cell_2_ptn.end())
        prc = mi->second;
      // if either plc and prc equals -1, that means
      // the corresponding lc or rc cells are not in this local
      // partition, or not a cell (e.g., a boundary)
      if(plc==-1 && prc==-1) {
        // something is wrong, they cannot be both invalid locally
        cerr << "thread face partition failed!!!" << endl;
        Loci::Abort();
      }

      if(plc == -1) {
        face_ptn[prc] += *ei;
        face_2_ptn[*ei] = prc;
      } else if(prc == -1) {
        face_ptn[plc] += *ei;
        face_2_ptn[*ei] = plc;
      } else if(plc == prc) {
        face_ptn[plc] += *ei;
        face_2_ptn[*ei] = plc;
      } else
        partition_boundary_faces += *ei;
    }

    // the way we will assign faces on the partition boundary to threads
    // is to try to balance the size of the face partition.  for example,
    // if face "f" sits on the partition boundary of partition 1 (p1) and
    // partition 2 (p2), then we will assign f to p1 if the current face 
    // partition size of t1 is smaller than that of p2, and vice versa.
    // if the face partition of p1 and p2 are of equal size, then f
    // is assigned to either partition.  note: the face partition size of
    // all partitions are checked dynamically as we move through the list
    // of partition boundary faces.
    for(ei=partition_boundary_faces.begin();
        ei!=partition_boundary_faces.end();++ei) {
      // note: we do not have to use the "find()" method on the map
      // "cell_2_ptn" here to check for search errors because we've
      // already done so in the previous step and we can be sure
      // that all the faces in this set will succeed in the search
      int plc = cell_2_ptn[cl[*ei]]; // owning partition of the cell on left
      int prc = cell_2_ptn[cr[*ei]]; // owning partition of the cell on right
      int sizel = face_ptn[plc].size();
      int sizer = face_ptn[prc].size();
      if(sizel < sizer) {
        face_ptn[plc] += *ei;
        face_2_ptn[*ei] = plc;
      } else {
        face_ptn[prc] += *ei;
        face_2_ptn[*ei] = prc;
      }
    }

    // next we will be partitioning nodes.  we will just use the same
    // exact method as we did for faces.
    //
    // first we build a node2face reverse map from the face2node map
    //
    entitySet nodes;
    map<Entity,entitySet> node2face;
    for(ei=faces.begin();ei!=faces.end();++ei) {
      int nn = face2node.num_elems(*ei);
      for(int i=0;i<nn;++i) {
        Entity nd = face2node[*ei][i];
        nodes += nd;
        node2face[nd] += *ei;
      }
    }
    nodes &= limit;

    vector<entitySet> node_ptn(p);
    entitySet partition_boundary_nodes;

    for(ei=nodes.begin();ei!=nodes.end();++ei) {
      const entitySet& fs = node2face[*ei];
      // check to see if all the faces belong to the same partition
      bool assign = true;
      int proc = -1;
      for(entitySet::const_iterator fi=fs.begin();fi!=fs.end();++fi) {
        map<Entity,int>::const_iterator mi = face_2_ptn.find(*fi);
        if(mi != face_2_ptn.end()) {
          if(proc == -1)
            proc = mi->second;
          else {
            if(proc != mi->second) {
              assign = false;
              break;
            }
          }
        } else {
          // if we couldn't find the face record in the face_2_ptn
          // map, then it means that this face is not in the local
          // face set (i.e., this is a remote face owned by other
          // MPI processes), in this case, we'll just ignore it and
          // do nothing
        }
      }
      if(proc == -1) {
        // this could not happen since at least one of the faces
        // that node belongs to should be in the local entity set
        cerr << "thread node partition failed!!!" << endl;
        Loci::Abort();
      }
      if(assign)
        node_ptn[proc] += *ei;
      else
        partition_boundary_nodes += *ei;
    }
    
    // then work on the partition boundary nodes to try to balance
    // the node partition size among all the threads.
    for(ei=partition_boundary_nodes.begin();
        ei!=partition_boundary_nodes.end();++ei) {
      const entitySet& fs = node2face[*ei];
      vector<pair<int,int> > node_ptn_size;
      for(entitySet::const_iterator fi=fs.begin();fi!=fs.end();++fi) {
        map<Entity,int>::const_iterator mi = face_2_ptn.find(*fi);
        // only work on faces that are owned by the local MPI process
        if(mi != face_2_ptn.end()) {
          int proc = mi->second; 
          node_ptn_size.push_back(make_pair(proc,node_ptn[proc].size()));
        }
      }
      if(node_ptn_size.empty()) {
        cerr << "thread node partition failed!!!" << endl;
        Loci::Abort();
      }
      sort(node_ptn_size.begin(), node_ptn_size.end(), SortNodePtnSize());
      // the first entry in the vector is the thread to give the node to
      node_ptn[node_ptn_size[0].first] += *ei;
    }

    // finally, we merge all the partitions
    ptn.resize(p);
    for(int i=0;i<p;++i)
      ptn[i] = cell_ptn[i]+face_ptn[i]+node_ptn[i];
  }

  void ThreadPartition_simple::create_thread_domain(fact_db& facts)
  {
    // Note: this function is run on the main control thread
    // initially, there are no work threads created yet at this
    // moment
    //
    // This is a test version that just uses customized version of
    // partition of the entities among the threads.  When the system
    // is ready, we will interface this part with the Loci partitioner
    // to gain access of the capability of that infrastructure.
    //
    // first get all of the entities assigned to the local MPI process
    // if fact database is distributed, then we will trim the entities
    // down to those that are owned by the local process
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info();
      whole_entities = d->my_entities;
    } else {
      // note, this probably needs to be fixed for big mode
      whole_entities = facts.get_init_ptn(0)[Loci::MPI_rank] ;
    }

    // we will partition the entities among all the work threads.
    // this partition works by first partitioning the cells using a
    // simple equal sequential cut (into p parts with each thread
    // taking one part).  then we will distribute all the faces and nodes
    // according to the cell partition.  for some faces and nodes,
    // they will be shared among multiple threads (since they are on
    // the boundaries of partition), we will then try to determine
    // a unique ownership of these faces and nodes by trying to
    // balancing the loads (see below for detailed comments).
    //
    // first get the revelant variables from the fact database
    constraint geom_cells = facts.get_variable("geom_cells");
    Map cl(facts.get_variable("cl"));
    Map cr(facts.get_variable("cr"));
    multiMap face2node(facts.get_variable("face2node"));
    store<string> boundary_tags(facts.get_variable("boundary_tags"));

    int p = thread_num * initial_block_num; 
    vector<entitySet> init_ptn;
    partition_topology
      (geom_cells,cl,cr,face2node,whole_entities,p,init_ptn);

    // what we do next is to make sure that no init_ptns are empty!
    // this is possible if the size of the initial entities to start with
    // is too small, or the number of thread blocks is way to large.
    // if this happens, then we will remove those empty partitions
    // and reassign the partitions to threads.
    int real_blocks = 0;
    vector<int> valid_ptn; // the index of non-empty init_ptn
    for(size_t i=0;i<init_ptn.size();++i)
      if(init_ptn[i].size() != 0) {
        real_blocks++; 
        valid_ptn.push_back(i);
      }

    if(real_blocks < thread_num) {
      // we are really in trouble here, the blocks partitioning
      // cannot even sustain one block per thread, in this case,
      // we will simply just give up partitioning at all, we will
      // use a simple straightforward partitioning method to
      // cut the whole_entities into "t" equal parts, where "t"
      // is the number of threads, also we will use one block
      // per thread.  note: some thread domain may even be empty
      // in this case, but this is okay presumably since we have
      // such a small initial set to start with.
      ptn = sequential_equal_cut(whole_entities, thread_num);
      for(int i=0;i<thread_num;++i) {
        real_block_num[i] = 1;
        vector<entitySet> b; b.push_back(ptn[i]);
        blocks[i] = b; 
      }
      total_blocks = thread_num;

      ReductionLock::locks = new pthread_spinlock_t*[blocks.size()];
      for(size_t i=0;i<blocks.size();++i) {
        ReductionLock::locks[i] = new pthread_spinlock_t[blocks[i].size()];
        for(size_t k=0;k<blocks[i].size();++k) {
          if(pthread_spin_init(&(ReductionLock::locks[i][k]), 
                PTHREAD_PROCESS_PRIVATE)) {
            // need better error handling later...
            Loci::Abort();
          }
        }
      }
      // give some warning at least
      Loci::debugout << "Thread domain creation given very small input..."
        << endl;
      Loci::debugout
        << "    creating really small partitions and single block per thread"
        << endl;

      cout << "Thread domain creation given very small input..." << endl;
      cout << "    creating really small partitions and single block per thread" << endl;
      return;
    }

    if(real_blocks < p) {
      // in this case, we do not meet the initial specification
      // (i.e, the given number of blocks per thread cannot be provided
      // due to small input entities size, or extremely large block number).
      // what we will do here is to adjust the number of blocks per
      // thread so that each thread will get approximately the same
      // partition size.

      // but we will issue some warning to let the user know about this
      int nb = real_blocks/thread_num;
      int nm = real_blocks%thread_num;
      for(int i=0;i<thread_num;++i) {
        if(i<nm)
          real_block_num[i] = nb+1;
        else
          real_block_num[i] = nb;
      }
      Loci::debugout << "Thread partition readjusted block numbers..."
        << endl;
      cout << "Thread partition readjusted block numbers..." << endl;
    } else {
      // this is the normal case where we can supply the initially
      // specified number of blocks to each thread with sufficient
      // entities inside
      for(int i=0;i<thread_num;++i)
        real_block_num[i] = initial_block_num;
    }
    total_blocks = real_blocks;

    // before continue, we will need to partition the remaining
    // entities in a round-robin fashion.
    entitySet ptned;
    for(size_t i=0;i<init_ptn.size();++i)
      ptned += init_ptn[i];
    
    entitySet remaining = whole_entities - ptned;

    vector<entitySet> remain_ptn(real_blocks);
    entitySet::const_iterator ei = remaining.begin();
    int proc = 0;
    while(ei != remaining.end()) {
      remain_ptn[proc] += *ei; ++ei;
      ++proc;
      if(proc==real_blocks)
        proc = 0;
    }

    // we are ready to create the final thread domain and blocks
    ptn.resize(thread_num);
    int idx = 0;
    for(int i=0;i<thread_num;++i) {
      int bn = real_block_num[i];
      vector<entitySet> tbs(bn);
      ptn[i] = EMPTY;
      for(int k=0;k<bn;++k,++idx) {
        tbs[k] = init_ptn[valid_ptn[idx]]+remain_ptn[idx];
        ptn[i] += tbs[k];
      }
      blocks[i] = tbs;
    }

    // once we know all the thread blocks, we will then initialize
    // the thread block reduction locks
    ReductionLock::locks = new pthread_spinlock_t*[blocks.size()];
    for(size_t i=0;i<blocks.size();++i) {
      ReductionLock::locks[i] = new pthread_spinlock_t[blocks[i].size()];
      for(size_t k=0;k<blocks[i].size();++k) {
        if(pthread_spin_init(&(ReductionLock::locks[i][k]), 
              PTHREAD_PROCESS_PRIVATE)) {
          // need better error handling later...
          Loci::Abort();
        }
      }
    }

  }

  ThreadPartition_simple::~ThreadPartition_simple()
  {
    // we will destroy the reduction locks initialized
    for(size_t i=0;i<blocks.size();++i) {
      for(size_t k=0;k<blocks[i].size();++k)
        pthread_spin_destroy(&(ReductionLock::locks[i][k]));
      delete[] ReductionLock::locks[i];
    }
    delete[] ReductionLock::locks;
  }

  const vector<entitySet>&
   ThreadPartition_simple::get_thread_block_list(int t) const
  { 
    if (t<0 || t>=thread_num)
      throw ThreadException("ThreadPartition_simple: invalid thread ID!");
    return blocks[t]; 
  }

  ThreadControl_pthread::~ThreadControl_pthread()
  {
    for(int i=0;i<tnum;++i) {
      sem_destroy(&worker_block[i]);
    }
    pthread_attr_destroy(&pattr);
    sem_destroy(&main_block);
    for(size_t i=0;i<finish_notification.size();++i) {
      pthread_spin_destroy(&(finish_notification[i].done_lock));
    }
    pthread_barrier_destroy(&barrier);
    delete tpn;
  }

  bool ReductionLock::try_lock(int tid, int bid)
  {
    const vector<pthread_spinlock_t*>& sub_locks = lock_pattern[tid][bid];
    for(size_t i=0;i<sub_locks.size();++i) {
      if(pthread_spin_trylock(sub_locks[i]) != 0) {
        // roll back all locks held and return
        for(size_t k=0;k<i;++k)
          pthread_spin_unlock(sub_locks[k]);
        return false;
      }
    }
    return true;
  }

  void ReductionLock::unlock(int tid, int bid)
  {
    const vector<pthread_spinlock_t*>& sub_locks = lock_pattern[tid][bid];
    for(size_t i=0;i<sub_locks.size();++i)
      pthread_spin_unlock(sub_locks[i]);
  }

  int ReductionLock::try_lock_entries(int tid, int bid)
  {
    return lock_pattern[tid][bid].size();
  }

  // these two functions build the termination tree
  
  void ThreadControl_pthread::
  build_term_tree(int pid, int depth, size_t& myid)
  {
    if(myid >= finish_notification.size())
      return;

    size_t me = myid;
    finish_notification[me].done = false;
    if(pthread_spin_init(&(finish_notification[me].done_lock), 
          PTHREAD_PROCESS_PRIVATE)) {
      // need better error handling later...
      Loci::Abort();
    }
    finish_notification[me].p = pid;
    if(pid >= 0) {
      finish_notification[pid].c.push_back(me);
    }
    ++myid;

    if(depth>0) {
      build_term_tree(me, depth-1, myid);
      build_term_tree(me, depth-1, myid);
    }
  }

  // a function that computes the nearest power of 2 for a given integer
  int pow2_roundup(int n)
  {
    int start = 1;
    while(start < n)
      start <<= 1;
    return start;
  }

  // build a list of partners for the parallel starting process
  void build_partner(int me, int limit, int stride,
                     vector<vector<int> >& partner)
  {
    stride >>= 1;
    if(stride < 1)
      return;
    int p = me + stride;
    if(p < limit)
      partner[me].push_back(p);
    build_partner(me, limit, stride, partner);
    if(p < limit)
      build_partner(p, limit, stride, partner);
  }
  
  ThreadControl_pthread::ThreadControl_pthread
    (int n, fact_db& facts, sched_db& scheds)
    :active(false),stop(false),work(0),factsP(0),schedsP(0),tnum(n)
  {
    thread_args.resize(n);
    worker_block.resize(n);
    if(sem_init(&main_block, 0, 0)) {
      // need to better handle errors...
      Loci::Abort();
    }

    // build the parallel start partner
    int stride = pow2_roundup(tnum);
    start_partners.resize(tnum);
    build_partner(0, tnum, stride, start_partners);

    // build the termination notification tree
    int size = tnum;
    finish_notification.resize(size);
    int depth = 0;
    while(size > 1) {
      ++depth;
      size /= 2;
    }
    size_t init_id = 0;
    build_term_tree(-1, depth, init_id);

    // build the reduction partner (just the reverse of the start_partners)
    reduction_partners.resize(tnum);
    for(int i=0;i<tnum;++i) {
      vector<int>::const_reverse_iterator b,e;
      b = start_partners[i].rbegin();
      e = start_partners[i].rend();
      reduction_partners[i] = vector<int>(b,e);
    }

    // system contention scope
    pthread_attr_init(&pattr);
    pthread_attr_setscope(&pattr, PTHREAD_SCOPE_SYSTEM);
    for(int i=0;i<n;++i) {
      if(sem_init(&worker_block[i], 0, 0)) {
        // need to better handle errors...
        Loci::Abort();
      }
    }

    stop = false;

    if(pthread_barrier_init(&barrier, 0, tnum)) {
      cerr << "ThreadControl barrier init failed!!!" << endl;
      Loci::Abort();
    }    

    // then we need to create a thread entitySet manager
    tpn = new ThreadPartition_simple(tnum);
    tpn->create_thread_domain(facts);
  }

  vector<sequence>
  ThreadControl_pthread::partition_seq(const sequence& s) const
  {
    return tpn->partition_by_domain(s);
  }

  vector<int>
  ThreadControl_pthread::block_thread_seq(int t, const sequence& s) const
  { return tpn->partition_by_blocks(t,s); }

  sequence
    ThreadControl_pthread::refine_seq_by_block
    (int t, int b, const sequence& s) const
    { return tpn->refine_seq_by_block(t,b,s); }
  entitySet
    ThreadControl_pthread::refine_seq_by_block
    (int t, int b, const entitySet& s) const
    { return tpn->refine_seq_by_block(t,b,s); }

  ReductionLock
    ThreadControl_pthread::create_empty_conflict_lock() const
    { 
      return tpn->generate_empty_conflict_lock();
    }
  ReductionLock
    ThreadControl_pthread::create_conflict_lock
    (const vector<vector<int> >& blocks, const vmap_info& map_info,
     fact_db& facts, sched_db& scheds) const
    { 
      vector<vmap_info> vm; vm.push_back(map_info);
      return tpn->generate_conflict_lock(blocks,vm,facts,scheds);
    }
  ReductionLock
    ThreadControl_pthread::create_conflict_lock
    (const vector<vector<int> >& blocks,
     const vector<vmap_info>& map_info, fact_db& facts, sched_db& scheds)
    const
    { return tpn->generate_conflict_lock(blocks,map_info,facts,scheds);}

  void ThreadControl_pthread::create_threads()
  {
    // first record the main control thread id
    main_pid = pthread_self();
    id_map[main_pid] = "tmain";
    // creating threads
    for(int i=0;i<tnum;++i) {
      thread_args[i].tid = i;
      thread_args[i].self = this;
      int rc = pthread_create(&(thread_args[i].ptid), &pattr, thread_fun,
                              (void*)(&thread_args[i]));
      if(rc) {
        // need to handle errors better later...
        Loci::Abort();
      }
      string id = "t";
      stringstream ss; ss << i;
      id += ss.str();
      id_map[thread_args[i].ptid] = id;
    }
    // record the leading work thread id
    lead_work_pid = thread_args[0].ptid;
  }

  string ThreadControl_pthread::get_local_id() const
  {
    pthread_t self = pthread_self();
    map<pthread_t,string>::const_iterator mi = id_map.find(self);
    if (mi != id_map.end())
      return mi->second;
    else
      return "tUNKNOWN";
  }

  void ThreadControl_pthread::global_barrier()
  {
    int rc = pthread_barrier_wait(&barrier);
    if(rc!=0 && rc!=PTHREAD_BARRIER_SERIAL_THREAD) {
      cerr << "ThreadControl barrier wait failed!!!" << endl;
      Loci::Abort();
    }
  }

  // determine if the calling thread is leading thread or not.
  // the leading thread in the pthread implementation is defined
  // to be either the main control thread or the leading work 
  // thread.  This is valid since in our design, the main control
  // thread will never work together with any of the work threads.
  bool ThreadControl_pthread::is_leading_thread() const
  {
    pthread_t sid = pthread_self();
    return (pthread_equal(sid,main_pid) || 
            pthread_equal(sid,lead_work_pid));
  }

  const vector<int>& ThreadControl_pthread::
  get_reduction_partners(int id) const
  {
    if(id < 0 || id >= tnum)
      throw ThreadException("get_reduction_partners: thread id not valid");
    return reduction_partners[id];
  }

  // this is the thread function to execute on each thread upon creation.
  // basically it just waits until there is work to do.
  void* ThreadControl_pthread::thread_fun(void* arg)
  {
    // retrieve the argument
    TArg* targ = (TArg*)arg;
    int tid = targ->tid;
    ThreadControl_pthread* self = targ->self;

    // // we want to set the CPU affinity of threads so that
    // // each thread is attached to a specific CPU and wont migrate (hopefully)
    // cpu_set_t cpuset;
    
    // CPU_ZERO(&cpuset);
    // CPU_SET(Loci::MPI_rank*self->tnum + tid, &cpuset);
    
    // if( (pthread_setaffinity_np(pthread_self(),
    //                             sizeof(cpu_set_t), &cpuset)) != 0)
    //   cerr << "WARNING: pthread CPU affinity setting failed!" << endl;

    // this holds the pointer to the next work module
    ThreadedEMP w = 0;

    while(true) {
      // first block on the semaphore waiting for the
      // go signal from the main thread
      sem_wait(&(self->worker_block[tid]));

      // parallel start...
      for(size_t i=0;i<(self->start_partners)[tid].size();++i) {
        int p = (self->start_partners)[tid][i];
        sem_post(&((self->worker_block)[p]));
      }

      if(self->stop)
        break;

      // do the module work
      w = (*(self->work))[tid];
      fact_db* f = self->factsP;
      sched_db* s = self->schedsP;
      if(w != 0 && f && s) {
        w->execute(*f, *s);
      }

      // propagate the done flag so that eventually the main thread
      // is notified that all are finished.

      // first of all, we will check are all children in the tree done yet?
      vector<Terminate>& fn = self->finish_notification;
      for(size_t i=0;i<fn[tid].c.size();++i) {
        int cid = fn[tid].c[i];
        bool done = false;
        while(!done) {
          pthread_spin_lock(&(fn[cid].done_lock));
          done = fn[cid].done;
          pthread_spin_unlock(&(fn[cid].done_lock));
        }
        // we will then reset the done flag in the child
        // we don't need to acquire the lock at this time since
        // this child should have finished by now.
        fn[cid].done = false;
      }
      // then we need to write the done flag
      // the root thread needn't do this since it just needs to
      // notify the main thread
      if(fn[tid].p < 0) {
        sem_post(&self->main_block);
      } else {
        pthread_spin_lock(&(fn[tid].done_lock));
        fn[tid].done = true;
        pthread_spin_unlock(&(fn[tid].done_lock));
      }
    }
    return 0;
  }

  void ThreadControl_pthread::shutdown()
  {
    // only the main thread is able to call this function
    // also, if the main thread is active, then all the
    // work threads must be blocked, so it is safe to update
    // any data at this moment.
    stop = true;
    vector<ThreadedEMP> w(tnum, ThreadedEMP(0));
    restart(w);
    // joining all the threads
    for(int i=0;i<tnum;++i) {
      if(pthread_join(thread_args[i].ptid,NULL)) {
        // better error handling later...
        Loci::Abort();
      }
    }
    active = false;
  }
  
  void ThreadControl_pthread::
  restart(std::vector<ThreadedEMP>& w)
  {
    if( (int)w.size() != tnum) {
      throw ThreadException("restart work size wrong!");
    }
    work = &w;
    active = true;
    sem_post(&worker_block[0]);
  }
  
  void ThreadControl_pthread::wait_threads()
  {
    if(!active)
      return;

    sem_wait(&main_block);
    active = false;
  }
  
  void ThreadControl_pthread::
  sequential_restart(std::vector<ThreadedEMP>& w)
  {
    if( (int)w.size() != tnum) {
      throw ThreadException("sequential_restart work size wrong!");
    }
    std::vector<ThreadedEMP> nw(tnum);
    for(int t=0;t<tnum;++t) {
      for(int i=0;i<tnum;++i)
        nw[i] = ThreadedEMP(0);
      nw[t] = w[t];
      restart(nw);
      wait_threads();
    }
  }

  //pthread_mutex_t copy_mutex = PTHREAD_MUTEX_INITIALIZER; 
  
} // end of namespace Loci

#else
#include "distribute.h"
namespace Loci {
  // this implements the leading execution predicate
  bool is_leading_execution()
  {
    return MPI_rank==0 ;
  }

  // these two functions are intended to provide a global atomic
  // region of execution (the ideal use of them is to use the 
  // Loci preprocessor to hide these calls from the users)
  void global_atomic_region_begin()
  {
  }

  void global_atomic_region_end()
  {
  }
}
#endif
