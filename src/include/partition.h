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
#ifndef PARTITION_H
#define PARTITION_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>
#include <data_traits.h>
#include <gstore.h>
#include <gmultimap.h>
#include <gmap.h>
#include <iterator>
#include <string>
#include <vector>
#include <utility>
#include <distribute_long.h>
namespace Loci {
  class fact_db;
  // //assume in1, in2 are sorted and have the same domain
  // //out will be allocalted on in2's image
  // template<class T>  void localJoin(const gStore<T> &in1, const gMultiMap &in2, gStore<T> &out) {
  //   // Find pairs where first entry are the same and create joined protomap
  //   out.clear() ;
  //   typename gStore<T>::const_iterator itr1 = in1.begin();
  //   gMultiMap::const_iterator itr2 = in2.begin(); 
  //   size_t j = 0 ;
  //   for(;itr1 != in1.end();++itr1) {
  //     while(itr2 != in2.end() && itr1->first < itr2->first )
  //       itr2++ ;
  //     gMultiMap::const_iterator itrs=itr2 ;
  //     while(itrs != in2.end() && itrs->first == itr1.first) {
  //       out.insert(itr1->second,itrs->second) ;
  //       itrs++ ;
  //     }
  //   }
  // }

  struct parStrategy{
    gKeySpaceP space1; //first space to partition
    gKeySpaceP space2; //second space to partition
    gKeySpaceP space3; //thrid space to partition
    std::string strategy; //the strategy to partition first space
    std::string map1; //the map used to create selfmap in space1 if metis partition is used
    std::string map2; //the map used in affinity partition when space2 is partitioned
    std::string map3; //the map used in affinity partition when space3 is partitioned
    bool from_space1;//is space3  using space1's partition 
  };
  
  void primary_partition_metis(fact_db& facts, std::vector<int>& e2p,  parStrategy& s);
  void primary_partition_orb(fact_db& facts,  std::vector<int>& e2p,  parStrategy& s );
  void affinity_partition(fact_db &facts, const vector<int>& procmap,  parStrategy& s );
  //given send_split, return recv_split;
  //or given recv_split, return send_split
  //allow overlap between processes
  std::vector<gEntitySet> transposePtn(const std::vector<gEntitySet> &ptn, MPI_Comm comm LOCI_DEFAULT_COMM);
}


#endif
