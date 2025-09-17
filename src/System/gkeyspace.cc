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
#include <gkeyspace.h>
#include <gkey_manager.h>
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

namespace Loci {

  //given send_split, return recv_split;
  //or given recv_split, return send_split
  //allow overlap between processes
  extern std::vector<gEntitySet> transposePtn(const std::vector<gEntitySet> &ptn, MPI_Comm comm);
  //static variable
  gKeySpace::gKeySpaceManager* gKeySpace::ksm=0;
  std::map<std::string, gKeySpaceP>  gKeySpace::gKeySpaceManager::space_map;

  //mutator for send_ptn and recv_ptn
  //ptn is send_ptn, since recv_ptn is always transposePtn of send_ptn
  //recv_ptn is modified whenever send_ptn is modified 
  void gKeySpace::set_send_recv_ptn(const std::vector<gEntitySet>& ptn){
    send_ptn = ptn;
    if( np>1 ) recv_ptn = transposePtn(ptn, comm);
    else recv_ptn = ptn;
  }
  

  gKeySpaceP gKeySpace::get_space(const std::string& spacename, const std::string& casename){
    return ksm->get_space(spacename, casename);
  }

  std::vector<gKeySpaceP> gKeySpace::get_all_spaces(){
    return ksm->get_all_spaces();
  }

  
  void gKeySpace::register_space(const std::string& spacename, const std::string& casename){
    ksm->add_space(spacename, casename, gKeySpaceP(this));
  }
  
} // end of namespace Loci



