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
#ifndef LOCI_GLOBS_H
#define LOCI_GLOBS_H

#include <fstream>
#include <string>
#include <vector>

namespace Loci {
  /////////////////////////////
  // flags to turn on/off the visualization feature
  extern bool verbose ;
  extern bool show_graphs ;
  extern bool show_decoration ;
  // flag to enable/disable the dynamic memory management
  extern bool use_dynamic_memory ;
  // flag to enable/disable output of dynamic memory
  // and multilevel graph decoration information
  extern bool show_dmm_verbose ;
  // flag to enable/disable chomping
  extern bool use_chomp ;
  // flag to enable outputing schedule to file
  extern bool schedule_output ;
  //flag for duplicating computations to save communication
  extern bool duplicate_work;
  extern bool multilevel_duplication;
  extern bool reduction_duplication;
  extern bool pointwise_duplication;
  extern bool extended_duplication;
  extern bool collect_timings;
  extern std::ofstream timeout;
  extern double time_duration_to_collect_data;
  extern bool use_duplicate_model;
  extern bool profile_memory_usage ;
  extern char * model_file;
  /////////////////////////////
  extern bool random_partition;
  
  extern bool measure_rule_timings;

  extern int printLevel ;
  inline std::ostream &printIndent(std::ostream &s) {
    for(int i=0;i<printLevel;++i)
      s << "  " ;
    return s ;
  }
      
}

#endif
