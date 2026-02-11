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
#include <Loci.h>
#include <parSampleSort.h>

#include <vector>
#include <string>
namespace VOG {

  struct BC_descriptor {
    std::string name ;
    int id ;
    bool Trans ;
  } ;

  extern std::vector<BC_descriptor> readTags(std::string filename) ;

  using std::vector ;

  // Optimize indicies of mesh to increase locality
  extern void optimizeMesh(store<vector3d<double> > &pos,
                           Map &cl, Map &cr, multiMap &face2node) ;
  // Establish geometrically consistent face orientation
  extern void orientFaces(store<vector3d<double> > &pos,
                          Map &cl, Map &cr, multiMap &face2node) ;
  extern void colorMatrix(store<vector3d<double> > &pos,
                          Map &cl, Map &cr, multiMap &face2node) ;
  extern std::vector<int> simplePartitionVec(int mn, int mx, int p) ;

  /// Return true if `value` ends with `suffix`.
  extern bool hasSuffix(const std::string &value, const std::string &suffix) ;
  /// Return `value` with trailing `suffix` removed, if present, otherwise
  /// return `value` unchanged.
  extern std::string stripSuffix(const std::string &value,
                                 const std::string &suffix) ;
  /// Return `value` with `suffix` appended, unless it already has it, in
  /// which case return `value` unchanged.
  extern std::string ensureSuffix(const std::string &value,
                                  const std::string &suffix) ;

}
