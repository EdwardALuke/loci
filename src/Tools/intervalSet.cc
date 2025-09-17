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
//##########################################################################

#include <Tools/intervalSet.h>
namespace Loci {
  // provide generic int pair<int,int> output
  std::ostream & operator <<(std::ostream &s, const std::pair<int,int> &p) {
    s << '(' << p.first <<',' << p.second << ')' ;
    return s ;
  }
  std::istream &operator >>(std::istream &s, std::pair<int,int> &p) {
    while(s.peek() == ' ')
      s.get() ;
    if(s.peek() == '(')
      s.get() ;
    s >> p.first  ;
    while(s.peek() == ' ')
      s.get() ;
    if(s.peek() == ',')
      s.get() ;
    s >> p.second ;
    while(s.peek() == ' ')
      s.get() ;
    if(s.peek() == ')')
      s.get() ;
    return s ;
  }
  
  // This code apparently has problems with the latest intel compiler
  // It seems that these initializations were occuring before the 
  // genIntervalSet template static members were initialized.
  // This change fixes it, though it is ugly.
  const genIntervalSet<int_type> EMPTY = genIntervalSet<int_type>() ;
  const int_type UNIVERSE_MAX = std::numeric_limits<int_type>::max() - 1 ;
  const int_type UNIVERSE_MIN = std::numeric_limits<int_type>::min() + 1 ;

  const genIntervalSet<gEntity> GEMPTY = genIntervalSet<gEntity>() ;
  const gEntity GUNIVERSE_MAX = std::numeric_limits<gEntity>::max() - 1 ;
  const gEntity GUNIVERSE_MIN = std::numeric_limits<gEntity>::min() + 1 ;
  

}
