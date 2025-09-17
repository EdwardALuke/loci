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
#ifndef VISITORABS_H
#define VISITORABS_H

// the abstract base class of visitor
// used as interface for other files (comp_dag.cc etc)
namespace Loci {
  
  class loop_compiler ;
  class dag_compiler ;
  class conditional_compiler ;
  class chomp_compiler ;
  class impl_recurse_compiler ;
  class recurse_compiler ;
  
  class visitor {
  public:
    virtual ~visitor() {}
    // the following three are mandatory methods
    virtual void visit(loop_compiler&) = 0 ;
    virtual void visit(dag_compiler&) = 0 ;
    virtual void visit(conditional_compiler&) = 0 ;
    // the recursive compiler may not needed in the
    // visit, therefore we supply two empty defaults
    virtual void visit(impl_recurse_compiler&) {}
    virtual void visit(recurse_compiler&) {}
  } ;
}

#endif
