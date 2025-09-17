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
#ifndef IDENTITY_H
#define IDENTITY_H

#include <typeinfo>

#include <Tools/debug.h>
#include <Tools/lmutex.h>

namespace Loci {
typedef unsigned long id_type ;

class IdentityIMPL {
  id_type id ;
  IdentityIMPL & operator=(const IdentityIMPL &i)
  { warn(true) ; return *this; }
  public:
    IdentityIMPL() { id = reinterpret_cast<id_type>(this) ; }
    IdentityIMPL(const IdentityIMPL &i) { id = reinterpret_cast<id_type>(this);} 
    ~IdentityIMPL() {} ;
    id_type ident() { return id ; }
} ;

class Identity : virtual public IdentityIMPL {} ;
}

#endif
