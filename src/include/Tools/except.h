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
#ifndef EXCEPT_H
#define EXCEPT_H

#include <string>
#include <iostream>

namespace Loci {

    class BasicException {
      public:
        int code ;
        BasicException() {code = 0 ;}
        BasicException(int val) : code(val) {}
        virtual ~BasicException() {}
        virtual std::ostream &Print(std::ostream &s) const {
            s << "Loci exception, code = " << code << std::endl ;
            return s ;
        }
    } ;

    struct StringError : public BasicException {
      public:
        std::string msg ;
        StringError() {}
        StringError(std::string m) :msg(m) {code = 0 ;}
        virtual ~StringError() {}
        virtual std::ostream &Print(std::ostream &s) const {
            s << msg << std::endl ;
            return s ;
        }
    } ;

    struct StringWarning : public BasicException {
      public:
        std::string msg ;
        StringWarning() {}
        StringWarning(std::string m) :msg(m) {code = 0 ;}
        virtual std::ostream &Print(std::ostream &s) const {
            s << msg << std::endl ;
            return s ;
        }
    } ;

}

#endif

    
