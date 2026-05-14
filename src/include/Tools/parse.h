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
#ifndef PARSE_H
#define PARSE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <istream>
#include <string>

namespace Loci {
namespace parse {
    
    void kill_white_space(std::istream &s) ;
    
    bool is_name(std::istream &s) ;
    std::string get_name(std::istream &s) ;

    bool is_int(std::istream &s) ;
    long get_int(std::istream &s) ;

    bool is_real(std::istream &s) ;
    double get_real(std::istream &s) ;

    bool is_string(std::istream &s) ;
    std::string get_string(std::istream &s) ;

    bool is_token(std::istream &s, const std::string &token) ;
    bool get_token(std::istream &s, const std::string &token) ;
}
}

#endif
