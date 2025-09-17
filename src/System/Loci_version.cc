//#############################################################################
//#
//# Copyright 2008, 2025, Mississippi State University
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
#include <Loci_version.h>
#include <Config/version_details.h>
namespace Loci {


  std::string version() {
    std::string rn = std::string(LOCI_BRANCH)+"-"+std::string(LOCI_VERSION) ;
    rn += " Compiled On " ;
    rn += __DATE__ ;
    rn += " " ;
    rn += __TIME__ ;
    return rn ;
  }
}
