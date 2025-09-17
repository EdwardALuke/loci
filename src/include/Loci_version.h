//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#ifndef LOCI_VERSION_H
#define LOCI_VERSION_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <string>

#define LOCI_VERSION_EQ(X,Y) (LOCI_VERSION_MAJOR == X && LOCI_VERSION_MINOR == Y)
#define LOCI_VERSION_LT(X,Y) (LOCI_VERSION_MAJOR < X || (LOCI_VERSION_MAJOR == X && LOCI_VERSION_MINOR < Y))
#define LOCI_VERSION_LE(X,Y) (LOCI_VERSION_MAJOR < X || (LOCI_VERSION_MAJOR == X && LOCI_VERSION_MINOR <= Y))
#define LOCI_VERSION_GT(X,Y) (LOCI_VERSION_MAJOR > X || (LOCI_VERSION_MAJOR == X && LOCI_VERSION_MINOR > Y))
#define LOCI_VERSION_GE(X,Y) (LOCI_VERSION_MAJOR > X || (LOCI_VERSION_MAJOR == X && LOCI_VERSION_MINOR >= Y))

namespace Loci {
  std::string version() ;
}

#endif
